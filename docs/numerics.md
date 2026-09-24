## Monte Carlo

**Terminal-distribution sampling instead of path-stepping:**
For European options, the payoff depends only on the terminal price
S_T, not on the path taken to get there. GBM's terminal distribution
is known exactly in closed form, so each simulated "path" is a single
random draw jumping directly from today to expiry:

  S_T = S · exp[(r − q − ½σ²)T + σ√T · Z],  Z ~ N(0,1)

This is mathematically equivalent to stepping through many small time
increments and compounding the moves, just computed in one draw per
path instead of many. Full path-stepping (and the Euler vs. Milstein
discretization choice) becomes necessary later for genuinely
path-dependent payoffs — specifically Longstaff-Schwartz American
Monte Carlo, which needs the price at every intermediate exercise
date, not just at expiry. Deferred to that phase, not skipped
entirely.

**Reporting uncertainty honestly:**
Unlike closed-form or the tree, Monte Carlo produces a statistical
estimate, not an exact number. Every price is reported alongside its
standard error, and validation against Black-Scholes uses a proper
statistical test (price within 3 standard errors of the known true
value) rather than an arbitrary fixed tolerance — the correct way to
validate a random estimate.

**Variance reduction — antithetic variates:**
For every random draw Z, also compute the path using −Z. Z and −Z are
perfectly negatively correlated, so averaging the resulting payoff
pairs cancels out part of the sampling noise, for the same number of
underlying random draws. Measured: reduced standard error from 0.0663
to 0.0464 (~30% reduction) on a reference case (S=100, K=100, T=1,
r=5%, vol=20%, 50,000 paths).

**Variance reduction — control variates:**
Uses the simulated terminal price S_T itself as a control, since its
exact expectation E[S_T] = S·e^((r−q)T) is known analytically with no
simulation error. For each path, the option payoff and S_T are
correlated (both come from the same random draw). Estimate the
optimal control coefficient c* = Cov(payoff, S_T) / Var(S_T) from the
sample, then adjust each payoff by c*·(simulated S_T − known exact
E[S_T]) before averaging — correcting each path's estimate using the
known error on a related, exactly-computable quantity.

Measured on the same reference case: standard error reduced from
0.0663 (plain) to 0.0251 (control variate) — a **6.97x variance
reduction ratio**, meaning plain Monte Carlo would need roughly 7x
more simulated paths to reach the same precision control variates
achieve for free.

**Why MC needs discretization error consideration at all (general
case):**
Unlike closed-form, path-stepping MC (used later for American
exercise) simulates a discrete-time approximation of a continuous-time
process — every step introduces bias that shrinks as step size → 0,
distinct from the Monte Carlo sampling error (reported via standard
error above) that shrinks as paths → ∞. The terminal-distribution
approach used for European pricing in this phase has no discretization
error at all, since it draws directly from the exact terminal
distribution — only sampling error applies here.

## Automatic Adjoint Differentiation (AAD)

**Why hand-rolled instead of a library (XAD/CoDiPack):**
Understanding the tape mechanics well enough to build it is the
resume-defining piece of this project — using a library would mean not
being able to defend the mechanism under questioning.

**Tape design:**
`Tape` is a flat `vector<Node>`, where each `Node` stores up to two
parent indices and the local derivative ("partial") with respect to
each parent. `ADouble` wraps a `double` value together with its index
on the tape. Constructing an `ADouble` from a raw `double` creates a
new leaf node (used only for genuine inputs — S, sigma, r, T — that a
derivative is wanted for); every arithmetic operator (+, -, *, /,
ad_exp, ad_log, ad_sqrt) records a new node with the correct local
derivative, without the caller needing to think about calculus at all.

Reverse-mode was chosen over forward-mode specifically because it
computes derivatives with respect to *all* inputs (all 5 potential
Greeks) in a single backward pass, regardless of how many inputs
there are — forward mode would need one full forward pass per input,
which is a different but structurally similar cost problem to
bump-and-revalue.

**Backward pass mechanics:**
`Tape::backward(root)` seeds the root node's adjoint at 1.0, then
walks every node from last-recorded to first. Because a node's parents
are always recorded before the node itself, walking in strict reverse
order guarantees each node's adjoint is fully accumulated (from every
child that used it) before it's propagated further back — this is why
a single linear backward pass correctly implements the chain rule
across the whole graph, including nodes with multiple children (like
S being used in both the drift and the S_T computation).

**Discounting folded into the tape:**
The discount factor e^(−rT) is computed as `ad_exp(-(r*T))` and
multiplied into the payoff *on the tape*, rather than applied as a
separate plain-double multiply after backward() runs. This matters
specifically for rho: r affects the price two distinct ways (through
the drift inside the payoff, and through discounting itself), and
folding both onto the same tape means the backward pass automatically
combines both effects via the chain rule — no manual reasoning about
how the two contributions add together was required.

**Handling the payoff's kink:**
max(S_T − K, 0) isn't differentiable at the single point S_T = K. The
implementation checks which branch applies using the plain `.value`
(not the tape), and returns whichever branch's smooth derivative is
correct for that path. This is valid because that exact knife-edge
point has zero probability of being hit in a continuous simulation —
a standard, correct way to handle this in practice.

**Validated against:**
Delta, vega, theta, and rho from `monte_carlo_greeks` (single
simulation, AAD) all matched the closed-form analytical Greeks from
Phase 1 within tight tolerance, on the standard reference case
(S=100, K=100, T=1, r=5%, vol=20%).

**Pathwise differentiation — the underlying technique:**
Differentiating the pricing formula while holding the random draw Z
fixed for each path is a named technique (pathwise differentiation /
infinitesimal perturbation analysis). AAD is the mechanism used to
compute those pathwise derivatives efficiently — not a separate
technique from pathwise differentiation.

**AAD vs. bump-and-revalue — measured result (an honest one):**
Benchmarked computing 4 Greeks (delta, vega, theta, rho) for 100,000
paths: bump-and-revalue (1 base run + 4 bumped runs, sharing a fixed
random seed across all 5 runs so the comparison isolates each input's
true effect) took 94.3ms; AAD (1 run) took 150.5ms. **At this scale,
AAD was slower, not faster.**

This is a real, informative result, not a failed implementation: the
tape allocates and walks a full computation graph on every single
path, and that per-path overhead (heap allocation, backward-pass
traversal) is real cost a plain bump-and-revalue loop never pays.
AAD's actual structural advantage is that its per-path cost stays
FIXED regardless of how many Greeks are requested, while
bump-and-revalue's cost grows linearly with the Greek count (one full
extra simulation per Greek). At 4 Greeks, bump-and-revalue's added
simulation cost apparently hasn't yet crossed AAD's fixed tape
overhead. The crossover point would likely appear at a higher Greek
count, for a more complex instrument needing many risk factors, or
with a lower per-path tape overhead — a possible future optimization,
not attempted here given the phase's remaining scope.

**A real bug found in the process:**
The first version of this benchmark showed bump-and-revalue's delta
disagreeing with AAD's by 3.2 — far too large to be sampling noise.
Root cause: `monte_carlo_price` created a fresh, unseeded RNG on every
call, so the "base" and "bumped" runs in the bump-and-revalue
comparison used completely unrelated random paths, not the same
underlying randomness with one input nudged — which defeats the whole
point of a finite-difference comparison. Fixed by adding an optional
`seed` parameter (default 0 preserves existing random behavior
everywhere else; a nonzero seed forces reproducible draws). After the
fix, bump-delta (0.6363) and AAD-delta (0.6376) agree closely, as
expected.

**Gamma via AAD — deferred, then built (Phase 11):**
Gamma is a second-order derivative (∂²price/∂S², the derivative of
delta), which the first-order reverse-mode tape above cannot compute
directly. At the time of Phase 4, this was deliberately deferred:
computing it would need a genuinely more complex "tape of tapes"
structure to differentiate the backward pass itself, and given
remaining project scope at that point (Longstaff-Schwartz, PDE
solver, hedging simulation, calibration), that was judged lower-value
than covering new ground. Gamma remained available via the exact
analytical Black-Scholes formula (Phase 1) for any downstream use,
such as the delta-hedging simulation.

Phase 11 revisited this and built it: `dual.hpp`/`tape2.hpp`/
`adouble2.hpp` implement forward-over-reverse AD (a genuinely separate,
isolated tape, not a modification of the original one), which does
compute a Hessian-vector product exactly in one backward pass. Applied
to the smooth Black-Scholes formula, it recovers analytical gamma to
~1e-9. Applied to a Monte Carlo path's payoff — the more natural place
to want it, since that's what the rest of this project's AAD work
(monte_carlo_greeks) targets — it returns exactly zero, because a
single path's payoff is piecewise LINEAR in spot (either S_T-K or 0),
and the second derivative of a piecewise-linear function is zero
almost everywhere. This isn't a bug in the new tape (verified
independently against hand-computed toy derivatives first, the same
way the original tape was in Phase 4); it's the real, well-known
reason pathwise differentiation doesn't extend to gamma the way it
does to delta. A working Monte Carlo gamma estimator needs a
different technique (e.g. the likelihood-ratio/score-function method,
or a smoothed payoff) — left out of scope here, documented rather
than silently shipped as a gamma that's always zero. See
`docs/phase-log.md` Phase 11 and `tests/test_main.cpp`'s `[gamma]`
tests.

## Delta-hedging simulation

**Setup:**
Models the standard hedging scenario: the option is SOLD (written),
receiving the theoretical Black-Scholes premium upfront. The writer is
the party exposed to risk and needs to hedge; the buyer has no further
obligation after paying the premium. At each of `num_rebalances`
evenly-spaced points between today and expiry, delta is recomputed at
the current spot and remaining time, and the BTC hedge position is
traded to match it. Cash in the hedging account earns the risk-free
rate between rebalances.

**Why the PnL calculation needed three separate quantities, not one
running balance:**
The premium received, the hedging account's trading activity, and the
final obligation owed at expiry are three logically distinct things.
Mixing the premium directly into the same running `cash` variable as
the trading activity made the final comparison ambiguous — see bug
note below. Keeping them separate (hedging account tracked from zero,
premium grown at the risk-free rate, obligation computed independently
at expiry, then combined once at the very end) is what makes the
final PnL calculation correct and auditable.

**Validated against:**
- Single-path sanity check: simulation runs, produces a finite PnL
- Distribution check across 2000 simulated paths: mean hedging PnL
  ≈ 0 (within a small fraction of the option's own price), confirming
  the hedge breaks even on average as theory predicts
- Rebalancing frequency comparison: hedging error (PnL standard
  deviation) shrinks as rebalancing gets more frequent — monthly (12
  rebalances) stdev=1.95, weekly (52) stdev=0.94, daily (252)
  stdev=0.44. Roughly a 4.5x reduction from monthly to daily.
  Qualitatively consistent with the theoretical result that
  discretization-driven hedging error scales with √(dt) — each ~4-5x
  increase in rebalancing frequency roughly halved the standard
  deviation, matching a square-root relationship rather than a linear
  one.

**Bug found and fixed:**
The first implementation initialized the hedging account's `cash`
directly with the theoretical premium, then mixed all subsequent
trading activity into that same variable, then subtracted the premium
(grown at the risk-free rate) again at the end to compute final PnL.
This double-counted the premium's role — once implicitly via the
initial value, once explicitly in the final subtraction. Symptom: a
single-path test showed PnL of -11.12 against a theoretical price of
10.45 (suspiciously close to -theoretical_price, the actual tell); a
2000-path distribution test then confirmed this wasn't noise — the
mean PnL was consistently ≈ -11 across all paths, a systematic bias
rather than random variation. Fixed by tracking the hedging account
from zero, and combining premium + hedging account + obligation as
three separate terms only once, at the very end. Post-fix, mean PnL
across 2000 paths landed at approximately 0.02, as expected.

**Why this connects back to Black-Scholes' assumptions:**
The original Black-Scholes derivation (Phase 1) assumes continuous
hedging — rebalancing infinitely often. The nonzero standard deviation
measured here, even at 252 rebalances (daily, roughly the most
frequent a real desk would practically do), is the direct, measured
gap between that theoretical idealization and what's achievable with
discrete rebalancing. This is the same PDE-derivation argument from
Phase 1 made empirically concrete: the tighter the rebalancing
interval, the closer the realized outcome tracks the theoretical
price, exactly as the continuous-hedging assumption predicts.

## Longstaff-Schwartz (American Monte Carlo)

**The problem this solves:**
European Monte Carlo only needs the terminal price. American options
require, at every point along every path, deciding "exercise now, or
keep holding?" — which needs an estimate of the expected value of
continuing to hold. That estimate seems to require knowing the future,
which appears circular mid-simulation.

**The resolution — regression on already-simulated paths:**
Simulate every path forward to expiry FIRST (so the future is already
known for every path), then walk backward through time. At each step,
fit a regression curve — continuation value as a function of current
price — using the discounted realized future cash flow of each
simulated path as training data. Compare that regression-estimated
continuation value against the immediate exercise value at each node,
and take whichever is larger, exactly the same decision rule as the
binomial tree (Phase 2), just estimated statistically instead of
computed exactly at a discrete node.

**Why regression only over in-the-money paths:**
Paths where exercise isn't a candidate (out-of-the-money) contribute
no useful information to the exercise/continue decision, and including
them dilutes the regression fit specifically in the price region where
the decision actually matters. This is a deliberate detail from the
original Longstaff-Schwartz (2001) formulation, not an arbitrary
simplification.

**Why quadratic (not linear or higher-order) basis functions:**
A quadratic captures curvature in the continuation-value relationship
(visible in the earlier diagram: continuation value isn't linear in
spot) with a small, cheap-to-fit 3-coefficient regression. Higher-order
polynomials risk overfitting with the amount of data available at each
time step, especially near the tails where fewer paths are in the
money.

**Why full path-stepping is required here (unlike European MC):**
European pricing only needs the terminal distribution (Phase 3), which
can be sampled directly in one draw. American exercise requires the
price at every intermediate exercise date, since the decision happens
at each of those dates — this is exactly the case flagged as deferred
in Phase 3's numerics notes.

**Regression solver:**
A hand-implemented least-squares quadratic fit via the normal
equations, solved with Gaussian elimination and partial pivoting
(row-swapping to avoid dividing by a near-zero pivot, a standard
numerical-stability safeguard). Verified independently, before use in
the full algorithm, by confirming it exactly recovers known
coefficients from zero-noise synthetic data — isolating correctness
of the linear algebra from any question about how well it fits noisy
simulated data.

**Validated against:**
- Regression solver alone: exact recovery of known coefficients from
  noise-free synthetic data
- Full LSM price converges to the (already-trusted) binomial tree
  price for an American put with a nonzero dividend yield
- A second, meaningfully different case (in-the-money call, different
  dividend yield, volatility, and time-to-expiry) also converges —
  confirming the first result wasn't specific to one parameter set

  ## Finite Difference (Crank-Nicolson)

**The core idea:**
Solve the Black-Scholes PDE numerically on a 2D grid of (price, time)
points, rather than symbolically (closed-form) or by random sampling
(Monte Carlo). Known values at expiry (the payoff) and at extreme
price boundaries seed the grid; the algorithm walks backward in time,
solving for each earlier column from the columns already known —
structurally the same backward-walking pattern as the tree and LSM,
applied to a PDE discretization instead.

**Why Crank-Nicolson over explicit or implicit:**
Explicit schemes are simple (each new value depends only on
already-known values, no system to solve) but only conditionally
stable — the time step must be small relative to the square of the
price-grid spacing, or the scheme diverges. Implicit schemes are
unconditionally stable but only first-order accurate in time.
Crank-Nicolson averages the explicit and implicit formulations,
achieving both unconditional stability and second-order accuracy.

**Stability demonstrated directly, not just asserted:**
Ran both an explicit scheme and Crank-Nicolson on the identical grid
(200 price steps, 20 time steps — deliberately coarse in time relative
to the fine price grid, the specific combination that breaks explicit
stability). Explicit scheme result: -4.06×10^22 (a nonsensical,
diverged value for an option price). Crank-Nicolson result: 10.4548,
against a true Black-Scholes price of 10.4506 — still essentially
correct on the same grid where the explicit scheme collapsed entirely.

**Solving each time step — tridiagonal system via Thomas algorithm:**
Discretizing the PDE's derivative terms (∂V/∂S, ∂²V/∂S², via central
differences) and averaging explicit/implicit versions produces, per
time step, a linear system where each unknown only depends on its two
immediate neighbors — a tridiagonal system. Solved via the Thomas
algorithm (forward elimination + back-substitution), which is O(n)
rather than the O(n³) general Gaussian elimination would cost —
significant given this runs once per time step across the whole grid.

**Regression solver — err, tridiagonal solver — verified independently
first:**
Same discipline as the LSM regression solver (Phase 6): before
trusting the Thomas algorithm inside the full PDE solver, verified it
against a small, hand-solvable tridiagonal system with a known exact
answer, isolating correctness of the linear algebra from correctness
of the PDE discretization.

**Boundary conditions:**
At S=0 and S=S_max (chosen as 3x the strike, far enough that the known
asymptotic behavior is essentially exact), the value is set directly
from known limiting behavior rather than the interior stencil, which
needs neighbors on both sides that don't exist at the grid's edges.

**Validated against:**
- Tridiagonal solver alone: exact match against a hand-solvable known
  system
- Full Crank-Nicolson price matches Black-Scholes closed-form for both
  European calls and puts
- Stability: explicit scheme diverges catastrophically on a grid where
  Crank-Nicolson remains accurate — the concrete proof behind the
  "unconditionally stable" claim, not just an assertion

## Implied Volatility & Real-Market Calibration

**The problem being solved:**
Every pricing method so far takes volatility as an INPUT. Real markets
quote option PRICES, not volatilities — volatility has to be backed
out from the price the market is actually willing to pay. This is the
inverse problem: given (S, K, T, r, market_price), find the σ that
makes black_scholes_price reproduce that price.

**Why Newton-Raphson, and why it needs vega:**
Black-Scholes has no closed-form inverse with respect to σ, so the
implied vol has to be found numerically. Newton-Raphson iteratively
refines a guess using the local slope of price with respect to σ —
which is exactly vega, already built and validated in Phase 1. Each
iteration: price the option at the current vol guess, measure how far
off the price is, and use vega to compute a smarter next guess. This
converges very fast when it works — the real-data test below
converged in 4 iterations.

**Handling near-zero vega (a real edge case, not hypothetical):**
Deep ITM/OTM options and options very close to expiry have vega near
zero (this was already discovered and tested explicitly in Phase 1's
edge case suite). Dividing by a near-zero vega in the Newton-Raphson
update would either blow up or produce a meaningless jump. The solver
detects this and returns cleanly with converged=false rather than
crashing or returning garbage — not every real market quote has a
well-defined, numerically stable implied vol, and the solver is
honest about that rather than forcing an answer.

**Validated against:**
- Round-trip test: given a known (price, vol) pair from Black-Scholes
  itself, the solver recovers the exact vol that produced that price
- Edge case: a deep OTM option with a tiny price does not crash or
  produce nan/inf, even when convergence isn't guaranteed
- REAL DATA: pulled a live BTC put (BTC-30OCT26-90000-P) from
  Deribit's public API (spot=$76,871.89, strike=$90,000, mark
  price=0.176 BTC = $13,529.50, ~43 days to expiry). The solver
  independently computed 34.13% implied vol; Deribit's own reported
  mark_iv was 34.36% — a gap of 0.23 percentage points, converging in
  4 iterations.

**Why the small gap to Deribit's own mark_iv is expected, not a bug:**
- BTC options are typically coin-margined (settled and margined in
  BTC, not USD) — Deribit's own pricing may include a quanto-style
  adjustment that vanilla Black-Scholes does not model
- The curl snapshot and Deribit's mark_iv were computed at slightly
  different instants; BTC's price (and therefore implied vol) moves
  continuously
- Deribit's mark_iv is itself likely smoothed/fitted across their own
  internal vol curve, not simply solved from that single mark price
  the same way this solver does it

**Data source:**
Deribit's public REST API (no authentication required), verified
directly (via curl) to be accessible from the US before committing to
it as the project's data source, per the standing rule of confirming
accessibility, cost, and data depth before adopting any external data
source. Prices are quoted in BTC, not USD, and must be converted
(price_usd = price_btc * underlying_price_usd) before use.

## Volatility Surface (real market data)

**What was built:**
`build_vol_surface_from_snapshot()` parses a full Deribit BTC options
chain (JSON, via nlohmann/json), decodes Deribit's instrument naming
convention (e.g. BTC-30OCT26-90000-P) into strike/expiry/type, converts
BTC-denominated prices to USD, and runs the implied vol solver
(independently, not just reading Deribit's own mark_iv) across every
contract in the chain.

**Reproducibility — snapshot-based, not live:**
The chain is pulled once via curl and saved to data/btc_chain_snapshot.json,
rather than the test suite hitting Deribit's live API on every run. A
live API dependency would make tests flaky (network failures, and
results that silently change over time as the real market moves) —
the same reproducibility principle applied throughout this project.
today_year/month/day are passed explicitly rather than read from the
system clock, so time-to-expiry calculations are pinned to when the
snapshot was taken, not whenever the code happens to run later.

**Full-market validation — not just one contract:**
Ran across the entire live BTC options chain: 904 total contracts
parsed, 860 (95.1%) converged. Across the converged contracts, this
engine's independently-solved implied vol averaged just 0.019 (1.9
percentage points) away from Deribit's own reported mark_iv — this
is the same result observed on the single hand-picked contract in the
previous section, now confirmed to hold at market scale, not a
coincidence of one lucky example.

**The volatility smile, visualized from real data:**
Plotted implied vol against strike for the most heavily-quoted expiry
in the snapshot (~99 days out, 118 contracts). The resulting curve is
a textbook volatility smile: implied vol bottoms out near-the-money
(~36% around $84,000-86,000 strikes, close to spot), and rises on
both wings — over 100% for deep OTM puts near $20,000, and up to ~68%
for deep OTM calls near $230,000. This is the real market pricing in
tail risk that Black-Scholes' constant-volatility assumption cannot
represent — the exact limitation flagged back in Phase 1's assumptions
list, now demonstrated empirically with live data rather than just
stated as a known theoretical gap.

**Non-convergence — investigated, not just counted:**
The 44 non-converged contracts (4.9%) were not a mystery left
unexamined. Inspecting them directly showed every one falls into the
already-known vega-collapse failure mode from Phase 1/8's edge-case
work: extreme moneyness (strike/spot ratios from 0.39 to 4.18) and/or
very short time-to-expiry (many under 8 days, some under 1 day).
Their solved_iv values confirm the specific mechanism: many are stuck
at exactly 0.5 (the solver's starting guess, meaning the near-zero-
vega guard triggered on the very first iteration) or exactly 0.01
(the clamp floor, meaning one Newton step overshot before hitting the
guard). This is the predicted pattern from the solver's design,
confirmed directly against real failure data rather than assumed.

## Heston calibration: single-expiry vs. multi-expiry (Phase 12)

**The identifiability problem, made concrete:**
`calibrate_heston()` (`heston_calibration.hpp`) was, from the start, generic
over which points it's given -- it prices each point at its own
`time_to_expiry`, with nothing in the code assuming a single expiry. What
changed in this phase is not the calibrator, but what's fed into it, and
an honest accounting of what that choice buys and costs.

Fit Heston to one expiry's smile, and kappa (mean-reversion speed) and
theta (long-run variance) are not separately identifiable: many (kappa,
theta) pairs price that one maturity almost identically, since what a
single-T price actually constrains is closer to an integrated combination
of the two. A synthetic test (`tests/test_main.cpp`, "Multi-expiry Heston
calibration resolves the single-expiry kappa/theta non-identifiability")
makes this concrete rather than asserted: starting from the same
deliberately-bad initial guess (kappa=6.0 against a true kappa=2.0),
calibrating against a single 90-day synthetic expiry recovers kappa over
200% off the true value -- while fitting that expiry's prices to within
0.01% RMSE. That's the identifiability failure in one number: an
excellent fit, a wrong kappa. Calibrating the same optimizer against the
same three synthetic expiries (30d, 90d, 270d) pooled, from the same bad
guess, recovers all five parameters to within floating-point noise of the
truth. Pooling expiries works because variance mean-reverts at a specific
rate: contracts at different maturities sample different stages of that
reversion, and fitting all of them at once with one shared parameter set
pins the rate down.

**On real data, the honest result is messier:**
`tools/calibrate_heston_multi.cpp` runs the same joint fit against the
real Deribit snapshot, pooling the ~43-day and ~99-day expiries (94 and 66
contracts respectively, same moneyness/liquidity filter as
`tools/calibrate_heston.cpp`). Fit quality: 1.35 percentage points RMSE in
vol space -- comfortably worse than the single-expiry tool's 0.43pp on its
own slice, though still better than the flat-vol baseline (headline result
above, ~1.9pp mean gap across the whole chain). Pooling a third or fourth
expiry (190d, 281d) was tried during development and makes this worse,
not better: 3 expiries lands around 2.2pp, 4 expiries around 2.9-3.1pp,
degrading as more maturities are added.

This is not a bug in the calibrator or the pooling logic -- it's
`calibrate_heston` correctly reporting that a single (kappa, theta, xi,
rho) parameter set genuinely cannot reproduce multiple real maturities'
smiles as well as a fresh fit reproduces any one of them. Real markets'
volatility dynamics have term structure that one-factor Heston, run with
fixed parameters, doesn't fully capture -- a well-known limitation in the
literature, not specific to this implementation. Properly closing that
gap needs a term-structure extension (piecewise-constant parameters, or a
multi-factor model like Bergomi), which is out of scope here, in the same
spirit as the other deliberate scope decisions in this document.

**Why `tools/calibrate_heston_multi.cpp` isn't wired into the live demo:**
The single-expiry tool feeds the demo's Heston card because it answers a
clean question (does Heston fit today's smile well?) with a clean, strong
number (0.43pp). The multi-expiry tool answers a more nuanced question
(can one Heston parameter set explain several maturities at once?) with
an honest but weaker number, and that nuance doesn't compress well into a
single demo stat without the surrounding explanation this section
provides. It's kept as a standalone research tool and test, documented
here rather than headlined.

## Exotic options: Asian and Barrier (Phase 13)

**Why these needed full path simulation, and the vanilla European pricer
didn't:** the Monte Carlo section above explains that vanilla European
options only need the terminal price `S_T`, so `monte_carlo_price()`
deliberately jumps straight from today to expiry in one draw per path.
Asian and Barrier payoffs depend on the *whole* path -- an Asian option's
payoff is a function of the average price along the way, and a Barrier
option's payoff depends on whether the path ever touches a level before
expiry -- so `exotic_options.cpp` steps through `num_steps` intermediate
dates per path instead, reusing the exact Euler-exact GBM stepping
`longstaff_schwartz.cpp` already uses for American Monte Carlo
(`S[step] = S[step-1] * exp((r-q-0.5*sigma^2)*dt + sigma*sqrt(dt)*Z)`).
This is the same "deferred to that phase" path-stepping the Monte Carlo
section above flags -- this is that phase, now for two more payoff types.

**Validating a pricer with no closed form, by building one that has one:**
There's no closed-form price for an arithmetic-average Asian option --
the sum of lognormal random variables isn't itself lognormal, so the
averaging breaks the algebra that makes Black-Scholes solvable. But the
*geometric* average of those same lognormal prices, sampled at equally
spaced dates, is itself exactly lognormal, which gives the geometric
Asian option a closed form (Kemna & Vorst, 1990): price it like a vanilla
Black-Scholes option, but with the drift and volatility replaced by

```
sigma_hat^2 = sigma^2 * (n+1)(2n+1) / (6n^2)
mu_hat      = (r - q - 0.5*sigma^2)*(n+1)/(2n) + 0.5*sigma_hat^2
```

for `n` monitoring dates, discounted at the real rate `r`. This doesn't
directly validate the arithmetic-average price anyone actually wants --
but it validates the *path-simulation machinery* that both the arithmetic
and geometric pricers share (`geometric_asian_price_mc()` reuses the same
`simulate_paths()` helper, only the averaging formula differs): if the
simulator reproduces a known closed form exactly, the same simulator's
arithmetic-average output can be trusted too. Measured:
`geometric_asian_price_closed_form` = 5.6411 vs. `geometric_asian_price_mc`
(300,000 paths) = 5.6680 ± 0.0175 -- within 1.5 standard errors, well
inside the 4-standard-error bound the test asserts.

**Two more absolute-correctness checks, both provable without a Monte
Carlo comparison at all:**
- AM-GM inequality: pathwise, the arithmetic average of a set of positive
  numbers is always ≥ their geometric average, and `max(x-K, 0)` is
  nondecreasing in `x` -- so the arithmetic Asian call price is provably
  ≥ the geometric Asian call price, in expectation as well as pathwise.
  Tested directly (same seed, same paths under the hood, only the
  averaging differs) rather than just asserted.
- `num_steps=1` collapses an Asian option to a vanilla European option
  exactly -- averaging one price is that price. Tested by comparing
  `asian_option_price_mc` with one monitoring date directly against
  `monte_carlo_price`, same seed, both stepping the same dynamics.

**Barrier options: method of images, then in/out parity for the rest:**
A down-and-out call, with the barrier `H` at or below both spot and
strike, has a closed form via the reflection principle (method of
images):

```
C_do(S,K,T) = C_BS(S,K,T) - (H/S)^(2*nu) * C_BS(H^2/S, K, T)
nu = (r-q)/sigma^2 - 1/2
```

This single-reflection formula only holds under that `H <= min(S,K)`
constraint -- above it, the vanilla payoff isn't zero at the barrier and
one reflection no longer cancels it exactly, so the general
Reiner-Rubinstein (1991) formula (16 cases across all direction/strike
configurations) would be needed to cover barrier levels above the
strike. That full formula was deliberately not implemented: it's easy to
get subtly wrong in exactly the case-selection logic that's hardest to
unit-test, and there was no independent reference implementation on hand
to check it against here. `down_and_out_call_price_closed_form()` is
scoped to the one case that's implementable with confidence.

The other three directions (down-and-in, up-and-out, up-and-in) are
validated without needing their own formulas at all, using a
model-independent identity: every path either touches the barrier or it
doesn't, so a knock-out and its same-direction knock-in partition every
path exactly, and

```
(same-direction knock-out price) + (same-direction knock-in price) == vanilla price
```

holds regardless of the underlying's dynamics -- it's a statement about
set partitioning, not about GBM specifically. Measured (200,000 paths,
same seed for the out/in pair so they share the same simulated paths and
only differ in which side of the partition each path falls on):
down-and-out + down-and-in = 10.4265 vs. vanilla = 10.4506; up-and-out +
up-and-in = 10.4265 vs. vanilla = 10.4506 -- both well within the
combined standard error the tests assert against.

**Why these aren't wired into the live WASM demo:** consistent with the
Phase 12 multi-expiry Heston tool above, Asian and Barrier pricing ship
as engine code, tests, and this documentation, but the demo page itself
is left unchanged. The demo's existing cards (Black-Scholes, binomial/
trinomial trees, Monte Carlo, Heston) each answer one clean question with
one headline number; Asian and Barrier pricing would need their own
inputs (monitoring frequency, barrier level and direction) and don't
compress into the existing card layout without a UI redesign that's out
of scope for this phase.
