# Phase Log

One entry per phase. Written the day the phase is finished, while it's
fresh — not backfilled later. Each entry: what was built, what it was
validated against, and a short "defend this" paragraph — the answer
you'd give if an interviewer picked this phase apart.

---

## Phase 0 — Scaffold & tooling (done)

**What was built:**
- CMake + Catch2 project skeleton (`deriv-engine`), Clang toolchain,
  Ninja generator
- `.gitignore`, initial git repo, pushed to
  github.com/roberto0607/deriv-engine
- GitHub Actions CI: configure → build → test on every push to `main`,
  running on a clean Ubuntu container
- `docs/design.md` and `docs/numerics.md` skeletons committed alongside
  the code, not bolted on after

**Validated against:**
- Local build succeeds (`cmake --build`, exit code 0, all 4 targets:
  `deriv_engine`, `Catch2`, `Catch2WithMain`, `tests`)
- CI run green on GitHub Actions (fresh container, no local-machine
  assumptions) — run completed in 47s

**Defend this:**
The project builds from a clean checkout with zero manual setup beyond
CMake and a compiler — no hand-installed dependencies, no
machine-specific paths. Catch2 is pulled via CMake FetchContent so
anyone cloning the repo gets a working build immediately. CI proves
this isn't just "works on my machine": it runs the same configure/
build/test sequence on a fresh container on every push.

---

## Phase 1 — Black-Scholes closed-form (done)

**What was built:**
- `MarketData` struct (spot, risk-free rate, dividend yield, volatility)
  and `EuropeanOption` struct (strike, time-to-expiry, call/put type)
- `black_scholes_price()` — closed-form pricing for European calls and
  puts
- `black_scholes_greeks()` — analytical delta, gamma, vega, theta, rho
- Internal `norm_cdf()` and `norm_pdf()` helpers
- Guard against division-by-zero at T<=0 in both pricing and Greeks
  (see bug note below)
- `year_fraction()` — day-count convention helper (ACT/365, ACT/360)

**Validated against:**
- Standard textbook reference case (S=100, K=100, T=1yr, r=5%, vol=20%):
  call = 10.4506, put = 5.5735 — matched within 0.1% tolerance
- Same reference case Greeks (call delta=0.6368, gamma=0.0188,
  vega=37.52, theta=-6.414, rho=53.23; put delta=-0.3632, theta=-1.658,
  rho=-41.89) — all matched within tolerance
- Put-call parity (C - P = S*e^(-qT) - K*e^(-rT)) holds for an
  out-of-the-money case with nonzero dividend yield, not just the
  original at-the-money reference case
- Zero-time-to-expiry edge case for both pricing and Greeks — see bug
  note
- Deep ITM call: price and delta converge to the forward-price/delta=1
  limit; deep OTM call: price and delta converge to ~0
- Negative interest rate: prices remain finite and positive, put-call
  parity still holds
- Day-count: ACT/365 vs ACT/360 produce correctly differing year
  fractions for the same calendar days; 0 days gives exactly 0

**Bug found and fixed:**
At T=0 (or very close to it), the original d1/d2 formula divides by
`sigma * sqrt(T)`, producing `nan`. Caught by an explicit test
(`REQUIRE(std::isfinite(price))`) rather than by inspection. Fixed by
adding an early-return guard: at T<=1e-8, price returns intrinsic value
directly (max(S-K,0) for a call), and Greeks return the boundary values
(delta = 0/1/-1 depending on moneyness, all other Greeks = 0) — this is
the mathematically correct limit as T->0, not just a defensive hack.

**Scope decision — discrete dividends deferred:**
Project scope is crypto-focused (BTC options), where continuous
dividend yield with q=0 is the mathematically correct assumption,
since BTC pays no dividends. Discrete dividend handling (adjusting
spot by the present value of scheduled cash payments) would only
matter if scope expanded to equity options — deliberately deferred,
not an oversight.

**Defend this:**
Can explain the formula at both levels: the intuition (a call's price
is the expected value of the asset in the exercise scenarios, minus
the discounted strike weighted by probability of exercise) and the
precise mechanics (N(d2) is the risk-neutral probability of finishing
in-the-money; d1 differs from d2 by one sigma*sqrt(T) term because it
weights the expected asset price rather than pure exercise probability,
correcting for the log-normal distribution's mean being pulled above
its median by volatility). Can explain each Greek in plain terms:
delta as price-sensitivity to the underlying, gamma as delta's own
rate of change (and why it flattens to 0 for a deep-ITM option once
delta saturates near 1), vega as sensitivity to volatility specifically
(not to be confused with theta), theta as time decay, rho as rate
sensitivity. Can derive the Black-Scholes PDE from the replicating-
portfolio / delta-hedging argument via Ito's Lemma, and explain why
the stock's drift term cancels out entirely (risk-neutral pricing).
Can explain why higher volatility always increases option value
regardless of call/put — capped downside, open-ended upside. Can
explain why day-count convention choice changes priced output and why
ACT/360 produces a larger year-fraction than ACT/365 for the same
calendar days. Phase considered closed: all identified gaps are either
resolved and tested, or explicitly and defensibly scoped out.

---

## Phase 2 — Tree pricers (done)

**What was built:**
- `AmericanOption` struct, added alongside `EuropeanOption` in
  instrument.hpp — separate types rather than a shared struct with an
  exercise-style flag, so the compiler enforces correct usage
- `binomial_tree_price()` — CRR binomial tree, overloaded for both
  `EuropeanOption` and `AmericanOption`, sharing a single internal
  `run_tree()` implementation controlled by an early_exercise flag
- Risk-neutral up/down move factors (u, d=1/u) and risk-neutral
  probability p, including the dividend yield q in the probability
  calculation

**Validated against:**
- Convergence to Black-Scholes closed-form for a European call: within
  2% at 10 steps, within 0.1% at 500 steps
- Convergence to Black-Scholes closed-form for a European put: within
  0.1% at 500 steps
- American put price ≥ European put price under a nonzero dividend
  yield (deliberately not tested on a call, since American calls on
  non-dividend-paying assets are never optimal to exercise early — the
  test would pass trivially without proving anything)

**Defend this:**
Can explain why d=1/u specifically makes the tree recombine — an
up-then-down path lands on the same price node as down-then-up, so
there are only N+1 distinct nodes at step N instead of 2^N distinct
paths, which is what makes the tree computationally tractable at all.
Can explain why trees can price American options and closed-form
cannot: at every node walking backward from expiry, the tree compares
continuation value against immediate exercise value and takes the
max, a node-by-node decision closed-form has no mechanism for. Can
explain why early exercise only matters for puts (or dividend-paying
calls), not non-dividend calls. Trinomial trees deliberately deferred
— binomial with sufficient steps meets this project's accuracy needs,
and remaining time is prioritized toward Monte Carlo and AAD, which
are higher-value differentiators for the resume goal.


---

## Phase 3 — Monte Carlo core + variance reduction (done)

**What was built:**
- `monte_carlo_price()` — European option pricing via direct sampling
  of GBM's terminal distribution (no path-stepping needed for
  European payoffs), with `MonteCarloResult` reporting both price and
  standard error
- Antithetic variates (default-on flag) — pairs each random draw Z
  with its negation −Z to reduce variance for the same draw count
- `monte_carlo_price_control_variate()` — uses the simulated terminal
  price S_T as a control variate, with its known exact expectation
  E[S_T] = S·e^((r−q)T), estimating the optimal control coefficient
  from sample covariance/variance

**Validated against:**
- Convergence to Black-Scholes closed-form for both European calls
  and puts, using the statistically correct test (price within 3
  standard errors of the known true value, not an arbitrary fixed
  tolerance) — 100,000 paths per test
- Antithetic variates confirmed to reduce standard error vs. plain MC
  on identical inputs
- Control variate confirmed to reduce standard error vs. plain MC,
  and its adjusted price still validated against Black-Scholes within
  statistical tolerance
- Measured variance reduction on a reference case (S=100, K=100, T=1,
  r=5%, vol=20%, 50,000 paths): plain SE=0.0663, antithetic
  SE=0.0464 (~30% reduction), control variate SE=0.0251 — a **6.97x
  variance reduction ratio** over plain Monte Carlo

**Defend this:**
Can explain why European options can be priced by sampling GBM's
terminal distribution directly, in one random draw per path, rather
than stepping through many small time increments — the payoff only
depends on S_T, and GBM's terminal distribution is known exactly in
closed form. Can explain why this changes for path-dependent payoffs
(deferred to the Longstaff-Schwartz phase, which needs prices at
intermediate exercise dates, not just expiry). Can explain why a
Monte Carlo price must be reported with its standard error, and why
validating it against Black-Scholes requires a statistical tolerance
(within N standard errors) rather than an arbitrary fixed epsilon —
this is a fundamentally different kind of correctness check than the
deterministic tree and closed-form tests. Can explain both variance
reduction techniques mechanically: antithetic variates exploit the
perfect negative correlation between Z and −Z; control variates
exploit the correlation between the option payoff and a
different, exactly-computable quantity (S_T) to correct each path's
estimate using its known error on that quantity. Can state and defend
the measured 6.97x variance reduction ratio as a concrete, reproducible
result, not just a qualitative claim.

---

## Phase 4 — AAD Greeks (done)

**What was built:**
- `Tape` — a flat reverse-mode computation-graph recorder (`push_leaf`,
  `push_unary`, `push_binary`, `backward`), one tape per thread via
  `thread_local get_tape()`
- `ADouble` — a number type wrapping a `double` value and a tape index,
  with operator overloads (+, -, *, /, unary -) and named functions
  (`ad_exp`, `ad_log`, `ad_sqrt`) that transparently record every
  operation onto the tape
- `monte_carlo_greeks()` — runs the GBM simulation using `ADouble` for
  S, sigma, r, and T, with discounting folded directly into the tape,
  producing price + delta + vega + theta + rho from a single
  simulation pass
- Optional `seed` parameter added to `monte_carlo_price()`, enabling
  reproducible randomness for finite-difference comparisons

**Validated against:**
- Tape mechanics verified against a hand-computed toy example
  (y = a·b + c) before any real pricing code touched it
- Delta and vega from AAD matched analytical Black-Scholes Greeks on
  the standard reference case
- Theta and rho, added in a second pass with discounting folded onto
  the tape, also matched analytical Greeks
- AAD vs. bump-and-revalue benchmarked directly: 4 Greeks, 100,000
  paths — bump-and-revalue 94.3ms, AAD 150.5ms. AAD was slower at this
  scale, not faster — a real, documented finding, not a discarded
  result (see numerics.md for the full reasoning on why, and where the
  crossover point would actually be)

**Bug found and fixed:**
The first version of the AAD-vs-bump-and-revalue benchmark showed the
two methods' delta disagreeing by 3.2 — far too large to be sampling
noise. Root cause: `monte_carlo_price` used a fresh unseeded RNG per
call, so bump-and-revalue's "base" and "bumped" runs used unrelated
random paths instead of shared randomness with one input nudged,
defeating the entire premise of a finite-difference comparison. Fixed
by adding an optional `seed` parameter (default preserves existing
behavior). Post-fix, bump-delta and AAD-delta agree closely (0.6363
vs. 0.6376).

**Scope decision — gamma deferred:**
Gamma requires second-order differentiation (a tape-of-tapes
structure to differentiate the backward pass itself), a genuinely
larger sub-project than first-order AAD. Given remaining phases
(Longstaff-Schwartz, PDE, hedging simulation, calibration) cover new
technical ground, gamma-via-AAD was judged lower-value than that
remaining breadth. Gamma is still available via the exact analytical
formula from Phase 1 wherever needed downstream (e.g. the hedging
simulation).

**Defend this:**
Can explain reverse-mode AAD mechanically: why walking the tape from
last-recorded node to first guarantees each node's adjoint is fully
accumulated before being propagated further back, and why this
correctly implements the chain rule across a graph where a value (like
S) has multiple children. Can explain why reverse mode was chosen over
forward mode: one backward pass yields derivatives w.r.t. every input,
where forward mode would need one full pass per input. Can explain
"pathwise differentiation" as the named technique AAD implements here
— differentiating while holding the random draw Z fixed per path. Can
explain the payoff kink handling (branching on the plain value, valid
since the kink has zero probability of being hit exactly). Can explain
why folding discounting onto the tape was necessary for rho
specifically, since r affects price through two channels that the
chain rule needs to combine automatically. Can state and defend the
counterintuitive benchmark result — AAD was slower at this scale — and
explain precisely why: fixed per-path tape overhead vs.
bump-and-revalue's cost scaling with Greek count, and where the
crossover would actually occur. Can explain the seed-sharing bug and
why finite-difference comparisons require shared randomness between
base and bumped runs. Can explain, without hedging, why gamma was
deliberately not attempted via AAD in this phase.

---

## Phase 5 — Delta-hedging PnL simulation (done)

**What was built:**
- `simulate_delta_hedge()` — simulates one BTC price path, recomputing
  delta and rebalancing the hedge at each of `num_rebalances`
  evenly-spaced points, tracking the hedging account's cash flow
  (trades + risk-free interest) separately from the premium and the
  final obligation
- Reused Phase 1 analytical Greeks (recomputed at each rebalancing
  point, at the current spot and remaining time) and Phase 3's GBM
  step logic

**Validated against:**
- Single-path run: finite PnL, correct price-path length
- 2000-path distribution: mean hedging PnL ≈ 0 (small relative to the
  option's own price), confirming the hedge breaks even on average
- Rebalancing frequency comparison: hedging error standard deviation
  shrinks from 1.95 (monthly) to 0.94 (weekly) to 0.44 (daily) —
  roughly a 4.5x reduction from monthly to daily, consistent with the
  theoretical √(dt) scaling of discretization-driven hedging error

**Bug found and fixed:**
The premium was double-counted in the PnL formula — folded into the
initial `cash` value and then subtracted again at the end. A
single-path test first showed a suspiciously large PnL (-11.12,
notably close to -theoretical_price); the 2000-path distribution test
confirmed this was systematic (mean ≈ -11 across all paths), not
random noise, which is what made it identifiable as a real bug rather
than expected variance. Fixed by tracking the hedging account
separately from zero and combining premium, hedging account, and
obligation as three distinct terms only once, at the end.

**Defend this:**
Can explain the hedging setup: the option writer receives the premium
upfront and is the party who needs to hedge, since the buyer's
obligation ends at payment. Can explain why cash earns the risk-free
rate between rebalances and why omitting this would make the
simulation wrong. Can explain why a single-path PnL check is
insufficient to validate a random simulation, and why a distribution
check across many paths was necessary to confirm both correctness (the
premium bug) and expected behavior (mean near zero). Can explain the
rebalancing-frequency result concretely: more frequent rebalancing
reduces hedging error because it more closely approximates the
continuous hedging assumed in the original Black-Scholes derivation,
and can connect the roughly-halving pattern at each ~4-5x frequency
increase to the theoretical √(dt) scaling rather than treating it as
just "smaller error, more rebalances." Can walk through the
double-counted-premium bug end to end: the symptom, why the
suspicious magnitude was the first clue, why a distribution test (not
just a single path) was needed to distinguish "systematic bug" from
"unlucky random path," and the actual fix.

---

## Phase 6 — Longstaff-Schwartz American Monte Carlo (done)

**What was built:**
- `fit_quadratic()` — hand-implemented least-squares quadratic
  regression via normal equations and Gaussian elimination with
  partial pivoting
- `longstaff_schwartz_price()` — full path simulation (storing
  complete price history per path, not just terminal values), followed
  by a backward walk that regresses continuation value against current
  price for in-the-money paths at each step, comparing against
  immediate exercise value and updating each path's realized cash flow
  and exercise time accordingly

**Validated against:**
- Regression solver verified independently first: exact recovery of
  known quadratic coefficients from zero-noise data, before trusting it
  inside the full algorithm
- LSM price converges to the binomial tree's American price for a put
  with a nonzero dividend yield
- A second, deliberately different case (in-the-money call, different
  dividend yield/vol/expiry) also converges, confirming the result
  generalizes rather than being a coincidence of one parameter set

**Defend this:**
Can explain the core resolution to the apparent circularity of
American Monte Carlo: simulate all paths to expiry first, so the
"future" needed for the continuation-value estimate is already known
data, then walk backward using regression to estimate continuation
value from that data. Can explain why only in-the-money paths are
used in the regression, and why quadratic (not linear or higher-order)
basis functions were chosen. Can explain why this phase specifically
required full path-stepping, unlike European Monte Carlo, connecting
back to the explicit deferral noted in Phase 3. Can explain why the
regression solver was verified independently (on synthetic,
zero-noise data) before being trusted inside the much more complex
full algorithm — isolating linear-algebra correctness from Monte
Carlo/pricing correctness as two separate claims, not one conflated
test. Can explain why validating against a single parameter set isn't
sufficient evidence the algorithm is correct, and why a second,
meaningfully different case was added.

---

## Phase 7 — Finite difference PDE (Crank-Nicolson) (done)

**What was built:**
- `solve_tridiagonal()` — Thomas algorithm (forward elimination +
  back-substitution) for solving tridiagonal linear systems in O(n)
- `crank_nicolson_price()` — full PDE grid solver: known expiry
  payoff and boundary conditions seed the grid, walked backward in
  time via the Crank-Nicolson stencil, each time step solved as a
  tridiagonal system
- `explicit_fd_price()` — a deliberately simple fully-explicit scheme,
  built specifically to demonstrate conditional stability by contrast,
  not for production pricing use

**Validated against:**
- Tridiagonal solver verified independently first, against a
  hand-solvable known system, before being trusted inside the full
  PDE solver
- Crank-Nicolson price matches Black-Scholes for both European calls
  and puts
- Stability demonstrated directly: on an identical grid (fine price
  steps, coarse time steps — the specific combination that breaks
  explicit stability), the explicit scheme diverged to -4.06×10^22
  while Crank-Nicolson remained accurate at 10.4548 against a true
  price of 10.4506

**Defend this:**
Can explain the grid setup and why it walks backward from expiry,
structurally mirroring the tree and LSM despite being a fundamentally
different numerical method (PDE discretization vs. tree branching vs.
simulated regression). Can explain why Crank-Nicolson averages
explicit and implicit formulations, and precisely what that trades
off: unconditional stability plus second-order time accuracy, versus
explicit's simplicity-but-conditional-stability and implicit's
stability-but-lower-accuracy. Can explain why solving each time step
reduces to a tridiagonal system, and why the Thomas algorithm (O(n))
is used instead of general Gaussian elimination (O(n^3)) for it. Can
walk through the concrete stability demonstration and explain
precisely why that specific combination of fine price steps and
coarse time steps is what breaks the explicit scheme's stability
condition, rather than treating "unconditionally stable" as a fact to
recite without a mechanism.

---

## Phase 8 — Calibration & real market data (done)

**What was built:**
- `implied_volatility()` — Newton-Raphson solver using vega for the
  update step, with a near-zero-vega guard and sigma clamping for
  numerical safety
- Verified Deribit's public API is accessible (no auth, no US
  geoblocking on public endpoints) via direct curl test before
  committing to it as the data source
- `build_vol_surface_from_snapshot()` — parses a full real BTC options
  chain (JSON, via nlohmann/json), decodes Deribit's instrument naming
  convention (e.g. BTC-30OCT26-90000-P) into strike/expiry/type,
  converts BTC-denominated prices to USD, and runs the implied vol
  solver independently across every contract in the chain

**Validated against:**
- Round-trip: recovers a known volatility exactly from its own
  Black-Scholes price
- Edge case: deep OTM option with near-zero vega handled gracefully,
  no crash or nan
- Single real contract (BTC-30OCT26-90000-P): solver's independently
  computed implied vol (34.13%) landed within 0.23 percentage points
  of Deribit's own reported mark_iv (34.36%), converging in 4
  iterations
- Full live chain, not just one contract: 904 contracts parsed from a
  saved snapshot, 860 (95.1%) converged. Across converged contracts,
  average |solved_iv - market_iv| was 0.019 (1.9 percentage points)
  — confirming the single-contract result holds at market scale, not
  a coincidence of one hand-picked example
- Volatility smile visualized for the most liquid expiry in the
  snapshot (~99 days out, 118 contracts): a clean, textbook smile
  shape emerges directly from real market prices — implied vol
  bottoms out near-the-money (~36% around $84,000-86,000 strikes,
  close to spot) and rises on both wings (over 100% for deep OTM
  puts near $20,000, ~68% for deep OTM calls near $230,000)
- The 44 non-converged contracts (4.9%) were individually inspected,
  not just counted: confirmed to cluster exactly where the solver's
  known vega-collapse failure mode predicts (moneyness ratios from
  0.39 to 4.18, many under 8 days to expiry, some under 1 day).
  Their solved_iv values directly confirm the mechanism: many stuck
  at exactly 0.5 (the solver's starting guess — the vega guard
  triggered on the very first iteration) or exactly 0.01 (the clamp
  floor — one Newton step overshot before hitting the guard)

**Reproducibility:**
The chain snapshot is pulled once via curl and saved to
data/btc_chain_snapshot.json rather than tests hitting Deribit's live
API — a live dependency would make tests flaky and non-deterministic
as the real market moves. today_year/month/day are passed explicitly
into the surface builder rather than read from the system clock, so
time-to-expiry is pinned to when the snapshot was taken.

**Defend this:**
Can explain why implied vol requires a numerical solver rather than a
closed-form formula, and specifically why vega is the right quantity
to drive Newton-Raphson's update step. Can explain the near-zero-vega
failure mode concretely, connecting it back to the exact edge cases
(deep ITM/OTM, near-expiry) already identified and tested in Phase 1 —
this wasn't a new discovery, it was an already-known risk being
handled defensively in a new context. Can explain, with specifics, why
a small gap against Deribit's own number is expected rather than a
correctness failure: coin-margined quanto effects, timing mismatch
between the snapshot and Deribit's internal calculation, and
Deribit's own smoothing across their internal vol curve. Can explain
why Deribit's public data was verified directly (via curl) before
being adopted, rather than assumed accessible. Can explain why
validating across the full 904-contract chain is meaningfully
stronger evidence than validating a single hand-picked contract. Can
explain the volatility smile's shape in terms of market-implied tail
risk, and connect it directly back to Black-Scholes' constant-
volatility assumption flagged as a known limitation in Phase 1. Can
walk through the non-convergence investigation end to end — not just
"44 failed," but why, demonstrated using the actual failure values
(0.5 and 0.01) as direct evidence of which specific guard triggered
for each case. Can explain why a snapshot-based approach was used
instead of live API calls in tests.

---

## Phase 9 — Heston stochastic volatility (done)

**What was built:**
- `HestonParams` (v0, kappa, theta, xi, rho) and `heston_price()` —
  European option pricing under Heston (1993) via the COS method
  (Fang & Oosterlee, 2008): a truncated cosine-series expansion of the
  terminal log-return density, built from the model's characteristic
  function, evaluated in the numerically stable "little trap"
  formulation (Albrecher, Mayer, Schoutens, Tistaert, 2007) rather
  than the original 1993 formula, which suffers branch-cut
  discontinuities in the principal complex log/sqrt as tau grows
- Prices European puts directly (COS payoff coefficients integrated
  over `[a, 0]`) and derives calls via put-call parity, rather than
  pricing calls directly — a deliberate direction choice, not an
  arbitrary one (see bug note below)
- 256-term COS summation, truncation range `[a, b]` sized generously
  from the model's own variance estimate rather than exact Heston
  cumulants, so the series converges without needing error-prone
  closed-form cumulant formulas transcribed from memory

**Validated against:**
- Deterministic-variance limit: with vol-of-vol xi forced to ~0, the
  variance process has (almost) no randomness and Heston must
  degenerate exactly to Black-Scholes with sigma = sqrt(v0) — checked
  across 9 scenarios spanning ATM/OTM, short/long-dated, and BTC-scale
  parameters, catching any sign error in the characteristic function
  that a numerical-stability artifact alone wouldn't reveal
- Independent Monte Carlo cross-check: full-truncation Euler
  simulation of the actual Heston SDEs (`dS`, `dv` stepped directly,
  not the closed-form CF) agrees with the COS-method price within
  statistical tolerance — confirms the closed-form isn't just
  internally consistent but actually solves the model it claims to
  price
- Heston's own put-call parity relation, checked independently of
  Black-Scholes, across parameter sets that actually exercise
  stochastic vol (nonzero xi and rho)
- Monotonicity: call price increases with initial variance v0, holding
  everything else fixed
- Edge cases: reduces to intrinsic value as expiry approaches; stays
  finite when the Feller condition (`2*kappa*theta > xi^2`) is
  violated, which real BTC-calibrated parameters routinely do

**Bugs found and fixed:**

This model went through two real numerical bugs before the pricer was
trustworthy across the full range of moneyness a real BTC options
chain actually has. Both were found the same way: a discrepancy
between Heston and Black-Scholes in the deterministic-variance limit,
where the two *must* agree exactly, that was invisible near the money
and only showed up far from it.

1. **Catastrophic cancellation in the characteristic function's log
   term (the real bug).** `heston_log_return_cf` computes a term
   `log((1 - g*E) / (1 - g))` where `g` is O(xi^2) — tiny for any
   small-to-moderate vol-of-vol, not just as xi→0. That makes the
   ratio inside the log extremely close to 1.0, and `std::log(1.0 +
   tiny)` loses most or all of its meaningful digits: `"1.0 - g*E"`
   and `"1.0 - g"` both round to (or very near) 1.0 in double
   precision *before* `log()` ever runs. This error is normally too
   small to notice — a tiny relative error in an intermediate
   quantity, lost in the noise for near-the-money options. It stopped
   being invisible for far-from-the-money options priced through
   put-call parity: parity derives a *small* target price as the
   difference of two *large* numbers (the put and the discounted
   forward), so the target's accuracy is bounded by the put's
   *absolute* accuracy, not its relative accuracy. Concretely: at
   spot=$100,000, strike=$300,000, 90 days, 60% vol (strike = 3x spot,
   a realistic scenario for a crypto options chain, found via manual
   testing of the live demo) this cancellation alone produced an 82%
   pricing error in the deterministic-variance limit, which must equal
   Black-Scholes exactly. **Fix:** rewrote the term as `1 +
   g*(1-E)/(1-g)` and evaluated its log via a hand-built complex
   `log1p` (`clog1p()` in `heston.cpp`) — the same technique
   `std::log1p` uses for reals, splitting the real part into
   `0.5*log1p(2*Re(z) + |z|^2)` (no subtraction of near-equal values)
   and the imaginary part into `atan2(Im(z), 1+Re(z))`. Reduced the
   error from 82% to ~0.001%, confirmed via the extreme-moneyness
   regression case added to the deterministic-limit test.

2. **A `long double` false start (not a fix — a mistake worth
   documenting).** The first attempt at fixing bug #1 was to switch
   `heston.cpp`'s complex type from `std::complex<double>` to
   `std::complex<long double>`, reasoning that more bits of precision
   would simply absorb the cancellation. This *appeared* to work: it
   fixed the extreme-moneyness test in this project's Linux/x86
   development sandbox, and even held up when cross-compiled to
   WebAssembly via Emscripten (which gives `long double` full IEEE
   quad precision). It was shipped on that evidence. It did not
   actually fix anything: on Apple Silicon (ARM64) Macs, `long double`
   is bit-for-bit identical to `double` — there's no extended-precision
   hardware to use — so the "fix" was a complete no-op there. This was
   only caught because the change was verified on the actual target
   hardware afterward, which reproduced the *exact* same wrong value
   (0.2907348679 instead of 1.6454796818) as before the "fix." The
   real fix (`clog1p`, above) uses only ordinary `double` throughout
   and is genuinely portable, because it fixes the *cancellation*
   rather than trying to out-precision it. **Lesson documented,
   deliberately, not smoothed over:** a numerical fix verified on one
   platform (even two, if both happen to have the same underlying
   float behavior) is not verified — `long double`'s width is
   implementation-defined per the C++ standard, and x86-64 and ARM64
   disagree on it in exactly the way that matters here. The fix that
   shipped doesn't rely on any platform-specific float behavior.

**Defend this:**
Can explain the COS method mechanically: it expands the terminal
log-return density as a truncated Fourier-cosine series on `[a, b]`,
with coefficients recovered directly from the characteristic function
via a cosine transform, avoiding the numerical Fourier inversion
integral that naive characteristic-function pricing would otherwise
need. Can explain why puts are priced directly and calls via parity,
not the reverse: the call payoff's COS coefficients integrate over
`[0, b]`, and `b` scales with total variance — for BTC-realistic
vol/tenor combinations `b` routinely lands in the 8-20 range, making
`exp(b)` astronomically large inside each coefficient, which an
N-term sum then has to cancel back down to an O(spot) answer;
observed failure mode was silently wrong prices (from -36% to +∞
depending on truncation width), not a crash, confirmed independent of
summation length (ruling out under-convergence as the cause). Puts
integrate over `[a, 0]` instead, where `exp(0)=1` and `exp(a)` is
tiny, so this instability doesn't arise regardless of truncation
width. Can explain both cancellation bugs at the mechanism level, not
just "floating point is hard": the `(A-d)` cancellation (an earlier,
separately-fixed issue in the same characteristic function, resolved
via the algebraic identity `A²-d² = -xi²(iu+u²)` so `(A-d)` is
computed as a sum in the denominator rather than a literal
difference) and the `log(1+tiny)` cancellation documented above are
two *different* hazards in the same formula, both invisible near the
money and both bounded by absolute, not relative, precision once
routed through put-call parity. Can explain precisely why `long
double` seemed to work and precisely why it didn't really: x86-64's
80-bit extended format and Emscripten's IEEE-quad `long double` both
happen to add real precision over `double`, while ARM64's ABI defines
`long double` as identical to `double` — so a fix that works by
adding unspecified extra bits is only ever accidentally portable. Can
state, without hedging, that this was caught by testing on real
target hardware, not by reasoning about the standard in advance — and
that the actual fix doesn't depend on `long double`'s width being
anything in particular, on any platform.

---

## Phase 10 — Trinomial tree (done)

**What was built:**
- `trinomial_tree_price()`, overloaded for `EuropeanOption` and
  `AmericanOption`, alongside the existing `binomial_tree_price()` in
  `tree_pricer.cpp` — the Boyle (1986) trinomial lattice: each node
  branches to three nodes (up, middle/unchanged, down) at the next
  time step instead of the binomial tree's two
- Up/down step size `dx = sigma*sqrt(3*dt)` and branch probabilities
  `pu`, `pm`, `pd` derived by moment-matching (choosing the three
  probabilities so the discretized process's mean and variance over
  one step match the true lognormal process's) — the same idea behind
  the binomial tree's `u`/`d`/`p`, just with a third free parameter
  (the middle branch) to match with
- Wired into the WASM demo (`wasm_bridge.cpp` gained
  `bridge_trinomial_price_european`/`_american`, `docs/index.html`
  gained a "Trinomial tree" result card next to "Binomial tree",
  respecting the same exercise-style toggle) — verified against the
  actual compiled WASM module in Node before shipping, not just the
  native build, given this codebase's own recent history with
  platform-dependent floating-point behavior (see Phase 9)

**Validated against:**
- Converges to Black-Scholes for both European calls and puts: ~2%
  gap at 10 steps (same slack the binomial tree test allows), ~0.1%
  at 500 steps
- American put priced above the corresponding European put under a
  nonzero dividend yield, same property the binomial tree is checked
  against
- **Agrees with the binomial tree directly**, not just via a shared
  Black-Scholes reference — checked across three scenarios (ATM call,
  OTM put with dividends, BTC-scale high-vol ITM call) at 800 steps
  each, within 0.2%. This is meaningfully stronger evidence than each
  tree separately matching Black-Scholes: `run_tree()` and
  `run_trinomial_tree()` are two completely independent
  implementations (different state-space indexing, different
  moment-matching derivation for their probabilities), so an error in
  one wouldn't be expected to also appear in the other. Two
  independently-derived numerical methods landing on the same answer
  rules out a shared mistake in a way that one method matching a
  third, unrelated closed-form reference doesn't fully rule out.

**Measured, not assumed:** does the trinomial tree actually converge
faster per step than the binomial tree, the way it's often described?
Checked directly (`[.manual]` test, same "measure it, don't just
claim it" approach as Phase 4's AAD benchmark) at steps = 10, 25, 50,
100, 200, 400 on the standard reference case (S=K=100, T=1, r=5%,
vol=20%). Result: trinomial error was lower at 5 of 6 step counts,
by roughly 5-10% — a modest, real improvement, not a different
convergence order. Both trees still shrink error at essentially the
same O(1/N) rate; the trinomial tree's extra branch buys a smaller
constant, not an asymptotically better exponent. Presented honestly
as that: a modest, measured improvement, not "3 branches is 50%
better" or similar oversold claim.

**Defend this:**
Can explain why a third (middle) branch needs a third condition to
pin down: two branches (binomial) need only match the process's mean
and variance, which is exactly two equations for two unknowns (`u`
and `p`, with `d=1/u` fixed by the recombination requirement); three
branches (trinomial) have three free probabilities, so a third
condition is needed — here, that's `dx` itself being chosen up front
(`sigma*sqrt(3*dt)`, not derived from moment-matching), leaving
`pu`/`pm`/`pd` to solve the two moment equations plus the
probabilities-sum-to-one constraint. Can explain why this tree was
verified against the actual compiled WASM artifact via Node, not just
the native Linux build, connecting directly back to the Phase 9
lesson: a numerical method "working" on one build target is not
evidence it works everywhere, and this project no longer treats a
single platform's test run as sufficient for anything touching
floating-point math. Can explain, honestly and without inflating the
result, that the trinomial tree here is not a strictly-better
replacement for the binomial tree — it's roughly 5-10% more accurate
per step at real, measured cost of running 2-3x more nodes per step
(2n+1 vs n+1 at depth n) — and can state plainly why both were kept
side by side in the live demo rather than the trinomial tree
replacing the binomial one: showing two independently-derived methods
agree is itself part of this project's correctness evidence (see
Validated against, above), and that evidence is weaker with only one
tree in the picture.

---

## Phase 11 — Second-order AAD (gamma) (done)

**What was built:**
- `Dual` (`dual.hpp`): a forward-mode dual number (`val + dot*eps`,
  `eps^2=0`), with arithmetic (`+ - * /`) and the two special
  functions (`dual_exp`, `dual_sqrt`) the Monte Carlo path/Black-
  Scholes formula below need
- `Tape2`/`ADouble2` (`tape2.hpp`/`adouble2.hpp`): a second, fully
  isolated sibling of the original `Tape`/`ADouble` (Phase 4) — same
  structure, same backward-walk logic, but every value and every local
  partial derivative is a `Dual` instead of a plain `double`. This is
  "forward-over-reverse" AD: seed ONE input's `Dual.dot = 1.0` (the
  direction to differentiate twice in), run the ordinary reverse-mode
  backward pass, and the resulting adjoint's own `.dot` component is
  the second derivative — computed exactly, in the same single pass,
  no finite differencing anywhere. Also gained `ad2_log` and
  `ad2_norm_cdf` (the standard normal CDF as a special function with a
  hand-supplied derivative, `N'(x) = phi(x)`, mirroring how `ad_exp`/
  `ad_sqrt` already worked in the original tape) — needed to run this
  machinery through the Black-Scholes formula itself, not just simple
  arithmetic
- Kept entirely separate from the original `Tape`/`ADouble`
  (duplicated rather than templated) so this new, less-battle-tested
  machinery could not silently change any already-shipped, already-
  tested first-order Greek if it had a bug
- Two real uses, wired up deliberately to show contrasting outcomes:
  `black_scholes_greeks_via_ad2()` (gamma on the smooth closed-form
  formula) and `monte_carlo_gamma()` (gamma on a Monte Carlo path's
  payoff) — see Validated against, below, for why these two land in
  very different places

**Validated against:**
- Hand-computed toy examples FIRST, before any real pricing code
  touched this machinery — same discipline Phase 4 used for the
  original tape (`y = a*b + c`). Checked `f(x)=x^3`, `f(x)=1/x`,
  `f(x)=exp(x)`, `f(x)=sqrt(x)`, `f(x)=log(x)`, `f(x)=N(x)` (the
  normal CDF), and a two-variable case `f(x,y)=x^2*y+e^x` — all first
  AND second derivatives matched their closed-form analytical values
  to machine precision (~1e-16), not just a statistical tolerance
- `black_scholes_greeks_via_ad2()`'s gamma matches
  `black_scholes_greeks()`'s independently-computed analytical gamma
  to ~1e-9 across 6 scenarios (ATM, OTM with a dividend, deep ITM,
  BTC-scale short-dated, and right at the edge of expiry) — the
  formula is smooth in spot, so this is exact agreement, not a
  statistical check
- `monte_carlo_gamma()`'s delta (a first-order quantity, computed on
  the same Dual-valued tape as an incidental byproduct) cross-checked
  against `black_scholes_greeks()`'s analytical delta, confirming the
  new machinery's ordinary first-order output is correct even where
  its second-order output (gamma) isn't usable — see the finding below

**The real finding — pathwise AAD gamma is exactly zero for Monte
Carlo, and that's not a bug:**
Running the exact same forward-over-reverse technique through a Monte
Carlo path's payoff instead of the smooth Black-Scholes formula
returns gamma = 0.0 — not approximately zero, not noisy-near-zero,
exactly `0.0`, on every run, at every path count tested. This is the
real, well-documented reason pathwise differentiation (the technique
Phase 4's `monte_carlo_greeks` already uses successfully for delta,
vega, theta, rho) does not extend to gamma. A single simulated path's
payoff, `max(S_T - K, 0)`, is a PIECEWISE LINEAR function of spot:
along any one path it's either exactly `S_T - K` (slope 1, zero
curvature) or exactly `0` (zero everything) — its second derivative is
identically zero everywhere except exactly at the kink point
(`S_T = K`), which has probability zero of ever being sampled under
continuous GBM. First-order pathwise differentiation works because the
payoff function itself is continuous even though it has a kink (Phase
4 already established this: "branching on the plain value, valid since
the kink has zero probability of being hit exactly"); second-order
differentiation of that same kinked function is fundamentally
different, because differentiating twice asks a question about local
curvature that a piecewise-LINEAR function structurally cannot answer
away from its kink. A working Monte Carlo gamma estimator exists in
the literature (the likelihood-ratio/score-function method, which
differentiates the probability density instead of the payoff, sidestepping
this problem entirely) but is a different technique, not an extension
of this one — deliberately out of scope here, and documented as a real
finding rather than silently shipping a `monte_carlo_gamma()` that
always returns zero without explaining why.

**Defend this:**
Can explain forward-over-reverse AD mechanically: layering a forward-
mode dual number underneath a reverse-mode tape computes, in the SAME
backward pass that produces first derivatives, a Hessian-vector
product — the second derivative in whatever single direction was
seeded — without needing to literally differentiate the backward pass
by hand (a real, separate undertaking Phase 4 correctly identified as
"genuinely larger"). Can explain precisely why this was built as a
fully separate `Tape2`/`ADouble2` rather than templating the original
`Tape`/`ADouble`: isolating blast radius, so a bug in newer, less-
proven machinery can't silently corrupt an already-shipped, already-
tested Greek. Can explain, with the actual numbers, why gamma works to
near machine precision on the Black-Scholes formula (smooth function,
no kink) and returns exactly zero on a Monte Carlo path (piecewise-
linear payoff, kink) — and can state clearly that this is a structural
fact about pathwise differentiation, not a limitation of this specific
implementation, connecting it to the general principle: an estimator
valid for one order of differentiation isn't automatically valid for
the next order up, and checking that (rather than assuming it) is
what this phase actually did. Can name the real fix (likelihood-ratio
method) and explain in one sentence why it works where pathwise
doesn't: it differentiates the path's PROBABILITY of occurring, not
the payoff's VALUE, so it never needs the payoff itself to be twice
differentiable.

## Phase 12 — Multi-expiry Heston calibration (done)

**What was built:**
- Discovered, rather than wrote: `calibrate_heston()` (Phase 8/9) was
  already generic over expiries — it prices every point at that point's
  own `time_to_expiry`, with no single-expiry assumption anywhere in its
  implementation. The single-expiry framing in its original header
  comment was a usage decision, not a technical limitation. Updated that
  header comment (`heston_calibration.hpp`) to say so precisely, instead
  of writing a second, parallel calibration function that would have
  duplicated working code
- A synthetic identifiability test
  (`"Multi-expiry Heston calibration resolves the single-expiry
  kappa/theta non-identifiability"`, `tests/test_main.cpp`) that makes
  the Phase 8 non-identifiability finding concrete on both sides: the
  same ground truth, the same deliberately-bad initial guess, calibrated
  once against one synthetic expiry and once against three pooled
  synthetic expiries
- `tools/calibrate_heston_multi.cpp`: a real-data counterpart to
  `tools/calibrate_heston.cpp`, pooling the ~43-day and ~99-day expiries
  from the real Deribit snapshot (same moneyness/liquidity filter) into
  one joint fit, instead of one expiry alone

**Validated against:**
- Synthetic ground truth (`v0=0.16, kappa=2.0, theta=0.09, xi=0.5,
  rho=-0.6` — `v0 != theta` deliberately, to create a genuine variance
  term structure for multiple expiries to carry different information
  about; if `v0 == theta`, variance never mean-reverts and kappa stays
  unidentifiable regardless of how many expiries are pooled). From the
  same bad guess (`kappa=6.0` against a true `2.0`): single-expiry
  calibration reprices its one synthetic expiry to <0.01% RMSE while
  recovering a kappa over 200% off the truth (the identifiability
  failure, isolated); multi-expiry calibration, same guess, three pooled
  expiries, recovers all 5 parameters to within floating-point noise of
  the truth
- The real Deribit chain, pooling the 43-day (94 contracts) and 99-day
  (66 contracts) slices under the same moneyness filter
  (`tools/calibrate_heston.cpp` uses): converges to a fit at 1.35
  percentage points RMSE in vol-space — worse than the single-expiry
  tool's 0.43pp on its own slice, but still better than the flat-vol
  headline result's ~1.9pp mean gap across the whole chain

**The honest finding — multi-expiry calibration resolves identifiability
on synthetic data, and exposes a real limitation on real data:**
Pooling expiries is not a strictly-better replacement for the
single-expiry + regularization approach Phase 8 used; it trades one
problem for a different, more informative one. On synthetic data
generated FROM Heston with fixed parameters, pooling expiries recovers
those exact parameters, because the data-generating process really is
constant-parameter Heston — there's nothing for the model to fail to
capture. On the real chain, the market's actual volatility dynamics are
not exactly constant-parameter Heston, and forcing a joint fit across
multiple real maturities makes that mismatch visible: fit quality
degrades as more expiries are pooled — 2 expiries (43d/99d): 1.35pp;
3 expiries (43d/99d/190d): ~2.2pp; 4 expiries (43d/99d/190d/281d):
~2.9-3.1pp, measured during development, not included in the shipped
tool. This is the expected, well-known behavior of a single-factor
stochastic-vol model asked to fit
a full term structure with one fixed parameter set — not a bug, and not
something a better optimizer or more iterations would fix. Properly
resolving it needs a term-structure extension (piecewise-constant
parameters per maturity bucket, or a multi-factor model like Bergomi),
deliberately out of scope here, consistent with every other scope
decision in this project.

**Why the real-data multi-expiry result didn't replace the demo's Heston
card:** the single-expiry number (0.43pp) is a clean, strong, honest
answer to "does Heston fit today's smile well?" The multi-expiry number
(1.35pp) is an honest answer to a different, more nuanced question ("can
one Heston parameter set explain several maturities at once?") that
doesn't compress into a single headline stat without the explanation
above. Rather than either hiding the weaker number or overselling it,
it's documented here and kept as a standalone tool + test, with the
finding written down plainly — see `docs/numerics.md`.

**Defend this:**
Can explain precisely why `calibrate_heston()` needed no code changes to
support multi-expiry calibration: it already prices each point against
its own `time_to_expiry`, so "single-expiry" was always just a property
of what the caller passed in, not of the function. Can explain the
identifiability mechanism itself: a single expiry's option prices
constrain kappa and theta only through an integrated combination of the
two (roughly, the average variance over that one time horizon), so many
(kappa, theta) pairs with the same integrated effect price that expiry
identically; multiple expiries at different horizons sample that average
at different stages of mean-reversion, which is what breaks the tie. Can
explain, with the actual measured numbers, why the real-data result gets
WORSE as more expiries are added rather than converging to a clean
answer — and can state directly that this is the correct, expected
behavior of a single-factor model meeting real term-structure dynamics
it cannot represent, not a defect in the fitting code. Can name the real
fix (piecewise-constant parameters or a multi-factor model) and explain
why it's out of scope: it's a genuinely different model, not a tuning
change to this one.