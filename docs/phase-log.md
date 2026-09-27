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
## Phase 13 — Exotic options: Asian and Barrier (done)

**What was built:**
- `AsianOption` and `BarrierOption` instruments (`instrument.hpp`), plus
  a `BarrierDirection` enum covering all four combinations (down/up,
  in/out)
- `exotic_options.hpp`/`exotic_options.cpp`: a full-path-simulation Monte
  Carlo pricer for arithmetic-average Asian options
  (`asian_option_price_mc`), a geometric-average sibling using the same
  simulator (`geometric_asian_price_mc`) built purely to validate against
  a closed form, a Kemna-Vorst closed-form geometric Asian pricer
  (`geometric_asian_price_closed_form`), a Monte Carlo barrier pricer
  covering all four directions (`barrier_option_price_mc`), and a
  method-of-images closed-form down-and-out call pricer
  (`down_and_out_call_price_closed_form`)
- Both Monte Carlo pricers reuse the Euler-exact GBM path-stepping
  pattern already established in `longstaff_schwartz.cpp`, rather than
  reinventing it — genuinely needed here (unlike vanilla European Monte
  Carlo, whose terminal-jump sampling is the deliberate, documented
  choice `docs/numerics.md` explains) because both payoffs depend on the
  whole path, not just `S_T`
- Seven new tests (`tests/test_main.cpp`, `[exotic]` tag)

**Validated against:**
- Geometric Asian: Monte Carlo (300,000 paths) vs. the Kemna-Vorst closed
  form — 5.6680 ± 0.0175 vs. 5.6411, within 1.5 standard errors
- AM-GM inequality: arithmetic Asian call price (5.8843) measured
  strictly above the geometric Asian call price on the same simulated
  paths — provable pathwise, not just observed
- Asian call price (5.8843) measured strictly below the vanilla European
  call price (10.4506) on the same market, consistent with averaging
  reducing effective volatility
- `num_steps=1` Asian option vs. vanilla European Monte Carlo, same seed,
  same dynamics — collapse to the same price within combined standard
  error, the provable single-monitoring-date edge case
- Down-and-out barrier call: Monte Carlo (300,000 paths) vs. the
  method-of-images closed form — 10.3526 ± 0.0329 vs. 10.3513
- In/out parity, both directions: down-and-out + down-and-in = 10.4265
  vs. vanilla = 10.4506; up-and-out + up-and-in = 10.4265 vs. vanilla =
  10.4506 — both within combined standard error, validating three of the
  four barrier directions without needing closed forms for any of them

**Why the general Reiner-Rubinstein barrier formula wasn't implemented:**
the full closed-form barrier pricing formula (Reiner & Rubinstein, 1991)
covers all 16 direction/strike-vs-barrier configurations, but the
case-selection logic across those 16 cases is exactly the kind of thing
that's easy to get subtly wrong and hard to catch by testing — and there
was no independent reference implementation on hand here to check it
against. Scoped down instead to the one case implementable with real
confidence (down-and-out, `H <= min(S,K)`, where a single reflection is
exact), and used the model-independent in/out parity identity — which
needs no formula derivation at all, just the fact that a knock-out and
its knock-in partition every path — to validate the other three
directions instead. A smaller, fully-trustworthy validation surface over
a larger, riskier one.

**Follow-up: wired into the live WASM demo in its own card.** Shipped
initially without touching the demo page (a deliberate first-cut scope
decision, keeping the live UI untouched and low-risk while the pricing
code itself was new). Once the engine side was verified working, a
second "Exotic options" card was added (`docs/index.html` section 07,
`wasm_bridge.cpp`'s `bridge_asian_price`/`bridge_barrier_price`/
`bridge_geometric_asian_closed_form`/`bridge_down_and_out_call_closed_form`),
with its own inputs (monitoring steps; barrier level and direction for
Barrier) rather than forcing them into the existing calculator's layout.
The card runs the same validation checks documented above live in the
browser — geometric Asian MC against its closed form, and for Barrier,
the selected direction against its same-barrier opposite side with the
sum checked against vanilla, plus the down-and-out closed form when
valid — so a visitor sees the correctness argument, not just a price.

**Defend this:**
Can explain why vanilla European Monte Carlo samples only the terminal
price while Asian and Barrier pricing need the full path: the payoff
itself is path-dependent (an average, or a touch condition), not just a
function of `S_T`, so there's no shortcut past simulating the
intermediate steps. Can explain why the geometric Asian option has a
closed form and the arithmetic one doesn't: a sum of lognormal random
variables isn't itself lognormal (arithmetic average), but a product of
lognormals is (geometric average, via the log-sum), which is the same
algebraic fact that makes Black-Scholes solvable in the first place. Can
explain the method-of-images formula's `H <= min(S,K)` constraint
precisely: the reflection only cancels the payoff exactly at the barrier
when the payoff is already zero there, which is guaranteed for a call
struck at or above a barrier below spot, and breaks down otherwise. Can
explain why in/out parity works without any model assumptions: it's true
by simple case exhaustion (every path either touches the barrier or
doesn't), not a Black-Scholes-specific or even GBM-specific result — it
would hold under any stochastic process. Can name the honest limitation:
no general (all-16-case) barrier formula, and no continuous-monitoring
correction (this prices discretely-monitored barriers at `num_steps`
dates, which is what's actually simulated — a continuously-monitored
barrier is a different, slightly more barrier-sensitive instrument, and
pricing it exactly needs a discrete-to-continuous correction this
implementation doesn't apply).

## Phase 14 — Bates model: Heston + Merton jumps (done)

**What was built:**
- `bates.hpp`/`bates.cpp`: `BatesParams` (a `HestonParams` plus
  `jump_intensity`, `jump_mean`, `jump_vol`) and `bates_price()`, pricing
  European options under Bates (1996) — Heston stochastic volatility with
  a compound Poisson jump process (Merton 1976's lognormal jump-size
  assumption) layered on the log-price. Prices via the same COS method
  (Fang & Oosterlee, 2008) `heston_price()` uses, not a different
  numerical technique
- `cos_pricing_internal.hpp` (new): `heston_log_return_cf`, `chi`, `psi`,
  and the `clog1p` numerical-stability helper, pulled out of `heston.cpp`
  verbatim into a shared `deriv::detail` header so `bates.cpp` could
  reuse Heston's numerically-hardened characteristic function directly —
  called with an adjusted drift (`r - lambda*k`) rather than re-derived —
  instead of copy-pasting that fragile complex-analysis code a second
  time. Confirmed behavior-preserving: `heston.cpp`'s own tests
  (50 assertions) pass unchanged after the refactor
- The jump component is a separate multiplicative term in the
  characteristic function (`jump_log_return_cf`, the standard
  compound-Poisson CF identity `exp(lambda*tau*(phi_X(u)-1))`), so
  `bates_log_return_cf(u) = heston_log_return_cf(u, r-lambda*k) *
  jump_log_return_cf(u)` — no re-derivation of Heston's own CF math, just
  composition
- `bates_mc.hpp`/`bates_mc.cpp`: an independent Monte Carlo cross-check,
  mirroring `heston_mc.cpp` exactly for the variance-process stepping
  (full-truncation Euler) and adding a compound Poisson jump
  (`std::poisson_distribution(lambda*dt)` jump count,
  `std::normal_distribution(jump_mean, jump_vol)` jump sizes) to the
  log-price increment each step
- Nine new tests (`tests/test_main.cpp`, `[bates]` tag): jump-intensity-
  zero collapses to Heston exactly, the deterministic-variance limit
  collapses to closed-form Merton jump-diffusion, and the independent SDE
  Monte Carlo agrees with the COS closed form

**Validated against:**
- Heston collapse (`jump_intensity = 0`): exact match to
  `heston_price()`, margin `1e-9`, across ATM/OTM/call/put/BTC-scale
  scenarios — with nonzero `jump_mean`/`jump_vol` left in the params
  deliberately, to confirm they're inert (not just that "all zeros"
  happens to work)
- Deterministic-variance limit (`xi = 1e-6`, `v0 = theta`) vs. an
  independent closed-form Merton (1976) series
  (`merton_jump_diffusion_price()` in the test file — a Poisson-weighted
  sum of Black-Scholes prices, sharing no code with `bates.cpp` or
  `heston.cpp`): agreement within 0.1% across five scenarios spanning
  ATM/OTM, positive/negative mean jump, and BTC-scale frequent jumps
- Independent SDE + jump Monte Carlo (60,000 paths, antithetic, fixed
  seed) vs. the COS closed form: agreement within 5 standard errors
  across three scenarios, including a BTC-like high-vol-of-vol,
  frequent-jump case (COS 9926.93 vs. MC 9840.09 ± 70.23 — comfortably
  inside the bound)

**A bug found in the test oracle, not in the model — worth documenting
honestly:** the first version of the closed-form Merton series
(`merton_jump_diffusion_price()`) disagreed with `bates_price()`'s
deterministic-variance-limit output by about 0.55% (13.357939 vs.
13.431517) — small enough to look like truncation error, too large to
actually be truncation error. Widening the COS truncation window (`L`
14→30) and term count (`N` 256→2048) left `bates_price()`'s output
unchanged to the last printed digit, ruling that out immediately. What
found the real bug was building a from-scratch, fully independent
verification chain rather than trusting either side: (1) both
`heston_log_return_cf` and the full `bates_log_return_cf` were checked
directly against their hand-derived theoretical characteristic functions
— diffs of `1e-7`-`1e-9`, i.e. correct; (2) a standalone COS pricer using
the exact theoretical CF (zero dependency on `bates.cpp`/`heston.cpp`)
reproduced `bates_price()`'s number exactly, proving the COS/quadrature
logic itself wasn't the problem; (3) the decisive check — a direct
20-million-path Monte Carlo simulation of the actual Merton jump-diffusion
SDE (`std::poisson_distribution` + `std::normal_distribution`, no CF, no
COS series anywhere) gave 13.360144, matching `bates_price()`'s
13.357939 within Monte Carlo noise and matching the "Merton series"
formula nowhere close. That pinned the bug to the series formula, not the
model. Root cause: `black_scholes_price()` called with `rate=r_n` (the
jump-adjusted rate each series term needs for its forward/terminal
distribution) also silently uses `r_n` as the discount rate — but the
true discount rate for every term is always `r`, not `r_n`. The fix
rescales each term's raw Black-Scholes output by `exp((r_n - r) * T)`
before summing, which swaps the wrong `e^{-r_n*T}` discount factor for
the correct `e^{-r*T}` one while leaving the `r_n`-driven forward alone —
after which the series matches `bates_price()` exactly (13.357939 both
ways). The lesson generalizes past this one function: a closed-form
"reference" formula is not automatically more trustworthy than the model
under test just because it's independent code, and the way to actually
resolve a small, real-looking discrepancy is to add a third, structurally
different check (here, the direct SDE Monte Carlo) rather than assume
which of the two disagreeing numbers is correct.

**Why jumps get a separate CF term instead of folding into Heston's own
variance process:** a jump is a discontinuous, instantaneous move in the
log-price — nothing in Heston's SDEs (`dS`, `dv`) can produce that no
matter how `xi` or `kappa` are tuned, since both are continuous
diffusions. The compound Poisson process is the standard, minimal way to
add genuine discontinuities on top of a diffusion, and because Poisson
jump arrivals and Brownian diffusion are independent processes by
construction, their joint characteristic function is exactly the product
of the two marginal characteristic functions — which is what makes
`bates_log_return_cf` a one-line multiplication rather than a new
derivation.

**Defend this:**
Can explain why `bates_price()` calls `heston_log_return_cf` with
`r - lambda*k` instead of the raw rate `r`: the jump compensator `k =
E[Y] = exp(jump_mean + 0.5*jump_vol^2) - 1` is subtracted from the drift
so that `E[S_T] = S0*exp((r-q)*T)` still holds exactly — the same
no-arbitrage martingale condition Heston/Black-Scholes satisfy without
jumps — even though the jump process on its own has a nonzero expected
log-price contribution. Can explain why puts are priced directly and
calls via put-call parity, same as Heston: the put's COS payoff
coefficients integrate over `[a, 0]` and stay `O(strike)` regardless of
truncation width, while a far-OTM call priced directly (or derived from a
far-ITM put) can lose precision to catastrophic cancellation — this
project already found and fixed that exact bug once, in Heston, and Bates
reuses the fix rather than reintroducing the risk. Can name the honest
limitation: jump sizes are assumed i.i.d. lognormal (Merton's own
assumption) with a single intensity/mean/vol triple — no
time-varying jump intensity, no separate up/down jump distributions, and
no jumps in the variance process itself (some later extensions, e.g.
Duffie-Pan-Singleton, add those; this implementation deliberately doesn't
chase that generality without a concrete need for it).

## Phase 15 — Historical delta-hedging backtest (done)

**What was built:**
- `price_history.hpp`/`price_history.cpp`: loads a plain "date,price"
  daily CSV into an ordered vector -- no JSON dependency needed for this,
  unlike the Deribit snapshot data
- `data/btc_price_history.csv`: 5,904 daily BTC/USD closes,
  2010-07-18 through 2026-09-15, sourced from a public GitHub dataset
  (`Habrador/Bitcoin-price-visualization`) since this project's own
  Deribit pull only ever captured a single live snapshot, not years of
  history, and no free historical-options API was reachable from this
  environment (checked: CoinGecko, Yahoo Finance, and Stooq were all
  blocked by either robots.txt or this environment's egress policy; the
  GitHub-hosted CSV was the one path that actually worked). Spot-checked
  against known price history before trusting it: Dec 2017 ATH ~$19,343,
  Nov 2021 ATH ~$67,292, the Mar 2020 COVID crash to ~$5,800, and the Nov
  2022 FTX-collapse low ~$15,787 all land in the right place; zero date
  gaps or duplicates across the full series
- `backtest.hpp`/`backtest.cpp`: a model-agnostic rolling-window
  delta-hedging backtest (`run_backtest()`) -- rolls a fresh
  at-the-money European call every `window_days` calendar days, hedges
  it daily using the REAL historical price path (not a simulated one,
  unlike `hedging.cpp`'s Monte Carlo hedge-error test), and records
  realized hedge P&L at expiry, normalized to basis points of that
  window's starting spot (raw dollar P&L is meaningless to average
  across a series that spans $0.08 to $78,000)
- `tools/run_backtest.cpp`: loads the price history and runs all three
  models (Black-Scholes, Heston, Bates) through the same backtest,
  emitting JSON -- same "never hand-type real numbers" convention as
  `gen_stats.cpp`/`calibrate_heston.cpp`
- 7 new tests (`tests/test_main.cpp`, `[backtest]` tag): `trailing_realized_vol`
  checked against a hand-computed value, window-boundary bookkeeping
  checked on a small synthetic series, an out-of-range-history edge
  case, a Bates-collapses-to-Heston-through-the-backtest-harness
  consistency check (mirroring the exact-collapse property from Phase
  14), and a `[.manual][real_data]` test running the actual full
  historical backtest

**A scope decision worth being explicit about — why this backtests
hedging, not a trading strategy:** the two are easy to conflate under
"backtesting," but they answer different questions. A strategy backtest
(P&L of an actual position-taking rule) depends heavily on choices --
sizing, entry/exit, transaction costs -- that have nothing to do with
whether the *pricing model* is correct, and it's a well-known way
backtests end up telling their author what they want to hear (tune the
strategy until the number looks good). This backtest instead asks a
narrower, more defensible question: given a model's calibrated
parameters, how well would its predicted hedge ratio have tracked *real*
BTC price moves, day by day, across the model's actual (not
model-assumed) full trading history? That's a direct extension of this
project's existing validation story (implied vol matched to a real
exchange, Heston calibrated to a real chain) into the time dimension,
rather than a new kind of claim about trading skill.

**Why there's no historical options data anywhere in this backtest:**
real historical implied-vol surfaces, day by day, going back years,
aren't available for free -- checked and ruled out before defaulting to
this design. The backtest instead uses the ONLY real historical input
that's freely available (spot price) and estimates each option's writing
volatility from trailing realized volatility (30 trading days by
default) at the moment it's written -- a standard, honestly-labeled
proxy, not a claim that this is what implied vol actually was on any
given day.

**How each model's parameters are set, and why they don't get refit
every window:** Black-Scholes takes the trailing realized vol directly
as sigma. Heston and Bates hold their structural shape parameters fixed
across the *entire* ~16-year backtest -- Heston's kappa/theta/xi/rho are
literally `tools/calibrate_heston.cpp`'s real output from fitting the
committed Deribit snapshot (the same numbers the live demo's Heston card
uses), and Bates' jump_intensity/jump_mean/jump_vol are estimated once,
directly from this same price history, via a simple 4-standard-deviation
threshold on daily log returns (days whose return sits more than 4
sample standard deviations from the mean are flagged as jump days; their
mean and standard deviation become jump_mean and jump_vol; count-per-year
becomes jump_intensity). Only each window's v0 (instantaneous variance)
updates window to window, from that window's own trailing realized
variance. This is a deliberate choice, not a missed opportunity to
refit: there's no historical implied-vol surface to refit kappa/xi/rho
against day by day (that's exactly the data this project doesn't have),
so the honest thing to test is whether ONE real, already-validated
parameter set generalizes prospectively across BTC's whole price
history -- not to quietly re-fit new parameters every window and call
the result "the model."

**Why finite differences, not analytical Greeks, for Heston/Bates
delta:** this engine's Heston and Bates pricers are COS-method
characteristic-function pricers (Fang & Oosterlee, 2008) -- fast to
evaluate but with no closed-form Greeks in this codebase (Black-Scholes
gets analytical Greeks via `greeks.cpp`; Heston/Bates never did, since
nothing before this phase needed them). A central finite difference
(bump spot by 0.1%, reprice, divide) is the pragmatic choice here: cheap
enough to run at backtest scale (195 windows x 30 daily rebalances x 2
reprices per delta x 2 models is under a second), accurate enough for a
hedge ratio, and it doesn't require deriving a new closed-form
sensitivity for a characteristic-function pricer -- a real project of
its own that this phase didn't need to take on to answer the question at
hand.

**Results (full history, 2010-07-18 to 2026-09-15, 195 non-overlapping
30-day windows):**

| Model | Mean hedge P&L | Std dev of hedge P&L |
|---|---|---|
| Black-Scholes | +103.5 bps of spot | 459.1 bps |
| Heston | +73.1 bps of spot | 446.2 bps |
| Bates | +225.8 bps of spot | 430.7 bps |

**Sensitivity check — same backtest restricted to 2016-01-01 onward**
(the era Deribit itself, and therefore every other real-market number in
this project, actually covers -- 2010-2013 BTC traded on razor-thin
volume and moved on liquidity events a modern options market never
would), 129 windows:

| Model | Mean hedge P&L | Std dev of hedge P&L |
|---|---|---|
| Black-Scholes | +28.7 bps of spot | 299.8 bps |
| Heston | +4.6 bps of spot | 295.9 bps |
| Bates | +173.8 bps of spot | 283.7 bps |

The ranking (Bates highest mean, lowest variance; Heston lowest mean) is
identical in both windows, not just an artifact of the extreme early
years -- some evidence this is a real property of the data rather than a
fluke of one noisy period. The honest reading: Bates' jump term makes it
write meaningfully higher-premium options (it prices in tail risk the
other two miss entirely), and on a price history that genuinely contains
large, sudden moves — BTC's whole history, not a cherry-picked
crash — collecting that extra premium tends to pay for itself, both in
higher average P&L and in lower variance of that P&L. This is exactly
the kind of behavior a jump model is supposed to produce; seeing it show
up on real data, not just in the synthetic collapse tests from Phase 14,
is the actual point of this backtest existing.

**Honest limitations:**
- No historical options data anywhere in this backtest -- trailing
  realized vol is a proxy for what an option would have been written at,
  not a record of what one actually was. A backtest against real
  historical implied vol, if that data becomes available, would be a
  stronger and different test.
- Heston/Bates structural parameters come from ONE real calibration (a
  single Deribit snapshot, one date), applied prospectively across 16
  years of history that mostly predates that snapshot -- this measures
  whether that one calibrated shape generalizes, not whether the model
  would have calibrated well at every point along the way.
- The jump-parameter estimation (4-sigma threshold) is a simple,
  transparent heuristic, not a rigorous jump-detection method (e.g.
  Lee-Mykland); a different threshold would give different jump
  parameters and could plausibly shift the exact bps numbers above, even
  if the qualitative ranking is unlikely to flip given how consistent it
  is across both the full-history and 2016+ windows.
- The price history itself comes from a public GitHub dataset, not a
  paid financial data vendor -- spot-checked against known price history
  (above) but not independently reconciled against exchange-level data
  bar by bar.
- Options are always calls, always struck exactly at-the-money at
  writing, on a fixed 30-day tenor -- a real hedging book runs a mix of
  strikes, tenors, and both calls and puts; this is a controlled,
  single-instrument test of relative model behavior, not a claim about
  what any real book's P&L would have been.

**Defend this:**
Can explain why this is a pricing/hedging backtest and not a trading
strategy backtest, and why that's the right scope for this project (see
above). Can explain exactly what a "window" is and how its three
model-specific option-writing volatilities are each derived (BS: trailing
realized vol directly; Heston/Bates: trailing realized variance as v0,
with everything else fixed structurally). Can explain why Heston/Bates
delta uses finite differences and what the tradeoff is (fast to
implement correctly, no need for a new closed-form Greeks derivation;
costs a small amount of numerical accuracy and roughly 2x the pricing
calls per rebalance versus an analytical formula). Can name, unprompted,
every honest limitation above -- in particular that "hedge P&L" here
compares a model's OWN theoretical price against its OWN hedge cost, not
against a market-observed premium, since no historical market premium
exists to compare against.

## Phase 16 — Portfolio-level risk: VaR and stress testing (done)

**What was built:**
- `portfolio.hpp`/`portfolio.cpp`: a generic `Portfolio` (a list of
  signed `Position`s, each either an option or a position in the
  underlying), `portfolio_value()` (sums quantity x price across every
  position, model-agnostic via the same `PortfolioPricer` function-object
  pattern `backtest.hpp` used for `PriceDeltaFn`), `portfolio_delta()`
  (a single finite-difference bump on the WHOLE portfolio's value, not a
  sum of per-position analytical deltas -- works identically whether the
  book is priced under Black-Scholes, Heston, or Bates), and
  `decay_time()` (returns a copy of a portfolio with every option's
  time-to-expiry reduced, floored at 0, for horizon-based P&L
  computations)
- `make_example_market_maker_book()`: a worked example shipped alongside
  the API, the same pattern `exotic_options.hpp` used -- short 5x 30-day
  ATM straddle + short 5x 60-day 10%-OTM strangle (a market maker who
  sold that flow to clients), partially offset with +0.5 BTC of spot.
  The unhedged combination nets to about -1.06 BTC of delta (computed,
  not hand-derived); +0.5 BTC brings it to about -0.56 BTC -- partially,
  not fully, hedged, so the book has real gamma risk left over for
  VaR/stress testing to actually measure
- `var.hpp`/`var.cpp`: `historical_var()` (replays the portfolio under
  every real overlapping horizon-day return in the committed BTC price
  history), `monte_carlo_var()` (simulates forward under the Bates SDE --
  pass `jump_intensity=0` for pure Heston), `delta_normal_var()` (the
  classic linear approximation, included as a labeled baseline), and
  `run_stress_test()` + `default_stress_scenarios()` (replays the
  portfolio through 3 real historical crash windows)
- 8 new tests (`tests/test_main.cpp`, `[portfolio]`/`[var]` tags):
  portfolio value/delta sanity checks against Black-Scholes' own
  analytical delta, `decay_time`'s floor-at-zero edge case, a direct
  check that the example book is partially (not fully) hedged, an exact
  reproduction of a real historical crash's spot move and elapsed days,
  a jumps-increase-99%-VaR property test, and a `[.manual][real_data]`
  test running the full real headline computation

**Why two independent VaR methods, cross-validated against each other,
instead of picking one:** this is the same discipline every prior
model-adding phase in this project used -- Heston got a COS closed form
*and* an independent SDE Monte Carlo check (Phase 9); Bates got the same
treatment (Phase 14). VaR is the same kind of claim ("this book's tail
risk is X"), so it gets the same treatment: `historical_var()` resamples
real BTC returns with zero distributional assumption, and
`monte_carlo_var()` simulates forward under the calibrated Bates SDE
with zero dependence on which specific historical days happened to
occur. They share no code and make different assumptions, so agreement
between them (or a well-understood gap) is real evidence rather than the
same computation checked against itself.

**Why delta-normal VaR is a labeled baseline, not a third real method:**
delta-normal VaR is linear in `portfolio_delta()` alone -- it has no way
to see gamma, vega, or any other real nonlinearity in an options book's
P&L. On the example book (short gamma by construction -- a market maker
who sold volatility), that blind spot isn't hypothetical:

| Method | 95% VaR | 99% VaR |
|---|---|---|
| Historical simulation | $9,921 | $37,757 |
| Monte Carlo (Bates) | $6,826 | $16,295 |
| Delta-normal | $2,795 | $3,954 |

At 99% confidence, delta-normal understates the loss both repricing-
based methods find by roughly 4-10x. This is exactly the textbook
failure mode of delta-normal VaR on a short-options book, reproduced on
real numbers rather than asserted: a large move hurts a short-gamma book
MORE than its current delta predicts (that's what negative gamma means),
and a method that only ever looks at delta can't see that getting worse
as the move gets bigger.

**An honest finding worth being explicit about — historical and Monte
Carlo VaR don't agree closely at 99%, and here's the likely reason:**
historical VaR's 99% number (\$37,757) is more than double Monte Carlo's
(\$16,295), even though both stayed within the loose 10x cross-validation
bound the tests assert. The likely explanation: `historical_var()`
resamples from FULL 16-year history, including 2010-2013's extremely
thin, extremely volatile early market (the same era flagged as a
limitation in Phase 15) -- its tail is fatter than what a single day's
Bates simulation, calibrated to one real (and comparatively calm)
Deribit snapshot, produces. This isn't presented as a bug in either
method; it's a real, informative disagreement between "what actually
happened over BTC's whole history" and "what one particular calibrated
model predicts going forward" -- precisely the kind of gap cross-
validating two independent methods is supposed to surface, not paper
over by only running one of them.

**Stress test results** (the example book, 3 real historical crash
windows):

| Scenario | Spot move | Portfolio P&L |
|---|---|---|
| 2020 COVID crash (Feb 19 - Mar 13, 2020) | -43.1% | -$185,111 |
| 2022 FTX collapse (Nov 6-21, 2022) | -24.6% | -$51,510 |
| 2021-2022 bear market, peak to trough (Nov 10, 2021 - Nov 21, 2022) | -76.5% | -$454,882 |

All three are real losses, not simulated ones -- `run_stress_test()`
looks up the window's actual start/end spot directly in the committed
price history and decrements every option's time-to-expiry by the
actual number of calendar days elapsed, rather than treating a
multi-month bear market as an instantaneous shock.

**Honest limitations:**
- The example book, its strikes, and its hedge ratio are illustrative --
  a real market-making book runs far more positions across many more
  strikes/expiries/underlyings, with real transaction costs and margin
  constraints this doesn't model.
- `monte_carlo_var()`'s quality is entirely bounded by the Bates
  calibration it's handed -- the same structural-parameter caveats from
  Phase 14/15 (one real snapshot, fixed structural shape) apply here too.
- The historical/Monte Carlo VaR gap at 99% (above) is flagged and
  explained but not resolved -- resolving it properly would mean
  re-examining whether the full 16-year history or the 2016+ subset is
  the more appropriate resampling window for VaR specifically, a
  decision this phase deliberately left open rather than picking one
  silently.
- `delta_normal_var()`'s z-scores are hardcoded for exactly the two
  confidence levels this project uses (95%/99%), not a general
  inverse-normal-CDF solver -- a real, stated scope limit, not an
  oversight.
- Stress scenarios shock spot only (holding volatility/model parameters
  fixed at today's level) -- a real crash typically also blows out
  implied vol, which would make a short-vega book's real stress loss
  worse than what's reported here.

**Defend this:**
Can explain why two structurally different VaR methods matter more than
one precise-looking number, and can point to the real disagreement
between them (above) as evidence this project reports what it finds
rather than picking the number that looks best. Can explain exactly why
delta-normal VaR fails on a short-gamma book: it's a first-order
(linear) approximation to a P&L function that has real, negative
curvature for this specific position, so the approximation gets worse
exactly where it matters most (large moves, which is what the 99% tail
actually is). Can explain why the example book is deliberately
partially, not fully, hedged, and can compute (not just quote) its
actual delta. Can name every honest limitation above without being
asked.

## Phase 17 — Wiring Bates/backtest/VaR into the live WASM demo (done)

**What was built:** Phases 14-16 (Bates, the historical backtest,
portfolio VaR/stress testing) existed only as C++ engine code and
tests until this phase — the live demo (`docs/index.html`) still only
showed Black-Scholes/tree/Monte Carlo/Heston/exotics. This phase wires
all three into it as three new live cards, using the same
"`wasm_bridge.cpp` forwards straight into the real engine code, no
logic reimplemented in JS" discipline the existing cards already
established — extended here, not abandoned, since it would have been
the easy shortcut for the backtest/VaR cards specifically (see below).

- `wasm_bridge.cpp`: 10 new exported functions —
  `bridge_bates_price` (thin wrapper, same pattern as
  `bridge_heston_price`); `bridge_load_price_history` (parses a CSV
  *string* the browser fetched over HTTP, populates a module-global
  `std::vector<PricePoint>` reused by every later call);
  `bridge_run_backtest`/`bridge_backtest_window_count`/
  `bridge_backtest_window_pnl_bps` (runs the real `run_backtest()` for
  a chosen model over the loaded history, stashes the per-window
  results for the sparkline chart); `bridge_run_var` (real
  `historical_var()` + `monte_carlo_var()` + `delta_normal_var()` on
  `make_example_market_maker_book()`); `bridge_run_stress_test` (real
  `run_stress_test()` against one of `default_stress_scenarios()`).
  Every function that can throw (`load_price_history_from_string`,
  `run_backtest`, `run_stress_test`) is wrapped in try/catch and
  returns a status code instead of propagating — needed because the
  bridge otherwise runs the same code paths the native test suite
  already exercises, just called from JS instead of Catch2.
- `price_history.hpp`/`.cpp`: added `load_price_history_from_string()`,
  sharing its actual parsing logic with the existing file-based
  `load_price_history()` via a small internal `parse_price_history()`
  helper — a WASM module has no filesystem to open
  `data/btc_price_history.csv` from, so the browser fetches the CSV as
  *text* and hands it to this function, keeping the parsing itself in
  real C++ rather than reimplementing a CSV reader in JavaScript. Two
  new tests (`[price_history]` tag) check it parses the same format the
  file-based loader does and throws on the same malformed input.
- `docs/index.html`: three new cards — Bates (Heston-vs-jumps
  comparison, editable jump parameters, live zero-jump-collapse check),
  a historical backtest (window-length picker, runs all three models
  live over the real 16-year history, bar chart of mean hedge P&L +
  sparkline of Bates' per-window P&L), and portfolio VaR/stress testing
  (bar chart of historical/Monte Carlo/delta-normal VaR at 95%/99%,
  plus a scenario picker replaying the three real crashes). The
  backtest/VaR cards' fixed Heston/Bates structural parameters come
  from the same `heston_calibration.json` CI already regenerates every
  push (extended to also read `theta`, previously fetched but unused by
  the existing Heston card, which derives its own theta from the
  calculator's flat-vol input instead) — not a second hand-typed copy.
- `.github/workflows/pages.yml`: added the 4 new source files to the
  `em++` build's source list and 10 new function names to
  `EXPORTED_FUNCTIONS`; added `-fexceptions` (the bridge's try/catch
  blocks above do nothing without it — Emscripten's default build
  aborts the whole module on an uncaught-by-default C++ exception
  rather than letting compiled code catch its own); added a
  `cp data/btc_price_history.csv docs/btc_price_history.csv` build step
  and added that path to the "commit rebuilt assets if changed" step,
  so the demo's copy stays in sync if the source data ever changes.

**Why live in-browser instead of precomputing the numbers with CI (the
same pattern `stats.json`/`heston_calibration.json` already use):**
that would have been less work, and was seriously considered — but it
would mean the backtest/VaR cards show a screenshot, not a
computation. Measured directly before choosing: the real backtest over
the full 16-year, 195-window history (all three models) runs in
1.5-2 seconds in the compiled WASM module (Node.js, close to what a
browser's own V8 does); VaR (20,000 Monte Carlo paths + the full
historical resample) in well under 100ms. Both are fast enough to run
synchronously on a button click without a Web Worker — so there was no
real performance reason to fall back to precomputed numbers, and doing
it live is strictly more convincing evidence: a visitor isn't shown a
number, they watch their own browser derive it from the same 16 years
of real price history and the same calibrated parameters the README
quotes.

**Verification:** every new bridge function was checked against this
project's own already-trusted numbers before touching the HTML —
`bridge_bates_price` at `jump_intensity=0` against `bridge_heston_price`
(exact match, the same identity `tests/test_main.cpp`'s `[bates]` suite
checks); `bridge_run_backtest` for all three models against
`tools/run_backtest.cpp`'s native output (BS 103.5059/459.1305bps,
Heston 73.1293/446.1977bps, Bates 225.7602/430.7260bps, 195 windows —
matched to 4 decimal places, run from Node against the actual compiled
`.wasm`); `bridge_run_stress_test` against the Phase 16 stress table
above (all three scenarios' P&L matched to the cent). Beyond the
numbers, the finished page was driven end-to-end with Playwright
against a local static server — engine load, all three new cards'
buttons/pickers, real network fetches of `btc_price_history.csv` —
with zero console/page errors traceable to this phase's code (the only
console error seen was Google Fonts being unreachable from this
sandbox's own network proxy, unrelated to and pre-existing before this
change).

**Honest limitations:**
- Monte Carlo VaR's exact numbers can differ slightly between this
  sandbox's native build and the WASM build (~0.2% on the 99% figure in
  testing) — `std::mt19937` itself is a fully specified, portable
  algorithm, but `std::normal_distribution`'s exact draw-to-draw
  behavior is implementation-defined, and Emscripten's C++ standard
  library isn't bit-identical to Linux's libstdc++. Not a bug, and not
  visible in the deterministic historical VaR or the stress tests
  (neither touches a random number generator).
- Running a fresh Bates backtest over all 195 windows blocks the
  browser's main thread for roughly a second — noticeable, not
  disruptive, and the run button visibly disables/relabels itself for
  the duration. A Web Worker would remove that pause entirely; not done
  here since the pause is well under the threshold where it reads as
  the page hanging, and moving the compiled module into a worker is
  real added complexity for a demo page.
- The backtest/VaR cards' Bates jump parameters are fixed (matching
  `tools/run_backtest.cpp`'s own hardcoded values) rather than editable
  — unlike the dedicated Bates card above, which does expose them.
  They're a one-time heuristic estimate (see Phase 15), not something
  CI regenerates, so there's no live "real calibrated value" to refetch
  the way kappa/theta/xi/rho are.

**Defend this:** can explain exactly why the backtest/VaR cards fetch a
CSV as text and parse it with a second entry point into the same C++
parser, rather than either shipping a second JS CSV reader or trying to
give a WASM module real filesystem access to `data/`. Can point to the
native-vs-WASM number comparison above as the actual verification this
phase relied on, not just "it looked right in the browser." Can explain
the MC VaR discrepancy's real cause (`normal_distribution` being
implementation-defined) without hand-waving it as generic
floating-point noise.

## Phase 18 — P&L attribution and basic CVA, in the engine and the live demo (done)

**What was built:** two new risk modules, reusing the existing
`Portfolio`/`PortfolioPricer` machinery from Phase 16 rather than
inventing new plumbing:

- **P&L attribution** (`include/deriv-engine/pnl_attribution.hpp`,
  `src/pnl_attribution.cpp`): `portfolio_gamma`/`portfolio_vega`/
  `portfolio_theta`, the same bump-and-reprice pattern
  `portfolio_delta()` already used, extended to the other three Greeks;
  and `attribute_pnl()`, which Taylor-decomposes a portfolio's realized
  P&L between two market snapshots into delta/gamma/vega/theta
  contributions plus an unexplained residual — the "P&L explain"
  process a sell-side desk runs every day to check a hedge is behaving
  the way its Greeks predict, not trusting them blindly.
- **Basic (unilateral) CVA** (`include/deriv-engine/cva.hpp`,
  `src/cva.cpp`): `compute_cva()` simulates forward exposure under the
  calibrated Bates SDE (a fourth independent implementation of the SDE
  stepping loop — see `var.cpp`'s own header comment for why this
  project doesn't share that code between files), reprices the
  portfolio at each future sample point on each path, and prices in
  counterparty default risk as `(1 - recovery) * sum[EE(t) * discount(t)
  * marginal default probability]` under a flat hazard rate.
  `hazard_rate_from_cds_spread()` converts a CDS spread to that hazard
  rate via the standard credit-triangle approximation.

Both were then wired into the live WASM demo as two new cards (sections
11-12), reusing the Phase 16 example book for P&L attribution (the same
short-gamma book the VaR/stress-test card already displays, so the new
card reads as "explaining a move on the book already on screen") and a
single long 90-day at-the-money call for CVA (deliberately NOT the
mostly-short example book, since CVA only cares about positive
exposure — a short-heavy book would mostly demonstrate the trivial
zero-exposure case). Even the CDS-spread-to-hazard-rate conversion is
its own bridge function rather than a JS one-liner, keeping "every
number on this page comes from the compiled engine" literally true.

**Why this shape and not a bigger scope:** both are named, specific
gaps identified when weighing this project against what a bank-side
(not buy-side) quant role actually looks for — P&L explain and CVA/XVA
are bread-and-butter sell-side desk work that a pure pricing-methods
portfolio doesn't otherwise demonstrate. A production CVA calculation
also needs wrong-way risk, netting/collateral (CSA) modeling, a real
bootstrapped credit curve, and the rest of the XVA family (DVA/FVA/
MVA/KVA) — each a multi-person desk's actual job, not a solo portfolio
addition. `cva.hpp`'s header comment lists exactly what's out of scope
rather than letting the simplification pass silently.

**Verification:** the C++ modules were tested first, independently of
the demo — 12 new test cases / 25 new assertions in `tests/test_main.cpp`,
including: an accounting-identity test that `attribute_pnl`'s five
components sum EXACTLY to total P&L (not approximately — a refactor
bug would show up as a hard failure, not a rounding difference); a
test that portfolio-level gamma/vega/theta (finite-difference) agree
with `black_scholes_greeks`' analytical values on a single option; a
test that `total_pnl` is bit-for-bit equal to `theta_pnl` for pure time
decay (dS=0, dVol=0), since `attribute_pnl`'s decay-and-reprice step is
in that case identical to what `portfolio_theta()` already computed
internally; and for CVA, exact-zero checks at recovery=100%, hazard=0,
and for a portfolio with only negative exposure (a short option),
plus a monotonicity check (higher hazard rate -> higher CVA) and a
"CVA can't exceed the position's own value" sanity bound. Full suite:
304 assertions / 88 test cases, up from 279/76.

The WASM build was then verified the same way Phase 17 was: `em++`
compiled locally with the exact flags `pages.yml` uses, the two new
bridge functions (`bridge_run_pnl_attribution`, `bridge_run_cva`,
`bridge_hazard_rate_from_cds_spread`) checked from Node against a
native C++ program computing the identical inputs (P&L attribution
matched to 6 decimal places — no RNG involved; CVA's Monte Carlo number
differed by ~1.4%, the same native-vs-WASM `std::normal_distribution`/
`std::poisson_distribution` implementation-defined discrepancy already
documented in Phase 17, not a new issue), and finally driven end-to-end
with Playwright against a local static server serving the real
`docs/` directory — both new cards' buttons clicked, real numbers
rendered, a zero-CDS-spread run confirmed CVA collapses to exactly
$0.00 live in the browser, and the pre-existing backtest/VaR cards
re-checked to confirm nothing regressed. Zero console/page errors
traceable to this phase (only the same pre-existing, unrelated Google
Fonts proxy block noted in Phase 17).

**Honest limitations:**
- CVA's number of time steps (12, roughly monthly resolution over the
  horizon) is fixed in the bridge function rather than exposed in the
  UI — deliberately, to keep the card's inputs to what actually changes
  the answer in an interesting way, but it means a user can't explore
  the sensitivity to sampling resolution from the demo itself.
- The CVA card prices a single fixed instrument (a long 90-day
  at-the-money call) rather than letting the user build an arbitrary
  book the way the VaR/stress cards' `make_example_market_maker_book`
  does — a deliberate scope choice (see above: a short-heavy book makes
  a worse demo of what CVA actually measures), but it means the card
  can't show CVA on a mixed long/short portfolio without a code change.
- Neither module is wrong-way-risk aware: CVA's exposure simulation and
  the counterparty's default probability are modeled as independent,
  which understates risk for a counterparty whose distress correlates
  with the same moves that increase exposure to them (see `cva.hpp`).

**Defend this:** can explain why P&L attribution's Greeks are computed
once, at the "before" market, rather than at the "after" market or
some average of the two — because that's what a desk actually does
(yesterday's Greeks predict today's P&L, not Greeks re-derived with the
benefit of hindsight). Can explain why CVA only counts POSITIVE
exposure (`max(value, 0)`) rather than raw signed portfolio value — a
counterparty's default only costs you money on the trades where they
owed you, not the ones where you owed them. Can explain the credit-
triangle approximation's own limitation (a flat single-point hazard
rate, not a bootstrapped term structure) without being asked twice.

## Phase 19 — Live Deribit data: real spot/vol autofill and a real order-book comparison (done)

**What was built:** the demo's calculators previously only ever priced
whatever numbers a user typed in — no connection to what BTC or its
options were actually doing right now. This phase adds two features
that call Deribit's public REST API (`https://www.deribit.com/api/v2/
public/...`) directly from the browser, client-side, with no backend
or API key:

- **Live autofill on the main calculator (section 06):** a "Fetch live
  BTC data" button calls `get_index_price` for the real current spot
  and `get_book_summary_by_currency` for the full live BTC options
  chain, picks the live contract closest to at-the-money and ~90 days
  out that has a usable `mark_iv`, and fills the Spot and Volatility
  fields with real numbers instead of requiring the user to look them
  up and type them in.
- **Live market check (new section 13):** given a target strike, days
  to expiry, and type, fetches the same live chain, finds the closest
  real, currently-tradeable contract, converts its BTC-denominated
  mark/bid/ask into USD using the real live index price, and reprices
  that exact contract with this engine's own calibrated Heston model
  using the *market's own* implied vol as the input — displaying the
  real market price and this engine's model price side by side with
  the percentage gap between them. This is the actual "does my model
  agree with the real market" check, not a hypothetical one built from
  assumed inputs.

**Why this shape:** before writing any code, CORS reachability was
verified for real rather than assumed — a live cross-origin `fetch()`
was run from an actual Safari session (not this sandbox, which is
itself proxied and cannot reach `deribit.com`) against
`get_index_price`, and the real response headers were inspected
(`Access-Control-Allow-Methods: GET, POST, OPTIONS` present, status
200, answered despite `Origin: https://claude.ai`) to confirm Deribit's
public endpoints are genuinely open to arbitrary origins, not just
Deribit's own site. Only after that confirmation was the feature built.
`get_book_summary_by_currency` was chosen over calling `ticker` once
per instrument because it returns bid/ask/mark price and mark IV for
every live instrument in a single request — the chain has hundreds of
strikes/expiries at any time, so one bulk call avoids hundreds of
round trips. Both features share one instrument-name parser
(`BTC-27DEC26-90000-C` → strike/expiry/type) and one spot+chain fetch
helper rather than duplicating that logic, since — unlike the
project's deliberately-independent SDE stepping loops (see `var.cpp`'s
header comment) — there's no cross-validation value in parsing the
same string format two different ways.

**Verification:** both inline scripts were syntax-checked
(`new Function(source)` on each extracted `<script>` block), then
driven end-to-end with Playwright against a local static server
serving the real `docs/` directory. Because this sandbox's own network
proxy blocks arbitrary outbound hosts (confirmed via the proxy's own
status endpoint — `deribit.com` is not on the allowlist), the live
success path could not be exercised from here; what was verified
instead is that: the pre-existing calculator (Black-Scholes/tree/MC/
Heston, sections 01-12) still produces identical output with zero
regressions after these changes: BS $12,392.30, Heston $11,950.52 on
the default inputs, unchanged from before this phase; both new
features fail closed and legibly (`Failed to fetch` surfaces as a
plain-English error message in the UI, not a hang, a crash, or a
silent fallback to fake numbers); both buttons correctly re-enable
themselves after a failed request; and zero uncaught JS exceptions
were raised (`page.on('pageerror', ...)` captured none — only expected
network-level console errors from the sandbox's own blocked egress).
The actual live-success path (a real 200 response with real BTC data)
was confirmed manually and interactively in a real browser session
before this was built, the same way the CORS-open finding above was
confirmed — not inferred from documentation alone.

**Honest limitations:**
- These two features are the only parts of the entire live demo with
  a runtime dependency on a third party being reachable. Every other
  card (sections 01-12) is fully self-contained, computing everything
  from the compiled WASM engine and the committed CSV/JSON data —
  these two will show a plain error message instead of a price if
  Deribit is down, rate-limits the request, or is blocked by the
  visitor's own network/ad-blocker/VPN.
- The live-autofill's "closest to at-the-money and ~90 days out"
  selection is a simple weighted-distance heuristic
  (`moneyness * 3 + tenorGap`), not a real liquidity or open-interest
  filter — it can pick a thinly-traded, wide-spread contract's IV over
  a more liquid one nearby.
- No caching or rate-limiting on repeated clicks — each click issues a
  fresh pair of requests. Fine for demo usage by one visitor at a time,
  not written for high click-through traffic.
- `mark_iv` is used as a single flat number for both the Heston
  `v0` and `theta` inputs in section 13's reprice, same simplification
  the main calculator already makes (see Phase 12's honest limitations)
  — the real market's actual vol-of-vol and mean-reversion aren't
  re-derived from the live chain, only the calibrated defaults are
  reused with the live IV swapped in for the flat-vol level.

**Defend this:** can explain why CORS was tested with a real browser
request and inspected response headers before writing a line of
fetch code, rather than assumed from the fact that "public APIs are
usually CORS-open" — an assumption that would have shipped a feature
with no way to know if it actually worked. Can explain why
`get_book_summary_by_currency` (one bulk call) was chosen over
per-instrument `ticker` calls (one call per contract), and the
BTC-to-USD conversion Deribit's own quoting convention requires that a
naive implementation would silently get wrong by two orders of
magnitude. Can explain exactly what "model vs. market gap" in section
13 does and doesn't prove: agreement confirms the calibrated
kappa/xi/rho still describe the real smile reasonably well when fed
today's real IV; disagreement doesn't necessarily mean the model is
wrong — it can just as easily mean the specific live contract matched
is thin, or that IV alone doesn't capture the smile's actual shape at
that strike the way a full recalibration against the live chain would.

## Phase 20 — Consolidating the live demo: tabbed calculators, one live-data fetch for the whole page (done)

**What was built:** the live demo had grown to thirteen separate
sections, several of which were really the same underlying idea shown
three different ways (vanilla/tree/MC/Heston, the two exotics, and
Bates were three separate full-page cards for "price a European-style
option"; VaR/stress-testing and P&L attribution were two separate
cards built on the same example short-gamma book). This phase merges
them into two tab-switched cards instead of removing anything: section
06 ("Try it live: pricing methods") now holds three pill-selected tabs
— Vanilla/Tree/MC/Heston, Exotic (Asian & Barrier), and Bates — and
what is now section 08 ("Portfolio risk: VaR, stress testing & P&L
attribution") holds two tabs — VaR/stress and P&L attribution. Every
calculator's own inputs, buttons, and output panel are unchanged
inside its tab; only the wrapper and a pill-row switcher were added.
The remaining, genuinely distinct cards were renumbered to fill the
gaps: Historical backtest is now section 07, Basic CVA is section 09,
and Live market check is section 10 — same content, no changes beyond
the number in the section head.

The other half of this phase: the "Fetch live BTC data" button, which
Phase 19 only wired into the vanilla pricer's Spot/Volatility fields,
now fills Spot and Volatility on every calculator on the page in one
click — the vanilla, exotic, and Bates tabs, plus VaR/stress, P&L
attribution, and CVA — instead of six separate live-fetch buttons or
requiring the same lookup to be repeated per card. The button, its
Deribit fetch, and its "closest to at-the-money, ~90 days out"
selection logic from Phase 19 are otherwise unchanged; only the fill
step at the end was widened from two field writes to twelve.

**Why this shape:** the user's own framing was that having "so many
different calculators" that don't share data felt disjointed —
several were variations on one question (price this option under
model X) rather than genuinely separate tools, and typing the same
spot/vol into six boxes by hand undercuts the point of a "live data"
feature. Tabs were chosen over deleting any pricer or flattening them
into one mega-form: every existing model's validation story (tree vs.
closed-form convergence, Heston vs. flat-vol, Bates vs. Heston) still
needs its own inputs and its own output panel to make sense on its
own, so collapsing them into one shared form would either hide inputs
some models don't use or silently reuse a value (like Bates' jump
parameters) in a context where it's meaningless. Reusing the existing
`.pill-row`/`.pill` CSS already used for the window-length and
stress-scenario pickers, rather than inventing new tab styling, keeps
the page visually consistent and meant zero new CSS was needed — only
a `data-tab` attribute on each button and a small `wireTabPills()`
helper that shows the matching tab and hides its siblings. Each
existing `calc-card` div was wrapped with a unique `id` and toggled via
`style.display`, with its internal HTML, IDs, and JS handlers left
completely untouched — so the price/exotic/Bates/VaR/P&L computation
functions needed zero changes; only the DOM wrapper and one click
handler were added. This was chosen over rewriting the calculation
logic into a single unified multi-mode function, which would have
touched far more surface area for the same visible result and risked
introducing a real numerical bug in code that was already correct and
tested.

For the live-data fill, filling every card's fields from one fetch —
rather than adding five more per-card fetch buttons — was the direct
ask ("have them all pull live BTC data from exchange"): one Deribit
lookup already returns the one number (spot) and one representative
implied vol that every card's Spot/Volatility inputs want, so there's
no reason to hit the API six times or make the user repeat the click.
Cards whose reprice needs more than Spot/Vol to run (VaR, P&L, CVA,
and the backtest all require a "run" click because they run Monte
Carlo simulations that take a moment) still require their own button
press after the fields are filled — auto-triggering a multi-second
simulation the instant a fetch completes, without the user asking for
it, would be a worse experience than a filled-in form waiting for
"Run".

**Verification:** the full merge was applied directly to
`docs/index.html` (no intermediate build step — this is a static
page), then checked three ways. First, a scripted diff of every
element `id` in the affected region before and after confirmed zero
IDs were dropped, renamed, or duplicated across the whole file — all
164 pre-existing IDs (all six calculators' fields, buttons, and output
spans) are still present, plus the new pill-row and tab-wrapper IDs.
Second, the page was served from a local static server and driven with
Playwright: clicking each pill in both tab groups correctly shows the
target tab and hides its siblings (verified via `is_visible()` on each
tab both before and after each click, for all three pricing tabs and
both risk tabs); the vanilla pricer's "Price it" button was clicked
after the merge and produced the same numbers as before this phase
(Black-Scholes $12,392.30, Heston $11,950.52 on default inputs — an
exact match, confirming the wrapping change didn't alter any
computation); and the WASM engine reported "loaded" with zero
JavaScript exceptions (`page.on('pageerror', ...)` captured none)
across the whole interaction. Third, the live-fetch button was
clicked; as expected from this sandbox's blocked egress to Deribit
(same constraint as Phase 19), it failed with a plain "could not fetch
live data" message rather than hanging or throwing, and every field it
would have filled was left untouched at its prior value — confirming
the failure path doesn't corrupt the page. The actual live-success
fill across all six cards could not be exercised from this sandbox for
the same reason Phase 19's live-success path couldn't be: this only
the user's own browser, with real network access to Deribit, can
confirm end-to-end (open the page, click "Fetch live BTC data" once,
and check that Spot/Volatility updated on the exotic, Bates, VaR/
stress, P&L, and CVA tabs — not just the vanilla one).

**Honest limitations:**
- The live-fetch fill still only writes Spot and Volatility. It does
  not re-run any calculator automatically except the vanilla pricer
  (which was already wired that way in Phase 19) — every other tab
  still needs its own "Price it"/"Run" click after the fields update,
  by design (see above), but a user who doesn't notice the fields
  changed and doesn't press the button will see stale output sitting
  next to fresh inputs until they do.
- Bates' vol field is still separately labeled "Base volatility (%)"
  and semantically distinct from a flat Black-Scholes vol (it's the
  long-run level the jump-diffusion process reverts to, not a constant
  vol), and VaR/P&L/CVA's "Annualized vol (%)" is a plain constant-vol
  assumption for their simulations — the live fetch fills all of them
  with the same one number (a single near-the-money mark IV) as a
  simplification, the same one Phase 19 already made for the vanilla
  and Heston reprice; it does not mean these fields' inputs are now
  interchangeable or equally rigorous.
- The pill/tab mechanism is pure client-side `display:none` toggling
  with no URL state or persistence — reloading the page or sharing a
  link always resets to the first tab in each group (Vanilla and VaR),
  which matters if someone bookmarks or shares a link expecting to
  land on Bates or P&L attribution specifically.
- As before, sections 06's live-fetch and section 10's live market
  check remain the only parts of the page with a runtime dependency on
  a third party; the newly-merged tabs inside section 06 and section
  08 inherit that same dependency only insofar as they can now be
  pre-filled by it — their own core pricing/simulation logic is still
  fully self-contained and works with zero network access if the
  fields are typed in by hand.

**Defend this:** can explain why tabs (hide/show existing cards) were
chosen over either deleting pricers or merging their calculation logic
into one function — the former loses each model's independent
validation story, the latter risks a real bug in previously-correct,
tested numerical code for a change that's purely about page
organization. Can explain why the live-fetch fill writes fields
without auto-running every calculator: VaR/P&L/CVA/backtest all run
Monte Carlo simulations that take real wall-clock time, and
auto-triggering all four the instant a fetch resolves would surprise
the user with unrequested computation rather than leaving them in
control of when to run each one. Can explain exactly what the "same
one mark IV feeds six differently-shaped vol inputs" simplification
does and doesn't cost: it's a fast, real, defensible starting point
that saves six manual lookups, not a claim that Bates' jump-vol,
Heston's vol-of-vol, or a constant-vol VaR assumption are the same
statistical object — each card's own documented limitations about its
vol assumption (Phase 14, Phase 16, Phase 18) still apply exactly as
before; only the number that seeds each field changed from
hand-typed to fetched.

## Phase 21 — SVI volatility surface fitting: one real curve per live expiry (done)

**What was built:** a new pricing methodology, not a variant of one already
here. `svi.hpp`/`svi.cpp` implement Gatheral's raw SVI parameterization
(2004) — five parameters (a, b, ρ, m, σ) describing one expiry's entire
implied-total-variance smile as a curve, fit directly to real market
quotes rather than derived from a stochastic-process story the way
Heston/Bates are. `svi_calibration.hpp`/`svi_calibration.cpp` fit those
five parameters to a real expiry's market points by minimizing squared
error in total-variance space, reusing the same Nelder-Mead simplex
optimizer `heston_calibration.cpp` already had — pulled out into a new
shared `optimizers.hpp` header so both calibrators call one implementation
instead of two subtly different copies. Two real, published no-arbitrage
checks come with it: a **butterfly check** (`svi_butterfly_condition`/
`svi_check_butterfly_arbitrage_free`), which catches a fitted curve that
implies a negative risk-neutral density somewhere in its wings, and a
**calendar check** (`svi_check_calendar_arbitrage_free`), which catches
two expiries' fitted curves implying that total variance decreases with
time — a real, constructable arbitrage against the fitted surface itself,
not just two curves that happen to disagree.

Three new WASM exports (`bridge_svi_fit`, `bridge_svi_vol`,
`bridge_svi_calendar_check`) wire this into a new eleventh section of the
live demo, "Volatility surface: fitting SVI to the whole live chain." One
button fetches the real live Deribit chain (the same fetch every other
live card on this page uses), groups it by real expiry, fits one SVI
curve per expiry with enough live contracts, and reports each expiry's
fit quality and both arbitrage checks — plus a chart (a new inline SVG
smile chart, `drawSviSmile`, following the same hand-drawn-SVG convention
`drawSparkline` already established for the backtest card) plotting the
fitted curve against the real market points, switchable across every
expiry that fit via a pill row.

**Why this shape:** the user's own framing, after being walked through why
Heston/Bates are respectable but not what a real crypto or equity vol desk
actually quotes off of day to day, was to add what those desks actually
use. SVI is that answer: unlike a stochastic-vol model, which has to
compromise across an entire smile with a handful of shared parameters, SVI
fits each expiry's real, observed skew/smile shape essentially exactly,
which is why it's the industry-standard curve for surface construction
rather than a competing pricing model. One curve per expiry (not one
shared fit across all of them, the way Heston's multi-expiry calibration
in Phase 12 deliberately pools maturities) is the correct design here, not
a missed opportunity to generalize: pooling would defeat SVI's entire
purpose, which is matching each maturity's real shape exactly rather than
explaining several maturities from one process. Fitting in total-variance
space (not price space, unlike `calibrate_heston`, and not vol space
directly) is the standard convention in the SVI literature (Gatheral &
Jacquier) and is what makes both arbitrage-check formulas exact rather
than approximate. Extracting the shared `optimizers.hpp` was a small,
low-risk refactor done specifically to avoid a second copy of a
71-line simplex implementation silently drifting out of sync with the
first — verified with a full test-suite run immediately after, before any
SVI-specific code was even written, to isolate that change from everything
that came after it.

**Verification:** four new unit tests, following this project's own
established pattern for a new pricer. A synthetic-recovery test generates
a full smile from known SVI parameters and confirms `calibrate_svi`
recovers a curve matching it to within 0.01 percentage points of implied
vol at every strike (SVI's 5 parameters, unlike Heston's kappa/theta,
*are* all identifiable from one expiry, so unlike the Heston synthetic
test this one does assert the fitted curve matches, not just the fit
quality). A direct formula test hand-checks `svi_total_variance` and both
arbitrage checks against simple, independently-computable cases, including
a deliberately-broken calendar case (a lower total variance at a longer
expiry) to confirm the check actually fails when it should, not just
passes when given reasonable input. An edge-case test confirms the
calibrator refuses to fit fewer than 5 real points rather than return an
overfit, meaningless curve. A real-data test fits every liquid expiry
(same [0.7, 1.4] moneyness filter the Phase 9 Heston real-data test
already established, for the same reason) in the committed Deribit
snapshot and found a genuine, monotonic, honestly-reported pattern: fit
quality degrades from 3.74 percentage points at 1 day to 1.76 at 2 days to
0.58 at 4 days, and the butterfly check actively **fails** at 1, 3, and 4
days — a raw SVI curve's limited curvature budget cannot track how steep
and narrow the real sub-two-week BTC smile gets without implying a
negative density somewhere. Every expiry from 15 days out to 281 days,
by contrast, fits under 0.22 percentage points with both arbitrage checks
clean. The test asserts against that measured boundary (excludes expiries
under 14 days, bounds the rest at 0.35pp with real headroom above the
observed 0.21pp worst case) rather than an assumed one. The live demo
itself was driven end-to-end with Playwright against mocked (not real,
since this sandbox's network cannot reach Deribit — the same limitation
Phase 19/20 already documented) Deribit responses built from a known SVI
curve: the fit recovered the ground truth to 0.026pp, the multi-expiry
case correctly reported "all adjacent pairs passed the calendar check,"
and a deliberately-broken second scenario (a far expiry with much lower
total variance) correctly reported "1 adjacent pair(s) failed the calendar
check" and flagged the specific expiry in the results table — confirming
the failure path surfaces honestly rather than silently passing. The full
374-assertion, 102-test-case suite (up from 327/91 before this phase) and
the entire rest of the live demo (all ten prior sections, all tab
switching, the global live-fetch button) were re-verified with zero
regressions after every change in this phase, and the real WASM module was
actually compiled with `em++` in this environment (not merely
syntax-checked) to confirm the new bridge exports build and run correctly
end-to-end, not just compile.

**Honest limitations:**
- Sub-two-week expiries are a genuine, measured weak point for raw SVI on
  this real chain, not a hypothetical caveat — the live demo will show a
  failed butterfly check or an elevated error on those expiries when
  they're live, and that's reported plainly rather than filtered out of
  the live UI the way the test suite filters them out of its own pass/fail
  bound. A smoother-in-time extension (e.g. SSVI, which ties nearby
  expiries' curves together) is the standard fix for exactly this, and is
  out of scope here.
- Each expiry's curve is fit completely independently, so nothing
  structurally prevents a calendar-arbitrage failure between two real,
  individually-reasonable-looking fitted curves on a noisy real
  snapshot — the calendar check exists specifically to catch this live
  rather than assume independence is harmless, and it will report a
  failure honestly when the real chain produces one, the same way the
  butterfly check does for short-dated expiries.
- SVI fits Deribit's own reported `mark_iv` per contract directly (the
  same number every other live card on this page already reads), not an
  independently re-solved implied vol the way this project's own
  Phase 8 headline result does for the committed snapshot — reproducing
  that independent Newton-Raphson solve against a live options chain
  client-side, in the browser, for every fetch, was judged not worth the
  added complexity for what the live SVI card is actually demonstrating
  (surface *shape*, not an independent audit of Deribit's own reported
  numbers, which the headline result already covers separately against
  the frozen snapshot).
- This is a curve fit, not a pricing model: it has no formula for pricing
  an option directly the way Black-Scholes, Heston, and Bates do — its
  output is a vol at a given strike/expiry, meant to feed one of this
  project's actual pricers, not replace them. It isn't wired into the
  vanilla/Heston/Bates pricer or the Live Market Check card's reprice
  step in this phase; doing so (feeding a fitted SVI vol into Heston's
  reprice instead of one contract's flat mark IV) is a natural, small
  follow-up, not part of this phase's scope.

**Defend this:** can explain why SVI is fit per-expiry rather than pooled
across expiries the way Heston's multi-expiry calibration deliberately is
— it's not an inconsistency between the two calibrators, it's that they're
solving different problems (one process explaining several maturities at
once vs. one curve matching one maturity's real shape exactly). Can
explain why total-variance space is the correct space to fit in rather
than price space or vol space directly, and what specifically breaks in
the arbitrage-check formulas if it weren't. Can explain exactly what the
butterfly and calendar checks each catch and why both are necessary (one
catches a single curve implying nonsense, the other catches two otherwise-
valid curves implying nonsense together), and can produce the real,
measured, checked-not-assumed evidence for both: the sub-two-week failure
pattern on the real committed snapshot, and the deliberately-constructed
calendar-violation scenario the live demo was tested against. Can also
explain why the shared Nelder-Mead extraction was verified in isolation
before any SVI code depended on it, rather than trusting that a
copy-paste-and-modify refactor didn't subtly change Heston's own
calibration behavior.

### Phase 21 fix — live SVI results were showing a 100x-inflated error, plus a real outlier-quote problem on short-dated expiries

First live run of Section 11 against real Deribit data showed rmse
figures that made no sense on their face for the short-dated expiries
(1d and 2d showed 620pp and 70pp -- an implied-vol error of hundreds of
percentage points isn't a "hard fit," it's broken). Rather than accept
that as more evidence for the already-documented short-dated-expiry
limitation, it was investigated as a bug, and it was two separate real
bugs:

- **A 100x display bug.** `calibrate_svi`'s `rmse_iv_pp` is already
  expressed in percentage points (`svi_calibration.cpp` multiplies the
  vol difference by 100 itself before squaring and averaging). The live
  demo's JS multiplied it by 100 *again* before rendering it. This alone
  explains almost the entire gap: the "good" long-dated expiries' true
  rmse (89d/180d/271d/362d: ~0.10/0.04/0.03/0.04pp) matched the native
  test suite's documented worst case almost exactly once the double
  scaling was removed -- it was never actually wrong, only mislabeled by
  two orders of magnitude on screen.
- **A real outlier-quote problem, separate from display.** Deribit still
  reports an extrapolated "mark IV" for contracts with essentially no
  real market -- a handful of illiquid, deep-wing quotes near expiry can
  carry mark IVs of several hundred percent. Reproduced directly: feeding
  `calibrate_svi` a synthetic 1-day bucket with a few such wild wing
  quotes mixed into an otherwise sane smile pushed a genuine, unscaled
  rmse to 100+pp, because an unweighted least-squares fit lets those
  untradeable points dominate the objective. Fixed with a standard
  robust-statistics trim (median absolute deviation, k=4) applied to each
  expiry's live quotes before fitting, in the JS grouping step -- no
  change to `calibrate_svi` itself, since the real, liquid part of the
  chain was never the problem.

Both were confirmed with a Playwright run against a mocked chain built
specifically to contain wild wing quotes: pre-fix, a clean 15-point
1-day smile fit at rmse 400pp; post-fix (both bugs corrected), the same
bucket fits at rmse ~4pp -- elevated, consistent with the documented
short-dated difficulty, but a real, honest number, not a broken one.
Nothing in `calibrate_svi`, `svi.cpp`, or the WASM bridge changed; this
was a live-demo-only fix (`docs/index.html`), so no rebuild of
`deriv_engine.js` was needed.

## Phase 22 — Wiring the fitted SVI surface into the Live Market Check card (done)

**What:** Section 10's Heston-vs-real-market comparison used to set Heston's
`v0`/`theta` from the single matched contract's own quoted mark IV, squared
— i.e., feeding the model the market's own answer for the exact contract
being checked, off one (possibly thin, possibly noisy) quote. This phase
wires it to Section 11's SVI machinery instead: before repricing, it fits
a real SVI curve to every live, liquid contract at that same expiry (same
5-parameter fit, same 0.7-1.4x moneyness band, same median/MAD outlier
trim as section 11 -- both cards now share one `madTrim()` and a new
`fitSviCurveForExpiry()` helper, rather than duplicating the logic) and
reads the ATM variance off the fitted curve for `v0`/`theta`. A second,
independent price is now also shown: plain Black-Scholes at the fitted
curve's vol for the exact target strike -- so the card shows three real
numbers side by side for the same live contract: the real market mark,
this engine's Heston price (now smile-informed), and a pure curve-fit
price. When an expiry doesn't have enough real, liquid quotes to fit a
curve (`calibrate_svi`'s own minimum is 5 points), the card says so
plainly and falls back to the old single-contract-IV behavior for Heston,
with the SVI-implied row showing "--" rather than a number built on too
little data.

**Why:** a single contract's own IV is the noisiest possible read on "the
vol level" for its expiry -- it's one data point, possibly on a thin
strike, with no information from the rest of that expiry's real chain.
The whole point of Section 11's SVI fit is a smoothed, arbitrage-checked
read of the *entire* live expiry; Section 10 had that machinery sitting
next to it, unused, since Phase 21 shipped. This also turns "does my
model agree with the market" into a genuine three-way comparison a real
desk would ask: does a stochastic-vol process (Heston) agree with the
market, does a pure curve fit (SVI) agree with the market, and do the two
models agree with *each other* -- three different, independently
checkable claims, not one.

**Verification:** confirmed with two Playwright runs against a mocked
live chain. A well-populated synthetic 90-day expiry (14 real points
after trimming) correctly fit an SVI curve, derived a 41.4% ATM vol, and
produced a Heston price (~$6,938 for a near-ATM 90-day call at that vol
and spot) matching a hand-computed ATM-call approximation
(`spot × vol × sqrt(T / 2π)`) almost exactly -- confirms the pipeline,
not just that it produced *a* number. A deliberately thin synthetic
5-day expiry (3 contracts, below `calibrate_svi`'s 5-point minimum)
correctly triggered the documented fallback: Heston used the single
contract's own IV, and the SVI-implied row showed "--" with an honest
explanation, instead of silently fitting a meaningless curve to 3 points
or crashing.

**Honest limitations:**
- Heston's `kappa`/`xi`/`rho` still come from the separate, real
  multi-expiry calibration in `heston_calibration.json` (unchanged by
  this phase) -- only `v0`/`theta` are now expiry-specific and
  SVI-derived. A single Heston process fit across many expiries and a
  per-expiry SVI-implied variance level are two different, real things;
  this phase doesn't reconcile them into one internally-consistent model,
  it presents both as what they are.
- The SVI-implied price uses plain Black-Scholes, not Heston, at the
  fitted vol -- it's a genuinely different pricing approach from the
  Heston row, which is the point of showing both, not an oversight.
- No dividend yield adjustment beyond what already existed elsewhere in
  the engine (q=0, matching every other card's BTC-scope assumption).

**Defend this:** can explain exactly why feeding a contract its own
quoted IV back into a model and comparing against that same contract's
market price is weaker evidence than it looks (the model is partly
reproducing information it was just handed), and why sourcing `v0`/theta
from the whole expiry's fitted curve instead is a real, if partial, fix
for that. Can explain why Heston's process parameters and SVI's per-
expiry level aren't merged into one calibration, and why that's an
honest choice rather than a shortcut. Can point to the exact hand-checked
number (the 90-day ATM call) that confirms the wiring computes something
real, not just something plausible-looking.
