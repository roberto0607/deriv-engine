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

## Phase 2 — Tree pricers

**What was built:**


**Validated against:**


**Defend this:**


---

## Phase 3 — Monte Carlo core + variance reduction

**What was built:**


**Validated against:**


**Defend this:**


---

## Phase 4 — AAD Greeks

**What was built:**


**Validated against:**


**Defend this:**


---

## Phase 5 — Delta-hedging PnL simulation

**What was built:**


**Validated against:**


**Defend this:**


---

## Phase 6 — Longstaff-Schwartz American MC

**What was built:**


**Validated against:**


**Defend this:**


---

## Phase 7 — Finite difference PDE (Crank-Nicolson)

**What was built:**


**Validated against:**


**Defend this:**


---

## Phase 8 — Calibration & real market data

**What was built:**


**Validated against:**


**Defend this:**


---

## Phase 9 (stretch) — CUDA Monte Carlo

**What was built:**


**Validated against:**


**Defend this:**