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

## Phase 1 — Black-Scholes closed-form (in progress)

**What was built:**
- `MarketData` struct (spot, risk-free rate, dividend yield, volatility)
  and `EuropeanOption` struct (strike, time-to-expiry, call/put type)
- `black_scholes_price()` — closed-form pricing for European calls and
  puts, using the standard d1/d2 formulation
- Internal `norm_cdf()` helper implementing the cumulative normal
  distribution via `std::erfc`

**Validated against:**
- Standard textbook reference case (S=100, K=100, T=1yr, r=5%, vol=20%):
  call = 10.4506, put = 5.5735 — both matched within 0.1% tolerance
- Still pending: put-call parity check, and the edge-case suite
  (near-zero vol, deep ITM/OTM, near-expiry, zero/negative rates)

**Defend this:**
Can explain the formula at both levels: the intuition (a call's price is
the expected value of the asset in the exercise scenarios, minus the
discounted strike weighted by probability of exercise) and the precise
mechanics (N(d2) is the risk-neutral probability of finishing
in-the-money; d1 differs from d2 by one sigma*sqrt(T) term because it
weights the expected asset price rather than pure exercise probability,
correcting for the log-normal distribution's mean being pulled above
its median by volatility). Can also explain why higher volatility
always increases option value regardless of call/put — capped downside,
open-ended upside. Not yet done: Greeks derivation, and full whiteboard
derivation of the BS PDE itself from the replicating-portfolio argument
— that's the remaining gap before this phase is genuinely closed out.

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