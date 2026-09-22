# deriv-engine

A C++20 derivatives pricing engine for Bitcoin options, built from
scratch: closed-form pricing, numerical methods, Monte Carlo
simulation, hand-built automatic differentiation, a Heston
stochastic-volatility model, and calibration against live market
data — validated against a real exchange, not just textbook values.

[![CI](https://github.com/roberto0607/deriv-engine/actions/workflows/ci.yml/badge.svg)](https://github.com/roberto0607/deriv-engine/actions)

## Live demo

**[roberto0607.github.io/deriv-engine](https://roberto0607.github.io/deriv-engine/)**
— the actual C++ (Black-Scholes, binomial tree, Monte Carlo, Heston)
compiled to WebAssembly via Emscripten and run directly in the
browser. Every number on that page comes out of the real compiled
engine, not a JS reimplementation, and the Heston card's kappa/xi/rho
are calibrated against the real Deribit chain by CI on every push
(see `tools/calibrate_heston.cpp`), not hand-typed.

## Headline result

Pulled the full live BTC options chain from Deribit (904 real
contracts) and independently solved implied volatility for every one,
using nothing but this engine's own Black-Scholes implementation and a
hand-built Newton-Raphson solver — no reference to Deribit's own
reported values during solving.

**Result: 96.8% convergence, a median gap of 0.3 percentage points
from Deribit's own reported implied volatility across the entire
chain** (mean gap 2.2pp — pulled up by a thin tail of near-expiry and
far-out-of-the-money contracts where Deribit's own quoted IV is likely
thin or stale to begin with; 90% of contracts land within 5.6pp — see
`docs/phase-log.md` for the full breakdown).

The resulting volatility smile, plotted directly from real market
prices (99-day expiry, 112 converged contracts):

![BTC options implied volatility smile, 99-day expiry](docs/vol_smile.png)

Implied vol bottoms out in the mid-30s% near the money and rises on
both wings — over 100% for the deepest out-of-the-money puts, into the
high 60s% for the deepest out-of-the-money calls — the same skew this
engine's Heston model is calibrated to reproduce.

## What's implemented

| Method | What it proves | Validated against |
|---|---|---|
| Black-Scholes closed-form | Exact European pricing + analytical Greeks | Textbook reference values, put-call parity |
| Binomial tree | American early exercise | Converges to Black-Scholes as steps→∞ |
| Monte Carlo | Simulation-based pricing, variance reduction | Converges to Black-Scholes; antithetic + control variates give a measured 6.97x variance reduction |
| Reverse-mode AAD | Hand-built automatic differentiation | Matches analytical Greeks; benchmarked against bump-and-revalue |
| Delta-hedging simulation | Realized P&L vs. theoretical price | Mean PnL ≈ 0 across 2000 simulated paths; hedging error shrinks ~4.5x from monthly to daily rebalancing |
| Longstaff-Schwartz | American option pricing via Monte Carlo | Converges to the binomial tree on two independent parameter sets |
| Crank-Nicolson PDE | Direct numerical solution of the Black-Scholes PDE | Converges to Black-Scholes; explicit-scheme instability demonstrated directly (diverges to -10^22 under the same conditions CN stays accurate) |
| Implied vol + calibration | Real market data, not synthetic inputs | See headline result above |
| Heston stochastic volatility | Closed-form pricing via the COS method (Fang & Oosterlee), capturing volatility skew a flat-vol model can't | Collapses to Black-Scholes in the deterministic-variance limit; independently cross-checked against a full Monte Carlo simulation of the Heston SDE; kappa/xi/rho calibrated (Nelder-Mead, regularized) against 66 real Deribit contracts, fitting the smile to 0.43pp RMSE in vol space |

Every method above independently arrives at the same price for the
same European option — four structurally different approaches
(formula, tree, simulation, PDE) agreeing is the core correctness
evidence for this engine.

## Notable engineering details

- **Every bug found is documented, not hidden.** Real bugs were found
  and fixed during development, each with a before/after test proving
  the fix — a division-by-zero at T=0, a double-counted premium in the
  hedging PnL calculation, a random-seed mismatch that corrupted a
  finite-difference comparison, and a catastrophic-cancellation bug in
  the Heston characteristic function that caused an 82% pricing error
  for far-out-of-the-money options. See `docs/phase-log.md`.
- **AAD was benchmarked honestly, not just claimed to be fast.**
  At the tested scale (4 Greeks, 100k paths), hand-built reverse-mode
  AAD was measured to be *slower* than naive bump-and-revalue — the
  reasoning for why, and where the actual crossover point is, is
  documented in `docs/numerics.md` rather than glossed over.
- **Deliberate scope decisions, not gaps.** Discrete dividends,
  trinomial trees, and second-order AAD (gamma) were all considered
  and explicitly deferred, with the reasoning written down — see
  `docs/numerics.md`.
- **Heston calibration is honest about non-identifiability.** kappa
  and theta aren't separately identifiable from a single expiry's
  smile — an unregularized fit converges to a valid but implausible
  parameter set. Real-data calibration uses a small L2 penalty toward
  a reasonable initial guess to resolve this; the synthetic
  parameter-recovery test deliberately doesn't, so its accuracy
  assertions stay meaningful. See `include/deriv-engine/heston_calibration.hpp`.
- **Reproducible market data.** Live data is pulled once and frozen
  to a snapshot (`data/btc_chain_snapshot.json`) rather than tests
  depending on a live API call, keeping the test suite deterministic.

## Build & test

```bash
git clone https://github.com/roberto0607/deriv-engine.git
cd deriv-engine
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build
./build/tests
```

CI runs the same sequence on a clean container on every push — see
the badge above.

## Documentation

- [`docs/design.md`](docs/design.md) — architecture and design decisions
- [`docs/numerics.md`](docs/numerics.md) — the math and numerics behind every method, including edge cases and deliberate scope decisions
- [`docs/phase-log.md`](docs/phase-log.md) — a build log across the project: what was built, how it was validated, and every bug found along the way

## Tech stack

C++20 · CMake · Catch2 · nlohmann/json · Deribit public API · Emscripten/WebAssembly

## Scope

Built around BTC options specifically (continuous dividend yield = 0,
correct for an asset with no dividends) rather than equities. See
`docs/numerics.md` for the reasoning behind this and other scope
decisions.
