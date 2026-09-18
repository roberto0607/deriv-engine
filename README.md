# deriv-engine

A C++20 derivatives pricing engine for Bitcoin options, built from
scratch across 8 phases: closed-form pricing, numerical methods,
Monte Carlo simulation, hand-built automatic differentiation, and
calibration against live market data — validated against a real
exchange, not just textbook values.

[![CI](https://github.com/roberto0607/deriv-engine/actions/workflows/ci.yml/badge.svg)](https://github.com/roberto0607/deriv-engine/actions)

## Headline result

Pulled the full live BTC options chain from Deribit (904 real
contracts) and independently solved implied volatility for every one,
using nothing but this engine's own Black-Scholes implementation and a
hand-built Newton-Raphson solver — no reference to Deribit's own
reported values during solving.

**Result: 95.1% convergence, averaging 1.9 percentage points from
Deribit's own reported implied volatility across the entire chain.**

The resulting volatility smile, plotted directly from real market
prices:

*(insert the smile chart image here — screenshot or re-render the
chart from Phase 8 and drop it in as `docs/vol_smile.png`)*

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

Every method above independently arrives at the same price for the
same European option — four structurally different approaches
(formula, tree, simulation, PDE) agreeing is the core correctness
evidence for this engine.

## Notable engineering details

- **Every bug found is documented, not hidden.** Three real bugs were
  found and fixed during development, each with a before/after test
  proving the fix: a division-by-zero at T=0, a double-counted
  premium in the hedging PnL calculation, and a random-seed mismatch
  that corrupted a finite-difference comparison. See `docs/phase-log.md`.
- **AAD was benchmarked honestly, not just claimed to be fast.**
  At the tested scale (4 Greeks, 100k paths), hand-built reverse-mode
  AAD was measured to be *slower* than naive bump-and-revalue — the
  reasoning for why, and where the actual crossover point is, is
  documented in `docs/numerics.md` rather than glossed over.
- **Deliberate scope decisions, not gaps.** Discrete dividends,
  trinomial trees, and second-order AAD (gamma) were all considered
  and explicitly deferred, with the reasoning written down — see
  `docs/numerics.md`.
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
- [`docs/phase-log.md`](docs/phase-log.md) — a build log across all 8 phases: what was built, how it was validated, and every bug found along the way

## Tech stack

C++20 · CMake · Catch2 · nlohmann/json · Deribit public API

## Scope

Built around BTC options specifically (continuous dividend yield = 0,
correct for an asset with no dividends) rather than equities. See
`docs/numerics.md` for the reasoning behind this and other scope
decisions.