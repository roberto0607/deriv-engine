# Design

## Scope

Instruments: European and American vanilla calls/puts only. No exotics
(Asian, Barrier) in the core scope — may revisit as a stretch goal once
everything below is solid and defensible.

## Architecture

Three-layer separation, deliberately decoupled so no layer knows about
the internals of another:

Instrument (what) ---> PricingModel (how) ---> Pricer<Model> (glue)
|
MarketData (inputs)


### `Instrument`
Base class describing the contract itself: strike, expiry, option type
(call/put), exercise style (European/American). Doesn't know how it's
priced — an instrument is just a description of the payoff.

- `EuropeanOption : Instrument`
- `AmericanOption : Instrument`

### `MarketData`
Plain struct: spot price, risk-free rate curve (flat rate to start),
dividend yield, volatility (flat to start, later a surface). Passed
into pricing calls rather than owned by the instrument — the same
instrument can be priced under different market conditions without
mutation.

### `PricingModel`
Strategy interface. Each concrete model (Black-Scholes closed-form,
CRR binomial tree, Monte Carlo, Longstaff-Schwartz, Crank-Nicolson PDE)
implements this interface. Adding a new model means adding a new class,
not touching existing ones.

### `Pricer<Model>`
Templated wrapper that binds an `Instrument` + `MarketData` to a
specific `PricingModel` at compile time. Avoids virtual dispatch
overhead in hot paths (e.g. Monte Carlo inner loops) while keeping the
model-selection code simple.

## Why this shape

- **Decoupling Instrument from PricingModel** means every new pricing
  method (tree, MC, PDE) is cross-validated against the same instrument
  definitions — no risk of subtly different payoff logic per method.
- **Strategy pattern over inheritance-heavy design** keeps each pricing
  method's numerics isolated and independently testable.
- **Template `Pricer<Model>`** over a virtual base avoids the vtable
  indirection cost in performance-critical paths (multithreaded/SIMD/
  CUDA Monte Carlo), where every avoided branch matters at scale.

## Tech stack

- CMake + FetchContent (Catch2, Eigen when needed) — no manual dependency install
- Eigen for linear algebra (PDE tridiagonal solves, calibration) — not hand-rolled
- Catch2 for tests
- Google Benchmark for performance harness (added when benchmarking starts)
- Hand-rolled reverse-mode AAD tape (not a library) — deliberate choice,
  see numerics.md for why
- CUDA for GPU Monte Carlo (stretch goal, only after CPU path is solid)
- pybind11 for Python bindings (added in calibration phase)

## Open design questions (revisit as they come up)

- Exact shape of the AAD tape (operator overloading on a custom `Number`
  type) — deferred until Phase in which AAD is implemented (weeks 9-11)
- Whether MarketData holds a flat rate or a curve from day one — starting
  flat, upgrading to a curve when calibration phase requires it