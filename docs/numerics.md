# Numerics

This file is the "why" behind every algorithm and numerical choice in
the engine — the answers an interviewer would probe for. Fill in each
section as that phase is completed, not retroactively.

## Black-Scholes closed-form

**Formula:**
(derivation + final formula go here once Phase 1 is done)

**Why closed-form is the ground truth:**
Every other pricing method (tree, MC, PDE) is validated against this
for European options, because it's the only method with a provably
exact answer under GBM assumptions.

**Assumptions and where they break:**
- Constant volatility (breaks: real markets have vol smiles/skew)
- Constant risk-free rate
- No dividends (or continuous yield only — discrete dividends need
  adjustment, see edge cases below)
- Log-normal price distribution (breaks: fat tails, jumps)
- European exercise only (no early exercise value — this is why trees
  are needed for American options)

**Edge cases and why they're numerically tricky:**
- **Near-zero volatility:** formula degenerates, d1/d2 blow up — needs
  explicit floor/guard
- **Deep ITM/OTM:** N(d1), N(d2) saturate near 0 or 1, losing precision
  in the tails
- **Near-expiry (T→0):** d1, d2 → ±infinity, numerical instability
  right at expiry
- **Zero/negative rates:** need to confirm formula still holds
  (it does, but worth stating explicitly why)

**Day-count and dividends:**
(fill in ACT/360 vs ACT/365 choice and discrete dividend handling once
implemented)

## Binomial / Trinomial trees

**Why trees converge to Black-Scholes:**
(fill in: as steps → ∞, the discrete random walk converges to GBM by
CLT-type argument)

**Why trinomial converges faster than binomial:**
(fill in once implemented — extra branch reduces discretization error
per step)

**Why trees can price American options and closed-form can't:**
At each node, compare continuation value vs. immediate exercise value —
closed-form has no mechanism for this decision at each time step.

## Monte Carlo

**Discretization: Euler vs. Milstein**
(fill in once implemented: Milstein adds a correction term for the
diffusion coefficient's derivative, reducing discretization bias
relative to Euler-Maruyama)

**Variance reduction**
- Antithetic variates: (why it works — negatively correlated pairs
  reduce variance of the average)
- Control variates: (using BS price as control — fill in the actual
  variance reduction ratio measured once implemented)

**Why MC needs discretization error consideration at all:**
Unlike closed-form, MC simulates a discrete-time approximation of a
continuous-time process — every step introduces bias that shrinks as
step size → 0, distinct from the Monte Carlo sampling error that
shrinks as paths → ∞.

## Longstaff-Schwartz (American Monte Carlo)

(fill in once implemented: regression-based continuation value
estimation, why regression rather than nested simulation, choice of
basis functions)

## Finite Difference (Crank-Nicolson)

**Why Crank-Nicolson over explicit/implicit:**
(fill in: unconditionally stable like implicit, but second-order
accurate in time like a centered scheme — explicit is conditionally
stable and can blow up for large time steps)

**Stability conditions:**
(fill in once implemented)

## Automatic Adjoint Differentiation (AAD)

**Why hand-rolled instead of a library (XAD/CoDiPack):**
Understanding the tape mechanics well enough to build it is the
resume-defining piece of this project — using a library would mean not
being able to defend the mechanism under questioning.

**Why AAD instead of bump-and-revalue for MC Greeks:**
Bump-and-revalue requires re-running the full simulation once per
Greek (N+1 simulations for N Greeks). Reverse-mode AAD computes all
Greeks in roughly one extra backward pass over the same computation
graph, regardless of how many Greeks are needed.

**Tape design:**
(fill in once implemented: operator overloading on a custom `Number`
type, recording the computation graph, backward pass mechanics)

## Delta-hedging simulation

(fill in once implemented: what the simulated hedging error shows
about the gap between theoretical Greeks and realized P&L, and why
that gap exists — discretization of the hedge rebalancing frequency)