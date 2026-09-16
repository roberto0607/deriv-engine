# Numerics

This file is the "why" behind every algorithm and numerical choice in
the engine — the answers an interviewer would probe for. Fill in each
section as that phase is completed, not retroactively.

## Black-Scholes closed-form

**Formula:**
Derived from a delta-hedging / replicating-portfolio argument: model
the underlying as geometric Brownian motion (dS = μS dt + σS dW),
build a portfolio long the option and short Δ=∂V/∂S shares of the
underlying, apply Ito's Lemma, and choose Δ so the random terms cancel
exactly. The resulting riskless portfolio must earn the risk-free
rate by no-arbitrage, which gives the Black-Scholes PDE:

  ∂V/∂t + ½σ²S²·∂²V/∂S² + rS·∂V/∂S − rV = 0

Notably, the stock's expected drift (μ) cancels out entirely during
the derivation — option price depends only on volatility and the
risk-free rate, not on any view of expected return. This is the core
of risk-neutral pricing.

Solving the PDE with the boundary condition V(S,T) = max(S−K, 0) for
a call gives the closed-form solution:

  C = S·e^(−qT)·N(d1) − K·e^(−rT)·N(d2)
  P = K·e^(−rT)·N(−d2) − S·e^(−qT)·N(−d1)

  d1 = [ln(S/K) + (r − q + ½σ²)T] / (σ√T)
  d2 = d1 − σ√T

N(d2) is the risk-neutral probability the option finishes
in-the-money. d1 differs from d2 by exactly one σ√T term because d1
weights the *expected asset price* (used in the S·N(d1) term), not
just exercise probability — this correction accounts for the
log-normal distribution's mean being pulled above its median by
volatility.

**Why closed-form is the ground truth:**
Every other pricing method (tree, MC, PDE) is validated against this
for European options, because it's the only method with a provably
exact answer under GBM assumptions.

**Assumptions and where they break:**
- Constant volatility (breaks: real markets have vol smiles/skew)
- Constant risk-free rate
- No dividends (or continuous yield only — see dividends note below)
- Log-normal price distribution (breaks: fat tails, jumps — relevant
  for BTC specifically, which has much fatter tails than the model
  assumes)
- European exercise only (no early exercise value — this is why trees
  are needed for American options)

**Edge cases and why they're numerically tricky:**
- **Near-zero volatility:** formula degenerates, d1/d2 blow up — needs
  explicit floor/guard
- **Deep ITM/OTM:** N(d1), N(d2) saturate near 0 or 1, losing precision
  in the tails. Tested: deep ITM call converges to the forward-price
  limit with delta→1 and gamma→0 (delta is saturated, so it stops
  changing); deep OTM call converges to price→0 and delta→0.
- **Near-expiry (T→0):** d1, d2 → ±infinity since they're divided by
  σ√T, which goes to 0. Found and fixed as a real bug — the original
  implementation produced `nan` at T=0 rather than degrading
  gracefully. Fixed with an early-return guard: at T≤1e-8, return
  intrinsic value directly (max(S−K,0) for a call), which is also the
  exact mathematical limit of the full formula as T→0 — not just a
  defensive patch, the mathematically correct answer.
- **Zero/negative rates:** confirmed via test — nothing in the
  derivation requires r>0. Negative rates (a real historical case in
  Europe/Japan 2014-2022) produce finite, valid prices, and put-call
  parity still holds under them.

**Day-count and dividends:**
Day-count convention is explicit via `year_fraction()` (ACT/365 or
ACT/360), rather than leaving `time_to_expiry` as an ambiguous raw
double — the two conventions produce meaningfully different T values
for the same calendar dates (ACT/360 gives a slightly larger fraction
than ACT/365 for the same actual days, since it divides by a smaller
denominator), which propagates directly into the priced output.

Dividend yield (q) is set to 0 for this project's scope (BTC options).
Bitcoin is not an equity and has no earnings to distribute, so q=0 is
the mathematically correct assumption, not a simplification. Discrete
dividend handling (adjusting spot by the present value of scheduled
cash payments) is explicitly out of scope.

Note for future extension: if the engine were extended to proof-of-
stake assets (e.g. ETH), q would represent staking yield — a real,
continuous yield analogous to a dividend yield, and already supported
by the existing MarketData.dividend_yield field without any code
changes. Funding rates (perpetual futures) and futures-spot basis are
related cost-of-carry concepts but apply to a different instrument
type, not vanilla options pricing directly.

## Binomial trees

**Why trees converge to Black-Scholes:**
As the number of steps → ∞, the discrete up/down random walk
converges to geometric Brownian motion by a CLT-type argument — more,
smaller steps approximate continuous-time random motion increasingly
well. Tested directly: tree price with 10 steps matches Black-Scholes
within ~2%, and with 500 steps matches within ~0.1%, for the same
European option.

**Why d = 1/u specifically:**
This choice makes the tree recombine — an "up then down" path lands
on exactly the same price node as "down then up." Without this
property, the tree would branch into 2^N distinct terminal nodes;
with it, there are only N+1 distinct nodes at step N. This is the key
efficiency trick that makes trees computationally practical.

**Why trees can price American options and closed-form can't:**
At each node, walking backward from expiry, compare the continuation
value (discounted expected value of holding) against the immediate
exercise value, and take whichever is larger:
`value = max(continuation, exercise_value_if_exercised_now)`.
Closed-form Black-Scholes has no mechanism for this node-by-node
decision — it only produces a single value assuming the option is
held to expiry. Tested: American put price ≥ European put price for
identical parameters, confirmed under a nonzero dividend yield case
(early exercise on a non-dividend-paying call is famously never
optimal, so a put or dividend-paying call is required to actually
exercise this code path meaningfully).

**Trinomial trees:**
(not yet implemented — deferred; binomial with sufficient steps is
adequate for this project's accuracy needs, and time is better spent
on Monte Carlo and AAD, which are higher-priority differentiators)

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