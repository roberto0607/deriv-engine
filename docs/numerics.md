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