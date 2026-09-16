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

## Automatic Adjoint Differentiation (AAD)

**Why hand-rolled instead of a library (XAD/CoDiPack):**
Understanding the tape mechanics well enough to build it is the
resume-defining piece of this project — using a library would mean not
being able to defend the mechanism under questioning.

**Tape design:**
`Tape` is a flat `vector<Node>`, where each `Node` stores up to two
parent indices and the local derivative ("partial") with respect to
each parent. `ADouble` wraps a `double` value together with its index
on the tape. Constructing an `ADouble` from a raw `double` creates a
new leaf node (used only for genuine inputs — S, sigma, r, T — that a
derivative is wanted for); every arithmetic operator (+, -, *, /,
ad_exp, ad_log, ad_sqrt) records a new node with the correct local
derivative, without the caller needing to think about calculus at all.

Reverse-mode was chosen over forward-mode specifically because it
computes derivatives with respect to *all* inputs (all 5 potential
Greeks) in a single backward pass, regardless of how many inputs
there are — forward mode would need one full forward pass per input,
which is a different but structurally similar cost problem to
bump-and-revalue.

**Backward pass mechanics:**
`Tape::backward(root)` seeds the root node's adjoint at 1.0, then
walks every node from last-recorded to first. Because a node's parents
are always recorded before the node itself, walking in strict reverse
order guarantees each node's adjoint is fully accumulated (from every
child that used it) before it's propagated further back — this is why
a single linear backward pass correctly implements the chain rule
across the whole graph, including nodes with multiple children (like
S being used in both the drift and the S_T computation).

**Discounting folded into the tape:**
The discount factor e^(−rT) is computed as `ad_exp(-(r*T))` and
multiplied into the payoff *on the tape*, rather than applied as a
separate plain-double multiply after backward() runs. This matters
specifically for rho: r affects the price two distinct ways (through
the drift inside the payoff, and through discounting itself), and
folding both onto the same tape means the backward pass automatically
combines both effects via the chain rule — no manual reasoning about
how the two contributions add together was required.

**Handling the payoff's kink:**
max(S_T − K, 0) isn't differentiable at the single point S_T = K. The
implementation checks which branch applies using the plain `.value`
(not the tape), and returns whichever branch's smooth derivative is
correct for that path. This is valid because that exact knife-edge
point has zero probability of being hit in a continuous simulation —
a standard, correct way to handle this in practice.

**Validated against:**
Delta, vega, theta, and rho from `monte_carlo_greeks` (single
simulation, AAD) all matched the closed-form analytical Greeks from
Phase 1 within tight tolerance, on the standard reference case
(S=100, K=100, T=1, r=5%, vol=20%).

**Pathwise differentiation — the underlying technique:**
Differentiating the pricing formula while holding the random draw Z
fixed for each path is a named technique (pathwise differentiation /
infinitesimal perturbation analysis). AAD is the mechanism used to
compute those pathwise derivatives efficiently — not a separate
technique from pathwise differentiation.

**AAD vs. bump-and-revalue — measured result (an honest one):**
Benchmarked computing 4 Greeks (delta, vega, theta, rho) for 100,000
paths: bump-and-revalue (1 base run + 4 bumped runs, sharing a fixed
random seed across all 5 runs so the comparison isolates each input's
true effect) took 94.3ms; AAD (1 run) took 150.5ms. **At this scale,
AAD was slower, not faster.**

This is a real, informative result, not a failed implementation: the
tape allocates and walks a full computation graph on every single
path, and that per-path overhead (heap allocation, backward-pass
traversal) is real cost a plain bump-and-revalue loop never pays.
AAD's actual structural advantage is that its per-path cost stays
FIXED regardless of how many Greeks are requested, while
bump-and-revalue's cost grows linearly with the Greek count (one full
extra simulation per Greek). At 4 Greeks, bump-and-revalue's added
simulation cost apparently hasn't yet crossed AAD's fixed tape
overhead. The crossover point would likely appear at a higher Greek
count, for a more complex instrument needing many risk factors, or
with a lower per-path tape overhead — a possible future optimization,
not attempted here given the phase's remaining scope.

**A real bug found in the process:**
The first version of this benchmark showed bump-and-revalue's delta
disagreeing with AAD's by 3.2 — far too large to be sampling noise.
Root cause: `monte_carlo_price` created a fresh, unseeded RNG on every
call, so the "base" and "bumped" runs in the bump-and-revalue
comparison used completely unrelated random paths, not the same
underlying randomness with one input nudged — which defeats the whole
point of a finite-difference comparison. Fixed by adding an optional
`seed` parameter (default 0 preserves existing random behavior
everywhere else; a nonzero seed forces reproducible draws). After the
fix, bump-delta (0.6363) and AAD-delta (0.6376) agree closely, as
expected.

**Gamma — deliberately deferred:**
Gamma is a second-order derivative (∂²price/∂S², the derivative of
delta), which the current first-order reverse-mode tape cannot
compute directly — it would require a genuinely more complex "tape of
tapes" structure to differentiate the backward pass itself. Given
remaining project scope (Longstaff-Schwartz, PDE solver, hedging
simulation, calibration), this was judged lower-value than building
that additional complexity: it would add depth to a technique already
proven working (AAD) rather than covering new ground. Gamma remains
available via the exact analytical Black-Scholes formula (Phase 1)
for any downstream use, such as the delta-hedging simulation.