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
(fill in once complete)

**Validated against:**
(fill in — known textbook BS values, put-call parity, edge case tests)

**Defend this:**
(fill in — should include: can derive the PDE and closed-form solution
on a whiteboard without notes)

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