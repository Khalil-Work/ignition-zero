# Verification evidence

Reproduced by building and running `tests/test_ignition.cpp`
(`make test`; exit code is the CI gate) and
`python examples/run_verification.py` for the recorded figures/GIF/CSV.
Verified with GCC 16.2.0 (WinLibs), `-std=c++17 -O2 -Wall -Wextra
-Wpedantic`, Windows 11; CI repeats the build+test on Ubuntu g++.

## What each check establishes

| Check | Property established |
| --- | --- |
| Isothermal exact benchmark | Y(t) matches Y₀·e^{−k t} to < 1e-9 relative at the tested dt. |
| RK4 order | Halving dt on the isothermal case reduces the error by ≈ 2⁴ (order > 3.5 asserted) — measured, not assumed. |
| Adiabatic invariant | cₚ·T + q·Y conserved to < 1e-12 relative over the runaway. |
| Long-time limit | T → T₀ + q·Y₀/cₚ (1233.333 K case) and Y → ~0 at long times. |
| Arrhenius | Matches its own formula exactly; monotone increasing in T. |
| Invalid inputs | Negative Ea, NaN parameters, and **±∞ parameters** rejected with flag 1; infinite `T_fixed` rejected. |
| Oversized-step contract (isothermal) | At dt·k = 10, raw RK4 would return Y ≈ 291·Y_old (positive but GROWING fuel); the monotonicity stage checks reject it with error 2 — no silent unphysical growth can be returned as success. |
| Output-buffer contract | A non-null history array with `len < n_steps` returns error 4 instead of silently truncating; null pointers (final-state-only) remain valid. |
| Stage checking | Every RK4 intermediate stage (not just the final state) must be finite, physical, and monotone. |
| Stiff/quench-safe | A deliberately coarse run (A = 1e10, dt = 1e-3) completes via step halving, conserves the invariant, and — since the monotonicity guard was added — cannot accept the unphysical fuel-increasing/temperature-dropping garbage an unstable step produces. |
| Energy audit | Released chemical energy q(Y₀−Y_f) equals sensible rise cₚ(T_f−T₀) to roundoff. |

## Recorded run (`docs/results/`)

- `ignition_runaway.png/.svg` + `ignition_history.csv` (raw C++ demo
  output): fuel collapse, temperature rise to the exact adiabatic limit,
  Arrhenius-rate feedback.
- `energy_audit.png` + `energy_audit.csv`: released = sensible balance to
  6e-7 J/kg on a 1.6e5 J/kg scale (~4e-12 relative).
- `ignition_evolution.gif`: 41 frames of the computed runaway.
- `run_summary.json`: metrics + compiler/runtime context.

## Limits of this evidence

- Illustrative single-step chemistry with constant properties; not a
  validated mechanism, not a reacting-flow solver, no emission/species
  claims.
- Explicit RK4 with safeguards is demonstrated only on the tested
  parameter range; genuinely stiff regimes need dedicated stiff solvers.
- The isothermal mode verifies the exponential benchmark only; no
  energy-conservation claim attaches to it.
- No experimental comparison anywhere.
