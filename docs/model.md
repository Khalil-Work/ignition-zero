# Model and method

## Problem

A reduced zero-dimensional reacting system with Arrhenius kinetics and
adiabatic energy coupling — an *illustrative* single-step mechanism, not a
detailed fuel mechanism and not a reacting CFD solver:

    dY/dt = −k(T)·Y,      dT/dt = (q/cₚ)·k(T)·Y,      k(T) = A·exp(−Ea/(R·T))

- Y: fuel mass fraction [-]; T: temperature [K]
- A: pre-exponential factor [1/s]; Ea: activation energy [J/mol];
  R: gas constant on the SAME molar basis [J/(mol·K)] (Ea and R must share
  that basis — the unit agreement is explicit in the API)
- q: heat release per unit fuel mass [J/kg]; cₚ: constant specific heat
  [J/(kg·K)]

## Exact properties

- **Adiabatic invariant:** for this constant-cₚ adiabatic model,
  cₚ·T + q·Y is conserved exactly along trajectories; the code exposes it
  as `energy_invariant` and the tests pin it to roundoff. Consequently the
  long-time limit is T → T₀ + q·Y₀/cₚ as Y → 0.
- **Isothermal benchmark:** with an external thermostat holding T fixed
  (a *separate* mode with its own exact solution Y(t) = Y₀·e^{−k t}), no
  energy conservation is claimed — the reservoir supplies/removes heat.

## Numerical method

Classical RK4 with a user step in `integrate_adiabatic` /
`integrate_isothermal`. The model is positive but can be stiff near
runaway; the integrators enforce the **physical domain** at every accepted
step: Y ≥ 0, T > 0, finiteness, AND monotonicity (Y non-increasing; T
non-decreasing in the exothermic adiabatic mode, q ≥ 0). The monotonicity
check matters: an unstable RK4 step can produce formally positive but
physically impossible results (fuel increasing, temperature dropping) —
stage garbage is rejected, never clipped.

`integrate_quench_safe` wraps the same RK4 with step halving (floor
2⁻⁶⁰·dt) for stiff cases; it still enforces all domain/monotonicity checks
and conserves the invariant. Error flags: 0 OK, 1 invalid input,
2 step left the physical domain, 3 non-finite state.

## Scope

Demonstrated range: the parameter sets in tests/examples (moderate
activation energy, ms–s scales). This is not a stiff-ODE solver and not a
chemical kinetics mechanism: explicit RK4 with safeguards, single-step
chemistry, constant properties.
