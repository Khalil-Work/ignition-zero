// Public interface of ignition_zero (see ignition_zero.cpp for the model).
#ifndef IGNITION_ZERO_HPP
#define IGNITION_ZERO_HPP

namespace ignition_zero {

struct Params {
    double A;        // pre-exponential factor [1/s]
    double Ea;       // activation energy [J/mol]
    double R;        // universal gas constant, same molar basis [J/(mol K)]
    double q;        // heat release per unit fuel mass [J/kg]
    double cp;       // specific heat (constant) [J/(kg K)]
    double T0;       // initial temperature [K]
    double Y0;       // initial fuel mass fraction [-]
    double dt;       // time step [s]
    int n_steps;     // number of steps to integrate
};

// Fixed-step RK4. Every intermediate stage and the final state must be
// finite, physical (Y >= 0, T > 0), and monotone (Y non-increasing; T
// non-decreasing in the adiabatic mode) - unstable steps are rejected,
// never clipped. Error codes: 0 OK, 1 invalid input, 2 step left the
// physical domain, 3 non-finite state, 4 output buffer too small (a
// non-null history array with len < n_steps is rejected rather than
// truncated; pass null pointers to ignore history).
int integrate_adiabatic(const Params& p, double* t_out, double* Y_out,
                        double* T_out, int len);

// Isothermal (thermostatted) benchmark: exact Y(t) = Y0 exp(-k t).
// No adiabatic energy conservation is claimed in this mode. The same
// stage checks apply; an oversized dt*k that would make RK4 grow Y is
// rejected with error 2.
int integrate_isothermal(const Params& p, double T_fixed, double* t_out,
                         double* Y_out, int len);

// RK4 with positivity/monotonicity step-halving retry (NOT adaptive error
// control: there is no truncation-error estimate). Contract: exactly
// n_steps ACCEPTED steps are taken; t_out records the ACTUAL accumulated
// times after any halving (use those, not (i+1)*dt, for every comparison
// and plot). Same buffer contract as above.
int integrate_quench_safe(const Params& p, double* t_out, double* Y_out,
                          double* T_out, int len);

double arrhenius(double T, const Params& p);
double energy_invariant(const Params& p, double Y, double T);

}  // namespace ignition_zero

#endif  // IGNITION_ZERO_HPP
