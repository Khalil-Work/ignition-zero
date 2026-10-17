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

// Fixed-step RK4; returns 0 on success (error codes in the .cpp).
int integrate_adiabatic(const Params& p, double* t_out, double* Y_out,
                        double* T_out, int len);

// Isothermal (thermostatted) benchmark: exact Y(t) = Y0 exp(-k t).
// No adiabatic energy conservation is claimed in this mode.
int integrate_isothermal(const Params& p, double T_fixed, double* t_out,
                         double* Y_out, int len);

// RK4 with step-halving protection against leaving the physical domain.
int integrate_quench_safe(const Params& p, double* t_out, double* Y_out,
                          double* T_out, int len);

double arrhenius(double T, const Params& p);
double energy_invariant(const Params& p, double Y, double T);

}  // namespace ignition_zero

#endif  // IGNITION_ZERO_HPP
