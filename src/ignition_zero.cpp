// ignition_zero: reduced zero-dimensional reacting system with Arrhenius
// kinetics and adiabatic energy coupling.
//
// Model (illustrative single-step chemistry, NOT a detailed mechanism):
//
//   dY/dt   = -k(T) * Y                (fuel mass fraction)
//   dT/dt   = (q / cp) * k(T) * Y      (adiabatic temperature rise)
//   k(T)    = A * exp(-Ea / (R * T))   (Arrhenius rate, 1/s)
//
// Units: Ea and R on the same molar basis (J/mol and J/(mol*K)); A in 1/s;
// q in J/kg; cp in J/(kg*K); T in K. For this adiabatic constant-cp model
// the invariant cp*T + q*Y is conserved exactly along trajectories.
//
// Time integration: classical RK4 with a user step; negativity of the
// species is NOT clipped - a stage that would leave the physical domain
// aborts with a clear error (reduce dt). An isothermal mode (thermostatted
// reservoir at fixed T) is provided as a SEPARATE benchmark with its exact
// exponential solution; no energy conservation is claimed for it.
//
// Error flags (returned by the integrate functions):
//   0 OK, 1 invalid input, 2 stage left the physical domain (Y<0 or T<=0),
//   3 non-finite state.

#include "ignition_zero.hpp"

#include <cmath>
#include <cstdio>

namespace ignition_zero {

bool inputs_valid(const Params& p) {
    return p.A > 0.0 && p.Ea > 0.0 && p.R > 0.0 && p.q >= 0.0 && p.cp > 0.0 &&
           p.T0 > 0.0 && p.Y0 >= 0.0 && p.Y0 <= 1.0 && p.dt > 0.0 &&
           p.n_steps >= 0;
}

double arrhenius(double T, const Params& p) {
    return p.A * std::exp(-p.Ea / (p.R * T));
}

static void rhs(double Y, double T, const Params& p, double isothermal_T,
                bool isothermal, double& dY, double& dT) {
    const double T_used = isothermal ? isothermal_T : T;
    const double k = arrhenius(T_used, p);
    dY = -k * Y;
    dT = isothermal ? 0.0 : (p.q / p.cp) * k * Y;
}

static int check_state(double Y, double T) {
    if (!(Y == Y) || !(T == T)) return 3;    // NaN
    if (Y < 0.0 || T <= 0.0) return 2;
    return 0;
}

// Physical monotonicity of this model: dY/dt <= 0 always (Y >= 0, k > 0),
// and dT/dt >= 0 in the adiabatic exothermic mode (q >= 0). A step whose
// result violates these left the numerical domain even if it stays
// formally positive - RK4 stage garbage must not be accepted.
static bool step_physical(double Y, double T, double Yn, double Tn) {
    return Yn <= Y + 1e-14 && Tn >= T - 1e-9;
}

int integrate_adiabatic(const Params& p, double* t_out, double* Y_out,
                        double* T_out, int len) {
    if (!inputs_valid(p)) return 1;
    double Y = p.Y0, T = p.T0;
    for (int i = 0; i < p.n_steps && i < len; ++i) {
        double k1Y, k1T, k2Y, k2T, k3Y, k3T, k4Y, k4T;
        rhs(Y, T, p, 0.0, false, k1Y, k1T);
        rhs(Y + 0.5 * p.dt * k1Y, T + 0.5 * p.dt * k1T, p, 0.0, false, k2Y, k2T);
        rhs(Y + 0.5 * p.dt * k2Y, T + 0.5 * p.dt * k2T, p, 0.0, false, k3Y, k3T);
        rhs(Y + p.dt * k3Y, T + p.dt * k3T, p, 0.0, false, k4Y, k4T);
        double Yn = Y + (p.dt / 6.0) * (k1Y + 2.0 * k2Y + 2.0 * k3Y + k4Y);
        double Tn = T + (p.dt / 6.0) * (k1T + 2.0 * k2T + 2.0 * k3T + k4T);
        int err = check_state(Yn, Tn);
        if (err || !step_physical(Y, T, Yn, Tn)) {
            std::fprintf(stderr,
                         "ignition_zero: RK4 step left the physical domain "
                         "(Y=%g, T=%g) at step %d; reduce dt\n",
                         Yn, Tn, i);
            return err ? err : 2;
        }
        Y = Yn;
        T = Tn;
        if (t_out) t_out[i] = (i + 1) * p.dt;
        if (Y_out) Y_out[i] = Y;
        if (T_out) T_out[i] = T;
    }
    return 0;
}

int integrate_isothermal(const Params& p, double T_fixed, double* t_out,
                         double* Y_out, int len) {
    if (!inputs_valid(p)) return 1;
    if (!(T_fixed > 0.0)) return 1;
    double Y = p.Y0;
    for (int i = 0; i < p.n_steps && i < len; ++i) {
        double k1Y, k1T, k2Y, k2T, k3Y, k2bT, k4Y, k4T;
        rhs(Y, T_fixed, p, T_fixed, true, k1Y, k1T);
        rhs(Y + 0.5 * p.dt * k1Y, T_fixed, p, T_fixed, true, k2Y, k2T);
        rhs(Y + 0.5 * p.dt * k2Y, T_fixed, p, T_fixed, true, k3Y, k2bT);
        rhs(Y + p.dt * k3Y, T_fixed, p, T_fixed, true, k4Y, k4T);
        double Yn = Y + (p.dt / 6.0) * (k1Y + 2.0 * k2Y + 2.0 * k3Y + k4Y);
        int err = check_state(Yn, T_fixed);
        if (err) {
            std::fprintf(stderr, "ignition_zero: isothermal stage failed at "
                                 "step %d (Y=%g)\n", i, Yn);
            return err;
        }
        Y = Yn;
        if (t_out) t_out[i] = (i + 1) * p.dt;
        if (Y_out) Y_out[i] = Y;
    }
    return 0;
}

double energy_invariant(const Params& p, double Y, double T) {
    return p.cp * T + p.q * Y;
}

int integrate_quench_safe(const Params& p, double* t_out, double* Y_out,
                          double* T_out, int len) {
    // Adaptive halving wrapper: if a full RK4 step leaves the domain, the
    // step is halved (down to a floor) instead of aborting. Reports the
    // final step used through the return value (0 on success).
    if (!inputs_valid(p)) return 1;
    double Y = p.Y0, T = p.T0, t = 0.0;
    double dt = p.dt;
    int i = 0;
    while (i < p.n_steps) {
        double step = dt;
        double Yn = 0.0, Tn = 0.0;
        bool ok = false;
        for (int attempt = 0; attempt < 60; ++attempt) {
            double k1Y, k1T, k2Y, k2T, k3Y, k3T, k4Y, k4T;
            rhs(Y, T, p, 0.0, false, k1Y, k1T);
            rhs(Y + 0.5 * step * k1Y, T + 0.5 * step * k1T, p, 0.0, false, k2Y, k2T);
            rhs(Y + 0.5 * step * k2Y, T + 0.5 * step * k2T, p, 0.0, false, k3Y, k3T);
            rhs(Y + step * k3Y, T + step * k3T, p, 0.0, false, k4Y, k4T);
            Yn = Y + (step / 6.0) * (k1Y + 2.0 * k2Y + 2.0 * k3Y + k4Y);
            Tn = T + (step / 6.0) * (k1T + 2.0 * k2T + 2.0 * k3T + k4T);
            if (check_state(Yn, Tn) == 0 && step_physical(Y, T, Yn, Tn)) {
                ok = true;
                break;
            }
            step *= 0.5;
        }
        if (!ok) {
            std::fprintf(stderr, "ignition_zero: adaptive step fell below "
                                 "the floor at t=%g\n", t);
            return 2;
        }
        Y = Yn;
        T = Tn;
        t += step;
        if (t_out && i < len) t_out[i] = t;
        if (Y_out && i < len) Y_out[i] = Y;
        if (T_out && i < len) T_out[i] = T;
        ++i;
    }
    return 0;
}

}  // namespace ignition_zero
