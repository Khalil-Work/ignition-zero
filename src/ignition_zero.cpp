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
    const bool finite = std::isfinite(p.A) && std::isfinite(p.Ea) &&
                        std::isfinite(p.R) && std::isfinite(p.q) &&
                        std::isfinite(p.cp) && std::isfinite(p.T0) &&
                        std::isfinite(p.Y0) && std::isfinite(p.dt);
    return finite && p.A > 0.0 && p.Ea > 0.0 && p.R > 0.0 && p.q >= 0.0 &&
           p.cp > 0.0 && p.T0 > 0.0 && p.Y0 >= 0.0 && p.Y0 <= 1.0 &&
           p.dt > 0.0 && p.n_steps >= 0;
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

// Domain check for ANY state (final or intermediate): finite, physical.
static int check_state(double Y, double T) {
    if (!std::isfinite(Y) || !std::isfinite(T)) return 3;
    if (Y < 0.0 || T <= 0.0) return 2;
    return 0;
}

// Physical monotonicity of this model: dY/dt <= 0 always (Y >= 0, k > 0),
// and dT/dt >= 0 in the adiabatic exothermic mode (q >= 0). Applied to
// every INTERMEDIATE stage as well as the final result: RK4 stage garbage
// (fuel increasing, temperature dropping) is rejected, never clipped.
static bool stage_physical(double Y, double T, double Y_stage, double T_stage,
                           bool isothermal) {
    if (Y_stage > Y + 1e-14) return false;
    if (!isothermal && T_stage < T - 1e-9) return false;
    return true;
}

// Evaluates one RK4 step for (Y, T); returns 0 when every intermediate
// stage and the final state are finite, physical, and monotone, else the
// error code. On rejection, Yn/Tn are left holding the last PROBED
// (rejected) stage values for diagnostics - the caller may print them but
// must never use them as an accepted state.
static int rk4_step(double Y, double T, const Params& p, double isothermal_T,
                    bool isothermal, double dt, double& Yn, double& Tn) {
    // Initialize outputs to the current state: never read unset memory.
    Yn = Y;
    Tn = T;
    double k1Y, k1T, k2Y, k2T, k3Y, k3T, k4Y, k4T;
    rhs(Y, T, p, isothermal_T, isothermal, k1Y, k1T);
    double Ys = Y + 0.5 * dt * k1Y;
    double Ts = T + 0.5 * dt * k1T;
    int e = check_state(Ys, Ts);
    if (e || !stage_physical(Y, T, Ys, Ts, isothermal)) { Yn = Ys; Tn = Ts; return e ? e : 2; }
    rhs(Ys, Ts, p, isothermal_T, isothermal, k2Y, k2T);
    Ys = Y + 0.5 * dt * k2Y;
    Ts = T + 0.5 * dt * k2T;
    e = check_state(Ys, Ts);
    if (e || !stage_physical(Y, T, Ys, Ts, isothermal)) { Yn = Ys; Tn = Ts; return e ? e : 2; }
    rhs(Ys, Ts, p, isothermal_T, isothermal, k3Y, k3T);
    Ys = Y + dt * k3Y;
    Ts = T + dt * k3T;
    e = check_state(Ys, Ts);
    if (e || !stage_physical(Y, T, Ys, Ts, isothermal)) { Yn = Ys; Tn = Ts; return e ? e : 2; }
    rhs(Ys, Ts, p, isothermal_T, isothermal, k4Y, k4T);
    Yn = Y + (dt / 6.0) * (k1Y + 2.0 * k2Y + 2.0 * k3Y + k4Y);
    Tn = T + (dt / 6.0) * (k1T + 2.0 * k2T + 2.0 * k3T + k4T);
    e = check_state(Yn, Tn);
    if (e || !stage_physical(Y, T, Yn, Tn, isothermal)) return e ? e : 2;
    return 0;
}

// Output-buffer contract: either all pointers are null (no history), or
// len must be >= n_steps; otherwise the requested integration cannot be
// recorded and error 4 (buffer too small) is returned without integrating.
static int buffers_ok(const double* t_out, const double* Y_out,
                      const double* T_out, int len, int n_steps) {
    const bool any = t_out != nullptr || Y_out != nullptr || T_out != nullptr;
    if (!any) return 0;
    return len >= n_steps ? 0 : 4;
}

int integrate_adiabatic(const Params& p, double* t_out, double* Y_out,
                        double* T_out, int len) {
    if (!inputs_valid(p)) return 1;
    int bc = buffers_ok(t_out, Y_out, T_out, len, p.n_steps);
    if (bc) return bc;
    double Y = p.Y0, T = p.T0;
    for (int i = 0; i < p.n_steps; ++i) {
        double Yn, Tn;
        int err = rk4_step(Y, T, p, 0.0, false, p.dt, Yn, Tn);
        if (err) {
            std::fprintf(stderr,
                         "ignition_zero: RK4 step left the physical domain "
                         "(Y=%g, T=%g) at step %d; reduce dt\n",
                         Yn, Tn, i);
            return err;
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
    if (!std::isfinite(T_fixed) || T_fixed <= 0.0) return 1;
    int bc = buffers_ok(t_out, Y_out, nullptr, len, p.n_steps);
    if (bc) return bc;
    double Y = p.Y0;
    for (int i = 0; i < p.n_steps; ++i) {
        double Yn, Tn;
        int err = rk4_step(Y, T_fixed, p, T_fixed, true, p.dt, Yn, Tn);
        if (err) {
            std::fprintf(stderr,
                         "ignition_zero: isothermal step left the physical "
                         "domain (Y=%g) at step %d; dt*k is too large\n",
                         Yn, i);
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

// Positivity/monotonicity retry with step halving (NOT adaptive error
// control: no truncation-error estimate exists here). Contract: exactly
// n_steps ACCEPTED steps are taken while advancing to the recorded actual
// times; the physical horizon reached is therefore >= n_steps * dt / 2^60
// and the per-step recorded times in t_out are the ACTUAL times - use
// them, not (i+1)*dt, for every comparison and plot.
int integrate_quench_safe(const Params& p, double* t_out, double* Y_out,
                          double* T_out, int len) {
    if (!inputs_valid(p)) return 1;
    int bc = buffers_ok(t_out, Y_out, T_out, len, p.n_steps);
    if (bc) return bc;
    double Y = p.Y0, T = p.T0, t = 0.0;
    double dt = p.dt;
    for (int i = 0; i < p.n_steps; ++i) {
        double step = dt;
        double Yn = 0.0, Tn = 0.0;
        bool ok = false;
        for (int attempt = 0; attempt < 60; ++attempt) {
            int err = rk4_step(Y, T, p, 0.0, false, step, Yn, Tn);
            if (err == 0) { ok = true; break; }
            step *= 0.5;
        }
        if (!ok) {
            std::fprintf(stderr, "ignition_zero: retry step fell below the "
                                 "floor at t=%g\n", t);
            return 2;
        }
        Y = Yn;
        T = Tn;
        t += step;
        if (t_out) t_out[i] = t;
        if (Y_out) Y_out[i] = Y;
        if (T_out) T_out[i] = T;
    }
    return 0;
}

}  // namespace ignition_zero
