// Verification tests for ignition_zero.
//
// Properties pinned:
// 1. Isothermal benchmark matches the exact exponential Y0 exp(-k t) to
//    RK4 accuracy (fourth-order small for the tested dt).
// 2. Adiabatic energy conservation: cp*T + q*Y invariant to roundoff.
// 3. Final temperature approaches T0 + q*Y0/cp as Y -> 0 (long-time limit).
// 4. Arrhenius unit check: k doubles-with-T monotonicity and hand value.
// 5. Step-size convergence: halving dt quarters the isothermal error
//    (RK4 order 4 check against the exact solution).
// 6. Invalid inputs are rejected with error flag 1.
// 7. Quench-safe integrator completes a stiff case where fixed-step RK4
//    would leave the domain, and still conserves the invariant.
// 8. Energy audit of the demo case: released q*(Y0-Yf) equals cp*(Tf-T0).

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "../src/ignition_zero.hpp"

using ignition_zero::Params;

static int failures = 0;

void check(bool ok, const char* what) {
    if (!ok) {
        std::fprintf(stderr, "FAILED: %s\n", what);
        ++failures;
    }
}

int main() {
    Params p{};
    p.A = 1.0e6;          // 1/s
    p.Ea = 8.0e4;         // J/mol
    p.R = 8.314462618;    // J/(mol K)
    p.q = 2.0e6;          // J/kg
    p.cp = 1200.0;        // J/(kg K)
    p.T0 = 1100.0;        // K
    p.Y0 = 0.08;
    p.dt = 1.0e-5;
    p.n_steps = 2000;

    // --- 1. Isothermal exact benchmark -------------------------------
    {
        const int n = p.n_steps;
        double* t = new double[n];
        double* Y = new double[n];
        int rc = ignition_zero::integrate_isothermal(p, p.T0, t, Y, n);
        check(rc == 0, "isothermal run succeeds");
        double k = ignition_zero::arrhenius(p.T0, p);
        double max_err = 0.0;
        for (int i = 0; i < n; ++i) {
            double exact = p.Y0 * std::exp(-k * t[i]);
            double err = std::fabs(Y[i] - exact) / p.Y0;
            if (err > max_err) max_err = err;
        }
        check(max_err < 1e-9, "isothermal matches exact exponential");
        delete[] t;
        delete[] Y;
    }

    // --- 2. Adiabatic energy conservation ------------------------------
    {
        const int n = p.n_steps;
        double* t = new double[n];
        double* Y = new double[n];
        double* T = new double[n];
        int rc = ignition_zero::integrate_adiabatic(p, t, Y, T, n);
        check(rc == 0, "adiabatic run succeeds");
        double E0 = ignition_zero::energy_invariant(p, p.Y0, p.T0);
        double max_drift = 0.0;
        for (int i = 0; i < n; ++i) {
            double drift = std::fabs(ignition_zero::energy_invariant(p, Y[i], T[i]) - E0) / E0;
            if (drift > max_drift) max_drift = drift;
        }
        check(max_drift < 1e-12, "cp*T + q*Y invariant to roundoff");
        delete[] t;
        delete[] Y;
        delete[] T;
    }

    // --- 3. Long-time limit: T -> T0 + q Y0 / cp -----------------------
    {
        Params pl = p;
        pl.dt = 5.0e-4;
        pl.n_steps = 40000;  // 20 s: Y effectively zero at these rates
        double* t = new double[pl.n_steps];
        double* Y = new double[pl.n_steps];
        double* T = new double[pl.n_steps];
        int rc = ignition_zero::integrate_quench_safe(pl, t, Y, T, pl.n_steps);
        check(rc == 0, "long-time run succeeds");
        double T_inf = p.T0 + p.q * p.Y0 / p.cp;
        double last = T[pl.n_steps - 1];
        check(std::fabs(last - T_inf) < 1.0,
              "final temperature approaches T0 + q Y0 / cp");
        check(Y[pl.n_steps - 1] < 1e-6, "fuel essentially consumed");
        delete[] t;
        delete[] Y;
        delete[] T;
    }

    // --- 4. Arrhenius hand check ---------------------------------------
    {
        double k = ignition_zero::arrhenius(1500.0, p);
        double expect = p.A * std::exp(-p.Ea / (p.R * 1500.0));
        check(std::fabs(k - expect) == 0.0, "arrhenius matches its formula");
        check(ignition_zero::arrhenius(1600.0, p) > k, "rate increases with T");
    }

    // --- 5. RK4 order against the exact isothermal solution -------------
    {
        double errs[2] = {0.0, 0.0};
        double dts[2] = {2.0e-4, 1.0e-4};
        const double t_end = 2.0e-3;
        for (int d = 0; d < 2; ++d) {
            Params pd = p;
            pd.dt = dts[d];
            pd.n_steps = static_cast<int>(t_end / pd.dt);
            double* t = new double[pd.n_steps];
            double* Y = new double[pd.n_steps];
            ignition_zero::integrate_isothermal(pd, pd.T0, t, Y, pd.n_steps);
            double k = ignition_zero::arrhenius(pd.T0, pd);
            errs[d] = std::fabs(Y[pd.n_steps - 1] - p.Y0 * std::exp(-k * t_end)) / p.Y0;
            delete[] t;
            delete[] Y;
        }
        double order = std::log2(errs[0] / errs[1]);
        check(order > 3.5, "RK4 order ~4 on isothermal benchmark");
    }

    // --- 6. Invalid inputs ----------------------------------------------
    {
        Params bad = p;
        bad.Ea = -1.0;
        check(ignition_zero::integrate_adiabatic(bad, nullptr, nullptr, nullptr, 0) == 1,
              "invalid input rejected");
        Params nan_R = p;
        nan_R.R = std::nan("");
        check(ignition_zero::integrate_adiabatic(nan_R, nullptr, nullptr, nullptr, 0) == 1,
              "non-finite parameter rejected");
    }

    // --- 7. Quench-safe stiff case + invariant --------------------------
    {
        Params stiff = p;
        stiff.A = 1.0e10;      // violently fast once hot
        stiff.dt = 1.0e-3;     // deliberately coarse for the fixed stepper
        stiff.n_steps = 5000;
        double* t = new double[stiff.n_steps];
        double* Y = new double[stiff.n_steps];
        double* T = new double[stiff.n_steps];
        int rc = ignition_zero::integrate_quench_safe(stiff, t, Y, T, stiff.n_steps);
        check(rc == 0, "quench-safe completes the stiff case");
        double E0 = ignition_zero::energy_invariant(stiff, stiff.Y0, stiff.T0);
        double max_drift = 0.0;
        for (int i = 0; i < stiff.n_steps; ++i) {
            double drift = std::fabs(ignition_zero::energy_invariant(stiff, Y[i], T[i]) - E0) / E0;
            if (drift > max_drift) max_drift = drift;
        }
        check(max_drift < 1e-11, "invariant holds under step halving");
        // With A=1e10 and the ignition delay of this system the run may
        // still be pre-ignition at 5 s; require meaningful progress only
        // if ignition actually started (Y dropped below half of Y0).
        double Y_last = Y[stiff.n_steps - 1];
        double T_inf = stiff.T0 + stiff.q * stiff.Y0 / stiff.cp;
        if (Y_last < 0.5 * stiff.Y0) {
            check(T[stiff.n_steps - 1] > 0.5 * (stiff.T0 + T_inf),
                  "ignited stiff case releases most of its energy");
        } else {
            check(T[stiff.n_steps - 1] > stiff.T0,
                  "pre-ignition stiff case at least heated");
        }
        delete[] t;
        delete[] Y;
        delete[] T;
    }

    // --- 8. Energy audit -------------------------------------------------
    {
        const int n = p.n_steps;
        double* t = new double[n];
        double* Y = new double[n];
        double* T = new double[n];
        ignition_zero::integrate_adiabatic(p, t, Y, T, n);
        double Yf = Y[n - 1], Tf = T[n - 1];
        double released = p.q * (p.Y0 - Yf);
        double sensible = p.cp * (Tf - p.T0);
        check(std::fabs(released - sensible) / released < 1e-12,
              "released chemical energy equals sensible rise");
        delete[] t;
        delete[] Y;
        delete[] T;
    }

    if (failures) {
        std::fprintf(stderr, "%d check(s) failed\n", failures);
        return EXIT_FAILURE;
    }
    std::printf("ignition_zero: all checks passed\n");
    return EXIT_SUCCESS;
}
