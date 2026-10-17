// Demo driver: runs the adiabatic ignition case and writes CSV history for
// plotting. Output: stdout summary + build/ignition_history.csv and
// build/energy_audit.csv.
#include <cmath>
#include <cstdio>

#include "ignition_zero.hpp"

int main() {
    ignition_zero::Params p{};
    p.A = 1.0e6;         // 1/s
    p.Ea = 8.0e4;        // J/mol
    p.R = 8.314462618;   // J/(mol K)
    p.q = 2.0e6;         // J/kg
    p.cp = 1200.0;       // J/(kg K)
    p.T0 = 1100.0;       // K
    p.Y0 = 0.08;
    p.dt = 1.0e-5;
    p.n_steps = 12000;

    double* t = new double[p.n_steps];
    double* Y = new double[p.n_steps];
    double* T = new double[p.n_steps];

    int rc = ignition_zero::integrate_adiabatic(p, t, Y, T, p.n_steps);
    if (rc != 0) {
        std::fprintf(stderr, "demo: integration failed with code %d\n", rc);
        return 1;
    }

    std::FILE* f = std::fopen("build/ignition_history.csv", "w");
    if (!f) {
        std::fprintf(stderr, "demo: cannot open build/ignition_history.csv\n");
        return 1;
    }
    std::fprintf(f, "time_s,Y_fuel,T_K,k_1_s,chem_energy_J_per_kg,sensible_J_per_kg\n");
    for (int i = 0; i < p.n_steps; ++i) {
        double k = ignition_zero::arrhenius(T[i], p);
        std::fprintf(f, "%.9e,%.12e,%.9f,%.12e,%.12e,%.12e\n", t[i], Y[i], T[i],
                     k, p.q * Y[i], p.cp * (T[i] - p.T0));
    }
    std::fclose(f);

    double Yf = Y[p.n_steps - 1], Tf = T[p.n_steps - 1];
    std::printf("final: Y=%.6e T=%.3f K\n", Yf, Tf);
    std::printf("released q*(Y0-Yf) = %.6e J/kg\n", p.q * (p.Y0 - Yf));
    std::printf("sensible cp*(Tf-T0) = %.6e J/kg\n", p.cp * (Tf - p.T0));
    std::printf("invariant drift: %.3e (relative)\n",
                std::fabs(ignition_zero::energy_invariant(p, Yf, Tf) -
                          ignition_zero::energy_invariant(p, p.Y0, p.T0)) /
                    ignition_zero::energy_invariant(p, p.Y0, p.T0));

    delete[] t;
    delete[] Y;
    delete[] T;
    return 0;
}
