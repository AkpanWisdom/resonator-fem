// =============================================================================
// Verification of the Kelvin-chain relaxation under full restraint.
//
// Single Kelvin unit, constant environment after an initial step. The coupled
// system has a closed-form solution:
//
//     sigma(t) = s_inf + (s_0 - s_inf) exp(-k t)
//     s_0   = eps0 / A0
//     s_inf = eps0 / (A0 + A1)
//     k     = (A0 + A1) / (A0 tau1)
// =============================================================================

#include <cstdio>
#include <cmath>
#include "environment.h"

int main() {
    double E = 1.0e10;
    double A0 = 1.0 / E;
    double alphaT = 5.0e-6;
    double T0 = 293.15;
    double h0 = 0.50;
    double dT = 10.0;
    double eps0 = alphaT * dT;

    double tau1 = 3600.0;
    double A1 = 2.0 * A0;

    // Single unit, microprestress switched off so tau is constant.
    EnvironmentModel env(1, E, alphaT, 0.0);
    env.SetMicroprestress(0.0, 0.0, 1.0, 1.0);
    env.ComputeKelvinProps(tau1, A0);
    env.Reset(T0, h0);

    // Override the unit compliance to the analytical value.
    env.SetUnitCompliance(0, A1);

    double s0 = -eps0 / A0;
    double sinf = -eps0 / (A0 + A1);
    double k = (A0 + A1) / (A0 * tau1);

    printf("s_0 = %.4e   s_inf = %.4e   k = %.4e 1/s\n\n", s0, sinf, k);
    printf("  t/tau      sigma_fem     sigma_exact      err %%\n");

    double dt = tau1 / 200.0;
    double t = 0.0;
    for (int step = 1; step <= 1000; step++) {
        env.Advance(T0 + dT, h0, dt);
        t += dt;
        double r = t / tau1;
        if (step % 100 == 0) {
            double exact = sinf + (s0 - sinf) * exp(-k * t);
            double err = std::abs(env.GetStress() - exact) / std::abs(exact) * 100.0;
            printf("%7.2f   %12.5e   %12.5e   %9.5f\n", r, env.GetStress(), exact, err);
        }
    }

    return 0;
}
