// Kelvin chain relaxation against the closed-form single-unit solution.

#include <cstdio>
#include <cmath>
#include "environment.h"

int main() {
    double E = 1.0e10, A0 = 1.0 / E;
    double alphaT = 5.0e-6, T0 = 293.15, h0 = 0.50, dT = 10.0;
    double eps0 = alphaT * dT;
    double tau1 = 3600.0, A1 = 2.0 * A0;

    EnvironmentModel env(1, E, alphaT, 0.0);
    env.SetMicroprestress(0.0, 0.0, 1.0, 1.0);
    env.ComputeKelvinProps(tau1, A0);
    env.Reset(T0, h0);
    env.SetUnitCompliance(0, A1);

    double s0 = -eps0 / A0;
    double sinf = -eps0 / (A0 + A1);
    double k = (A0 + A1) / (A0 * tau1);

    double dt = tau1 / 200.0, t = 0.0, maxerr = 0.0;
    for (int step = 1; step <= 1000; step++) {
        env.Advance(T0 + dT, h0, dt);
        t += dt;
        double exact = sinf + (s0 - sinf) * exp(-k * t);
        double err = std::abs(env.GetStress() - exact) / std::abs(exact);
        if (err > maxerr) maxerr = err;
    }

    printf("Kelvin relaxation, max relative error = %.3e  %s\n",
           maxerr, maxerr < 5e-3 ? "PASS" : "FAIL");

    // Constant environment produces no stress.
    EnvironmentModel env2(10, E, alphaT, 0.0);
    env2.ComputeKelvinProps(tau1, A0);
    env2.Reset(T0, h0);
    for (int i = 0; i < 24; i++) env2.Advance(T0, h0, tau1);
    printf("constant environment, sigma = %.3e  %s\n",
           env2.GetStress(), std::abs(env2.GetStress()) < 1e-6 ? "PASS" : "FAIL");

    return (maxerr < 5e-3 && std::abs(env2.GetStress()) < 1e-6) ? 0 : 1;
}
