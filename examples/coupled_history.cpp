// =============================================================================
// Coupled environmental-modal analysis of a clamped-clamped beam.
//
// Environmental history drives restrained axial stress. The stress enters the
// geometric stiffness, so the natural frequency follows the environmental
// history rather than the instantaneous environment.
//
// Three temperature paths, all reaching the same final temperature.
// =============================================================================

#include <cstdio>
#include <cmath>
#include "beam.h"
#include "environment.h"

static double TemperatureAt(int mode, double t_hours, double T0, double dT, double thold) {
    if (mode == 0) {
        double r = t_hours / 1.0;                    // fast ramp, then hold
        return T0 + dT * (r < 1.0 ? r : 1.0);
    }
    if (mode == 1) {
        double r = t_hours / thold;                  // slow ramp
        return T0 + dT * (r < 1.0 ? r : 1.0);
    }
    double r = t_hours / thold;                      // up, down, up
    if (r < 0.33) return T0 + dT * (r / 0.33);
    if (r < 0.66) return T0 + dT * (1.0 - (r - 0.33) / 0.33);
    if (r < 1.0)  return T0 + dT * ((r - 0.66) / 0.34);
    return T0 + dT;
}

int main() {
    int nelem = 32;
    double length = 1.0;
    double EI = 1.0;
    double rhoA = 1.0;
    double area = 1.0e-4;

    BeamModel beam(nelem, length, EI, rhoA);
    double f0 = beam.ComputeNaturalFrequencies(0.0)(0);

    int nkelv = 10;
    double E = 1.0e10;
    double alphaT = 5.0e-6;
    double alphaH = 2.0e-3;
    double T0 = 293.15;
    double h0 = 0.50;
    double dT = 10.0;

    double hour = 3600.0;
    double thold = 200.0;
    int nsteps = 2000;
    double dt_hours = thold / nsteps;

    double sigma_tol = 0.0;   // solve every step for the plotted figure

    printf("f1 at zero stress = %.6f Hz\n\n", f0);

    const char* names[3] = {"A fast ramp", "B slow ramp", "C up-down-up"};

    FILE* out = fopen("history.csv", "w");
    if (!out) {
        printf("could not open history.csv\n");
        return 1;
    }
    fprintf(out, "case,t_hours,T,sigma,f1\n");

    for (int mode = 0; mode < 3; mode++) {
        EnvironmentModel env(nkelv, E, alphaT, alphaH);
        env.SetMicroprestress(1.0e-6, 1.0e-4, 1.0, 1.0);
        env.ComputeKelvinProps(hour, 1.0 / E);
        env.Reset(T0, h0);

        double sigma_last = 0.0;
        double f1 = f0;
        int nsolve = 0;

        for (int k = 1; k <= nsteps; k++) {
            double t_hours = k * dt_hours;
            double T = TemperatureAt(mode, t_hours, T0, dT, thold);
            double sigma = env.Advance(T, h0, dt_hours * hour);

            // Re-solve the eigenproblem only on a significant stress change.
            if (std::abs(sigma - sigma_last) > sigma_tol) {
                f1 = beam.ComputeNaturalFrequencies(sigma * area)(0);
                sigma_last = sigma;
                nsolve++;
            }

            if (k % 20 == 0) {
                fprintf(out, "%s,%.4f,%.4f,%.6e,%.8f\n", names[mode], t_hours, T, sigma, f1);
            }
        }

        double drift = (f1 - f0) / f0 * 1.0e6;
        printf("%-14s  final sigma = %11.4e Pa   f1 = %.6f Hz   drift = %8.1f ppm   eigensolves = %d\n",
               names[mode], env.GetStress(), f1, drift, nsolve);
    }

    fclose(out);
    printf("\nwrote history.csv\n");
    return 0;
}
