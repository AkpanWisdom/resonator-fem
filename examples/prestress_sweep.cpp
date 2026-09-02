// Fundamental frequency of a clamped-clamped beam under axial prestress.
// Verifies the buckling limit and the high-tension string limit.

#include <cstdio>
#include <cmath>
#include "beam.h"

int main() {
    int nelem = 32;
    double length = 1.0, EI = 1.0, rhoA = 1.0;

    BeamModel beam(nelem, length, EI, rhoA);
    double Ncr = 4.0 * M_PI * M_PI * EI / (length * length);
    double f0 = beam.ComputeNaturalFrequencies(0.0)(0);

    printf("N_cr = %.4f     f1(N=0) = %.6f\n\n", Ncr, f0);
    printf("   N/N_cr           N         f1     f1/f0\n");

    double ratios[9] = {-0.99, -0.75, -0.50, -0.25, 0.0, 1.0, 10.0, 100.0, 1000.0};
    for (int i = 0; i < 9; i++) {
        double N = ratios[i] * Ncr;
        double f1 = beam.ComputeNaturalFrequencies(N)(0);
        printf("%9.2f  %10.2f  %9.6f  %8.4f\n", ratios[i], N, f1, f1 / f0);
    }

    printf("\nstring limit at high tension:\n");
    printf("      N        f1_fem   f1_string     ratio\n");
    double Nbig[3] = {1.0e4, 1.0e5, 1.0e6};
    for (int i = 0; i < 3; i++) {
        double f1 = beam.ComputeNaturalFrequencies(Nbig[i])(0);
        double fstr = sqrt(Nbig[i] / rhoA) / (2.0 * length);
        printf("%7.0e  %12.4f  %10.4f  %8.4f\n", Nbig[i], f1, fstr, f1 / fstr);
    }
    return 0;
}
