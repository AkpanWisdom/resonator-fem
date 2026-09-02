// Convergence of the first three natural frequencies of a clamped-clamped
// beam against the Euler-Bernoulli analytical solution.

#include <cstdio>
#include <cmath>
#include "beam.h"

int main() {
    double length = 1.0, EI = 1.0, rhoA = 1.0;
    int meshes[5] = {4, 8, 16, 32, 64};

    BeamModel ref(4, length, EI, rhoA);
    double fex[3];
    for (int n = 0; n < 3; n++) fex[n] = ref.ComputeExactFrequency(n + 1);

    printf("Clamped-clamped beam, L=1, EI=1, rhoA=1\n");
    printf("exact: f1=%.6f  f2=%.6f  f3=%.6f\n\n", fex[0], fex[1], fex[2]);
    printf("nelem        f1     err1 %%         f2     err2 %%         f3     err3 %%\n");

    for (int m = 0; m < 5; m++) {
        BeamModel beam(meshes[m], length, EI, rhoA);
        Eigen::VectorXd f = beam.ComputeNaturalFrequencies(0.0);
        printf("%5d", meshes[m]);
        for (int n = 0; n < 3; n++) {
            printf("  %9.6f  %9.4f", f(n), std::abs(f(n) - fex[n]) / fex[n] * 100.0);
        }
        printf("\n");
    }
    return 0;
}
