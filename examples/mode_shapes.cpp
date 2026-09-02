// First three mode shapes of a clamped-clamped beam.

#include <cstdio>
#include "beam.h"

int main() {
    int nelem = 64;
    double length = 1.0;
    BeamModel beam(nelem, length, 1.0, 1.0);

    Eigen::VectorXd f = beam.ComputeNaturalFrequencies(0.0);
    Eigen::MatrixXd phi = beam.ComputeModeShapes(0.0, 3);

    FILE* out = fopen("modes.csv", "w");
    if (!out) {
        printf("could not open modes.csv\n");
        return 1;
    }
    fprintf(out, "x,mode1,mode2,mode3\n");
    for (int i = 0; i <= nelem; i++) {
        double x = i * length / nelem;
        fprintf(out, "%.6f,%.6f,%.6f,%.6f\n", x, phi(i, 0), phi(i, 1), phi(i, 2));
    }
    fclose(out);

    printf("f1 = %.6f  f2 = %.6f  f3 = %.6f Hz\n", f(0), f(1), f(2));
    printf("wrote modes.csv\n");
    return 0;
}
