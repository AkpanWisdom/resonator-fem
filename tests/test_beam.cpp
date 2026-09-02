// Verification tests. Returns nonzero if any check fails.

#include <cstdio>
#include <cmath>
#include "beam.h"

static int check(const char* name, double got, double want, double tol) {
    double err = std::abs(got - want) / std::abs(want);
    printf("%-40s got %.8e  want %.8e  rel %.2e  %s\n",
           name, got, want, err, err < tol ? "PASS" : "FAIL");
    return err < tol ? 0 : 1;
}

int main() {
    int fail = 0;
    double L = 1.0, EI = 1.0, rhoA = 1.0;

    // Static deflection under a midspan point load: exact for cubic elements.
    {
        BeamModel beam(20, L, EI, rhoA);
        Eigen::MatrixXd K = beam.AssembleStiffness();
        std::vector<int> free = beam.GetFreeDofsClampedClamped();
        Eigen::MatrixXd Kr = beam.Reduce(K, free);
        int middof = 2 * 10;
        Eigen::VectorXd F = Eigen::VectorXd::Zero(free.size());
        for (int i = 0; i < (int)free.size(); i++)
            if (free[i] == middof) F(i) = 1.0;
        Eigen::VectorXd u = Kr.ldlt().solve(F);
        double d = 0.0;
        for (int i = 0; i < (int)free.size(); i++)
            if (free[i] == middof) d = u(i);
        fail += check("static midspan deflection", d, 1.0 / 192.0, 1e-10);
    }

    // Natural frequencies converge to the analytical values.
    {
        BeamModel beam(64, L, EI, rhoA);
        Eigen::VectorXd f = beam.ComputeNaturalFrequencies(0.0);
        for (int n = 0; n < 3; n++) {
            char name[64];
            snprintf(name, sizeof(name), "natural frequency mode %d", n + 1);
            fail += check(name, f(n), beam.ComputeExactFrequency(n + 1), 1e-5);
        }
    }

    // Total element mass equals rho*A*L.
    {
        BeamModel beam(8, L, EI, rhoA);
        Eigen::MatrixXd M = beam.AssembleMass();
        double total = 0.0;
        for (int i = 0; i < M.rows(); i += 2)
            for (int j = 0; j < M.cols(); j += 2)
                total += M(i, j);
        fail += check("total mass", total, rhoA * L, 1e-12);
    }

    // Buckling: f1 collapses as compression approaches N_cr.
    {
        BeamModel beam(32, L, EI, rhoA);
        double Ncr = 4.0 * M_PI * M_PI * EI / (L * L);
        double f0 = beam.ComputeNaturalFrequencies(0.0)(0);
        double f = beam.ComputeNaturalFrequencies(-0.99 * Ncr)(0);
        printf("%-40s f1/f0 = %.4f  %s\n",
               "buckling limit", f / f0, f / f0 < 0.15 ? "PASS" : "FAIL");
        if (f / f0 >= 0.15) fail++;
    }

    printf("\n%s\n", fail == 0 ? "all tests passed" : "TESTS FAILED");
    return fail;
}
