// =============================================================================
// Static verification: clamped-clamped beam under midspan point load.
// Exact solution: delta = P L^3 / (192 EI)
// =============================================================================

#include <iostream>
#include <iomanip>
#include <cmath>
#include "beam.h"

int main() {
    int nelem = 20;          // even, so a node sits at midspan
    double length = 1.0;
    double EI = 1.0;
    double rhoA = 1.0;
    double P = 1.0;

    BeamModel beam(nelem, length, EI, rhoA);

    Eigen::MatrixXd K = beam.AssembleStiffness();
    std::vector<int> free = beam.GetFreeDofsClampedClamped();
    Eigen::MatrixXd Kr = beam.Reduce(K, free);

    // Load applied at the transverse DOF of the midspan node.
    int middof = 2 * (nelem / 2);
    Eigen::VectorXd F = Eigen::VectorXd::Zero(free.size());
    for (int i = 0; i < (int)free.size(); i++) {
        if (free[i] == middof) {
            F(i) = P;
        }
    }

    Eigen::VectorXd u = Kr.ldlt().solve(F);

    double computed = 0.0;
    for (int i = 0; i < (int)free.size(); i++) {
        if (free[i] == middof) {
            computed = u(i);
        }
    }

    double exact = P * pow(length, 3) / (192.0 * EI);
    double err = std::abs(computed - exact) / exact * 100.0;

    std::cout << std::scientific << std::setprecision(10);
    std::cout << "computed = " << computed << "\n";
    std::cout << "exact    = " << exact << "\n";
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "error    = " << err << " %\n";

    return 0;
}
