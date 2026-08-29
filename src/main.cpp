#include <iostream>
#include <Eigen/Dense>

int main() {
    Eigen::Matrix2d A;
    A << 2.0, 1.0,
         1.0, 2.0;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> solver(A);
    std::cout << "Eigenvalues: " << solver.eigenvalues().transpose() << "\n";
    return 0;
}
