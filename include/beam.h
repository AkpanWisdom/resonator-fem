// =============================================================================
// Beam element matrices for prestressed modal analysis
//
// 2-node Euler-Bernoulli beam. Nodal DOFs: {v, theta}.
// Element DOFs: {v1, th1, v2, th2}.
// =============================================================================

#ifndef BEAM_H
#define BEAM_H

#include <vector>
#include <Eigen/Dense>

class BeamModel {
  public:
    BeamModel(int nelem, double length, double EI, double rhoA);

    int GetNdof() const { return 2 * (m_nelem + 1); }
    double GetElemLength() const { return m_length / m_nelem; }

    /// Element matrices for a single element of length Le.
    Eigen::Matrix4d ComputeElementStiffness(double Le);
    Eigen::Matrix4d ComputeElementMass(double Le);
    Eigen::Matrix4d ComputeElementGeometric(double N, double Le);

    /// Global assembly, unconstrained.
    Eigen::MatrixXd AssembleStiffness();
    Eigen::MatrixXd AssembleMass();
    Eigen::MatrixXd AssembleGeometric(double N);

    /// Free DOFs for clamped-clamped end conditions.
    std::vector<int> GetFreeDofsClampedClamped();

    /// Extract the submatrix at the given DOF indices.
    Eigen::MatrixXd Reduce(const Eigen::MatrixXd& A, const std::vector<int>& dofs);

    /// Solve K phi = omega^2 M phi. N > 0 in tension. Frequencies in Hz.
    Eigen::VectorXd ComputeNaturalFrequencies(double N);

    /// Analytical clamped-clamped frequency, mode n = 1,2,3, zero axial force.
    double ComputeExactFrequency(int n);

  private:
    int m_nelem;
    double m_length;
    double m_EI;
    double m_rhoA;
};

#endif
