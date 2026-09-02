// =============================================================================
// Beam element matrices for prestressed modal analysis
// =============================================================================

#include "beam.h"

BeamModel::BeamModel(int nelem, double length, double EI, double rhoA)
    : m_nelem(nelem), m_length(length), m_EI(EI), m_rhoA(rhoA) {
}

Eigen::Matrix4d BeamModel::ComputeElementStiffness(double Le) {
    double L2 = Le * Le;
    double L3 = Le * L2;
    Eigen::Matrix4d K;
    K <<  12.0,    6.0*Le,  -12.0,    6.0*Le,
           6.0*Le, 4.0*L2,   -6.0*Le, 2.0*L2,
         -12.0,   -6.0*Le,   12.0,   -6.0*Le,
           6.0*Le, 2.0*L2,   -6.0*Le, 4.0*L2;
    return (m_EI / L3) * K;
}

Eigen::Matrix4d BeamModel::ComputeElementMass(double Le) {
    double L2 = Le * Le;
    Eigen::Matrix4d M;
    M << 156.0,     22.0*Le,   54.0,    -13.0*Le,
          22.0*Le,   4.0*L2,   13.0*Le,  -3.0*L2,
          54.0,     13.0*Le,  156.0,    -22.0*Le,
         -13.0*Le,  -3.0*L2,  -22.0*Le,   4.0*L2;
    return (m_rhoA * Le / 420.0) * M;
}

// N > 0 in tension. Total stiffness K + Kg: tension stiffens, compression softens.
Eigen::Matrix4d BeamModel::ComputeElementGeometric(double N, double Le) {
    double L2 = Le * Le;
    Eigen::Matrix4d Kg;
    Kg <<  36.0,    3.0*Le,  -36.0,    3.0*Le,
            3.0*Le, 4.0*L2,   -3.0*Le, -1.0*L2,
          -36.0,   -3.0*Le,   36.0,   -3.0*Le,
            3.0*Le, -1.0*L2,  -3.0*Le,  4.0*L2;
    return (N / (30.0 * Le)) * Kg;
}

Eigen::MatrixXd BeamModel::AssembleStiffness() {
    int ndof = GetNdof();
    double Le = GetElemLength();
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(ndof, ndof);
    Eigen::Matrix4d Ke = ComputeElementStiffness(Le);
    for (int e = 0; e < m_nelem; e++) {
        int base = 2 * e;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                K(base + i, base + j) += Ke(i, j);
            }
        }
    }
    return K;
}

Eigen::MatrixXd BeamModel::AssembleMass() {
    int ndof = GetNdof();
    double Le = GetElemLength();
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(ndof, ndof);
    Eigen::Matrix4d Me = ComputeElementMass(Le);
    for (int e = 0; e < m_nelem; e++) {
        int base = 2 * e;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                M(base + i, base + j) += Me(i, j);
            }
        }
    }
    return M;
}

Eigen::MatrixXd BeamModel::AssembleGeometric(double N) {
    int ndof = GetNdof();
    double Le = GetElemLength();
    Eigen::MatrixXd Kg = Eigen::MatrixXd::Zero(ndof, ndof);
    Eigen::Matrix4d Kge = ComputeElementGeometric(N, Le);
    for (int e = 0; e < m_nelem; e++) {
        int base = 2 * e;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                Kg(base + i, base + j) += Kge(i, j);
            }
        }
    }
    return Kg;
}

// Clamped-clamped: v and theta fixed at both ends, so DOFs 0,1 and ndof-2,ndof-1.
std::vector<int> BeamModel::GetFreeDofsClampedClamped() {
    std::vector<int> free;
    int ndof = GetNdof();
    for (int d = 2; d < ndof - 2; d++) {
        free.push_back(d);
    }
    return free;
}

Eigen::MatrixXd BeamModel::Reduce(const Eigen::MatrixXd& A, const std::vector<int>& dofs) {
    int k = dofs.size();
    Eigen::MatrixXd R(k, k);
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            R(i, j) = A(dofs[i], dofs[j]);
        }
    }
    return R;
}
