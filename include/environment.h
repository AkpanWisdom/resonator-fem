// =============================================================================
// Environmental history integrator.
//
// Restrained thermal and hygral strain generates axial stress in a clamped
// member. The stress relaxes through a Kelvin chain whose viscosity evolves
// with microprestress S, driven by the rate of environmental change:
//
//     dS/dt + c0 S^2 = c1 |d(T ln h)/dt|
//
// =============================================================================

#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <vector>

class EnvironmentModel {
  public:
    EnvironmentModel(int nkelv, double E, double alphaT, double alphaH);

    /// Microprestress parameters.
    void SetMicroprestress(double c0, double c1, double cp, double S0);

    /// Kelvin chain retardation times, logarithmically spaced from tau0.
    void ComputeKelvinProps(double tau0, double A0);

    /// Advance one step. T in K, h in [0,1], dt in s.
    /// Returns the current axial stress.
    double Advance(double T, double h, double dt);

    double GetStress() const { return m_sigma; }
    double GetMicroprestress() const { return m_S; }
    double GetViscosity() const { return 1.0 / (m_cp * m_S); }

    /// Reset to the reference state (T0, h0) with zero stress.
    void Reset(double T0, double h0);

    /// Override the compliance of one Kelvin unit, used for verification.
    void SetUnitCompliance(int i, double A) { m_A[i] = A; }

  private:
    int m_nkelv;
    double m_E;
    double m_alphaT;
    double m_alphaH;

    double m_c0;
    double m_c1;
    double m_cp;
    double m_S;

    double m_T0;
    double m_h0;
    double m_Tprev;
    double m_hprev;

    double m_sigma;
    double m_A0;                     // instantaneous compliance
    std::vector<double> m_tau;       // retardation times
    std::vector<double> m_A;         // Kelvin unit compliances
    std::vector<double> m_gamma;     // internal strain of each unit
};

#endif
