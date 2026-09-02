// =============================================================================
// Environmental history integrator.
// =============================================================================

#include <cmath>
#include "environment.h"

EnvironmentModel::EnvironmentModel(int nkelv, double E, double alphaT, double alphaH)
    : m_nkelv(nkelv), m_E(E), m_alphaT(alphaT), m_alphaH(alphaH) {
    m_c0 = 0.0;
    m_c1 = 0.0;
    m_cp = 1.0;
    m_S = 1.0;
    m_sigma = 0.0;
    m_A0 = 1.0 / E;
    m_T0 = 293.15;
    m_h0 = 0.5;
    m_Tprev = m_T0;
    m_hprev = m_h0;
}

void EnvironmentModel::SetMicroprestress(double c0, double c1, double cp, double S0) {
    m_c0 = c0;
    m_c1 = c1;
    m_cp = cp;
    m_S = S0;
}

// Retardation times spaced by decades. Unit compliances follow a power law
// in tau, which gives the log-time creep response typical of these materials.
void EnvironmentModel::ComputeKelvinProps(double tau0, double A0) {
    m_A0 = A0;
    m_tau.resize(m_nkelv);
    m_A.resize(m_nkelv);
    m_gamma.assign(m_nkelv, 0.0);

    m_tau[0] = tau0;
    for (int i = 1; i < m_nkelv; i++) {
        m_tau[i] = 10.0 * m_tau[i - 1];
    }
    for (int i = 0; i < m_nkelv; i++) {
        m_A[i] = A0 * pow(m_tau[i] / tau0, 0.25);
    }
}

void EnvironmentModel::Reset(double T0, double h0) {
    m_T0 = T0;
    m_h0 = h0;
    m_Tprev = T0;
    m_hprev = h0;
    m_sigma = 0.0;
    m_gamma.assign(m_nkelv, 0.0);
}

double EnvironmentModel::Advance(double T, double h, double dt) {
    if (dt <= 0.0) {
        return m_sigma;
    }

    // Microprestress evolution. The driving term is the rate of change of
    // T*ln(h), so the state depends on how fast the environment moved.
    // Implicit on the decay term keeps S positive for any step size.
    double drive = std::abs(T * log(h) - m_Tprev * log(m_hprev)) / dt;
    m_S = (m_S + m_c1 * drive * dt) / (1.0 + m_c0 * m_S * dt);
    if (m_S < 1.0e-12) {
        m_S = 1.0e-12;
    }

    // Free environmental strain relative to the reference state.
    double eps_env = m_alphaT * (T - m_T0) + m_alphaH * (h - m_h0);

    // Kelvin chain. Each unit integrates exactly over the step:
    //     gam_new = A_i sigma (1 - beta_i) + gam_old beta_i
    // Substituting into the restraint condition
    //     A0 sigma + sum(gam_new) + eps_env = 0
    // gives sigma implicitly. Viscosity scales as 1/(cp*S), so a higher
    // microprestress shortens the retardation times.
    double eta_factor = 1.0 / (m_cp * m_S);
    double num = eps_env;
    double den = m_A0;
    std::vector<double> beta(m_nkelv);
    for (int i = 0; i < m_nkelv; i++) {
        double tau = m_tau[i] * eta_factor;
        beta[i] = exp(-dt / tau);
        num += m_gamma[i] * beta[i];
        den += m_A[i] * (1.0 - beta[i]);
    }

    m_sigma = -num / den;

    for (int i = 0; i < m_nkelv; i++) {
        m_gamma[i] = m_A[i] * m_sigma * (1.0 - beta[i]) + m_gamma[i] * beta[i];
    }

    m_Tprev = T;
    m_hprev = h;
    return m_sigma;
}
