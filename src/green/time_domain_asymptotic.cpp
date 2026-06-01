
/*
 * Copyright (c) 2025 Sergio Fernández Ruano / IHCantabria
 *
 * This file is part of SeaMotions Software.
 *
 * SeaMotions is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SeaMotions is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

// External libraries
#include <cmath>

// Local modules
#include "../config.hpp"
#include "../math/math_constants.hpp"
#include "time_domain_asymptotic.hpp"


// ============================================================================
//  Internal implementation
// ============================================================================
namespace {

// ----------------------------------------------------------------------------
//  Precomputed coefficient tables (built once at program start)
// ----------------------------------------------------------------------------

/**
 * Holds the two families of coefficients that appear in the stationary-phase
 * expansion of the transient Green function derivatives:
 *
 *   alpha[n] = (2n+2)! / n!
 *   d[n][m]  = (2m+2n-2)! / ((2n-2)! · 4^m · m!)   for n ≥ 1, 0 ≤ m ≤ n
 *              1                                       for m = 0 (any n)
 *
 * All values are stored as double to preserve precision during the
 * initialisation products; they are cast to cusfloat at the call sites.
 *
 * Recurrences used:
 *   alpha[0]   = 2
 *   alpha[n+1] = alpha[n] · (2n+3)(2n+4) / (n+1)
 *
 *   d[n][0]    = 1
 *   d[n][m]    = d[n][m-1] · (2n+2m-3)(2n+2m-2) / (4m)
 */
struct TDTables {
    double alpha[TD_ASYMP_NMAX];
    double d[TD_ASYMP_NMAX][TD_ASYMP_NMAX];   // d[n][m], valid for m ≤ n

    TDTables() noexcept {
        // --- alpha[n] ---
        alpha[0] = 2.0;
        for (int n = 0; n < TD_ASYMP_NMAX - 1; ++n)
            alpha[n+1] = alpha[n]
                       * static_cast<double>((2*n + 3) * (2*n + 4))
                       / static_cast<double>(n + 1);

        // --- d[n][m] ---
        for (int n = 0; n < TD_ASYMP_NMAX; ++n) {
            d[n][0] = 1.0;
            for (int m = 1; m <= n; ++m)
                d[n][m] = d[n][m-1]
                         * static_cast<double>((2*n + 2*m - 3) * (2*n + 2*m - 2))
                         / static_cast<double>(4 * m);
        }
    }
};

static const TDTables TABLES;

// Convenience accessors (cast to cusfloat at the point of use)
inline cusfloat A(int n)      noexcept { return static_cast<cusfloat>(TABLES.alpha[n]);    }
inline cusfloat D(int n, int m) noexcept { return static_cast<cusfloat>(TABLES.d[n][m]);  }


// ----------------------------------------------------------------------------
//  Shared helper: f2A oscillatory envelope
// ----------------------------------------------------------------------------

/**
 * f2A(β, μ) = (−4i / √(2η)) · exp(−β²μ/4) · exp(i·(β²η/4 − arccos(μ)/2 + π/4))
 *
 * Common to all five derivative functions; depends on (β, μ, η = √(1−μ²)).
 */
inline cuscomplex f2a_envelope(cusfloat beta, cusfloat mu, cusfloat eta) noexcept {
    const cusfloat amp   = cusfloat(-4.0) / std::sqrt(cusfloat(2.0) * eta);
    const cusfloat decay = std::exp(cusfloat(-0.25) * beta * beta * mu);
    const cusfloat phase = cusfloat(0.25) * beta * beta * eta
                         - cusfloat(0.5)  * std::acos(mu)
                         + cusfloat(0.25) * static_cast<cusfloat>(PI);
    // (0 + i·amp) · decay · exp(i·phase)
    return cuscomplex(cusfloat(0), amp)
           * (decay * std::exp(cuscomplex(cusfloat(0), phase)));
}


// ----------------------------------------------------------------------------
//  Shared helper: Legendre polynomials + derivatives via 3-term recurrence
// ----------------------------------------------------------------------------

/**
 * Fills P[0..nmax-1] and dP[0..nmax-1] where
 *   P[n]  = P_n(μ)          (Legendre polynomial)
 *   dP[n] = P_n'(μ) = n · (P[n-1] − μ·P[n]) / η²
 *
 * P and dP must each point to an array of at least nmax elements.
 * dP[0] = 0 always; P[1] = μ, dP[1] = 1.
 */
static void compute_legendre(int nmax, cusfloat mu, cusfloat eta2,
                              cusfloat* P, cusfloat* dP) noexcept {
    if (nmax < 1) return;
    P[0]  = cusfloat(1);
    dP[0] = cusfloat(0);
    if (nmax >= 2) {
        P[1]  = mu;
        // dP[1] = 1·(P[0] − μ·P[1]) / η² = (1−μ²)/η² = 1
        dP[1] = (cusfloat(1) - mu * mu) / eta2;
    }
    for (int n = 2; n < nmax; ++n) {
        P[n]  = (cusfloat(2*n - 1) * mu * P[n-1]
                 - cusfloat(n - 1) * P[n-2]) / cusfloat(n);
        dP[n] = cusfloat(n) * (P[n-1] - mu * P[n]) / eta2;
    }
}

} // anonymous namespace


// ============================================================================
//  Public functions
// ============================================================================
//
//  Inner-loop optimisation strategy
//  ---------------------------------
//  All five functions share the stationary-phase inner double sum
//
//    f2B(β,μ) = Σ_{n=0}^{N-1} (i/η)^n Σ_{m=0}^{n} d[n][m]
//               · (β·z1/2)^{1−2m−2n} · z2^m
//
//  with  z1 = iμ+η,  z2 = μ−iη,  hbz1 = β·z1/2.
//
//  Rather than calling std::pow inside the loops, all β- and μ-dependent
//  powers are advanced by a single complex multiplication per step:
//
//    (j/η)^n  : r0_pow     ×= j/η
//    hbz1^p   : hbz1_p     ×= 1/hbz1²  (p decreases by 2 for each m)
//    z2^m     : z2_pow     ×= z2
//
//  The initial value hbz1_n_pow = hbz1^{1−2n} at the start of each n
//  iteration is maintained by multiplying hbz1_n_pow by 1/hbz1² after
//  each outer step.  This eliminates every std::pow call from the loops.
//
//  All five derivative functions accumulate different subsets of the inner
//  sums S0, S1, S2, dS0, dS1 in a single pass over the same nested loop:
//
//    S0       = Σ r0·r1                          (f2B itself)
//    S1       = Σ r0·(Σ_m p·base)               (β · ∂f2B/∂β)
//    S2       = Σ r0·(Σ_m p(p−1)·base)          (β² · ∂²f2B/∂β²)
//    dS0      = Σ r0·(n·μ/η²·r1 + Σ_m base·c_mn)   (∂f2B/∂μ)
//    dS1      = Σ r0·(n·μ/η²·s1 + Σ_m p·base·c_mn) (β · ∂²f2B/(∂β∂μ))
//
//  where  p = 1−2m−2n  and  c_mn = −p·z2/(η·z1) + m·z1/(η·z2).
// ============================================================================


// ----------------------------------------------------------------------------
//  dGdt
// ----------------------------------------------------------------------------

cusfloat dGdt_asymptotic(cusfloat beta, cusfloat mu, TDAsympConfig cfg) {

    const int      nmax = cfg.nmax;
    const cusfloat eta  = std::sqrt(cusfloat(1) - mu * mu);
    const cuscomplex J  { cusfloat(0), cusfloat(1) };
    const cuscomplex z1 = J * mu + eta;
    const cuscomplex z2 = mu - J * eta;

    const cuscomplex f2a         = f2a_envelope(beta, mu, eta);
    const cuscomplex hbz1        = cusfloat(0.5) * beta * z1;
    const cuscomplex inv_hbz1_sq = cusfloat(1) / (hbz1 * hbz1);
    const cuscomplex j_over_eta  = J / eta;

    // Accumulate S0 = f2B
    cuscomplex r0_pow     { cusfloat(1), cusfloat(0) };
    cuscomplex hbz1_n_pow = hbz1;          // (hbz1)^{1−2n} at n=0 : hbz1^1
    cuscomplex S0 { cusfloat(0), cusfloat(0) };

    for (int n = 0; n < nmax; ++n) {
        cuscomplex hbz1_p = hbz1_n_pow;
        cuscomplex z2_pow { cusfloat(1), cusfloat(0) };
        cuscomplex r1     { cusfloat(0), cusfloat(0) };

        for (int m = 0; m <= n; ++m) {
            r1     += D(n, m) * hbz1_p * z2_pow;
            hbz1_p *= inv_hbz1_sq;
            z2_pow *= z2;
        }

        S0         += r0_pow * r1;
        r0_pow     *= j_over_eta;
        hbz1_n_pow *= inv_hbz1_sq;
    }

    cusfloat result = (f2a * S0).real();

    // f0 = −4 Σ_{n=0}^{N−1} α[n] · β^{−2n−3} · P_n(μ)
    if (cfg.use_f0) {
        cusfloat P[TD_ASYMP_NMAX]{};
        cusfloat dP[TD_ASYMP_NMAX]{};
        compute_legendre(nmax, mu, eta * eta, P, dP);

        const cusfloat inv_beta2 = cusfloat(1) / (beta * beta);
        cusfloat beta_pow        = cusfloat(1) / (beta * beta * beta);  // β^{−3}
        cusfloat f0              = cusfloat(0);

        for (int n = 0; n < nmax; ++n) {
            f0       += A(n) * beta_pow * P[n];
            beta_pow *= inv_beta2;
        }
        result += cusfloat(-4) * f0;
    }

    return result;
}


// ----------------------------------------------------------------------------
//  dGdtt  (∂/∂β of dGdt)
// ----------------------------------------------------------------------------

cusfloat dGdtt_asymptotic(cusfloat beta, cusfloat mu, TDAsympConfig cfg) {

    const int      nmax = cfg.nmax;
    const cusfloat eta  = std::sqrt(cusfloat(1) - mu * mu);
    const cuscomplex J  { cusfloat(0), cusfloat(1) };
    const cuscomplex z1 = J * mu + eta;
    const cuscomplex z2 = mu - J * eta;

    const cuscomplex f2a         = f2a_envelope(beta, mu, eta);
    // df2A/dβ = f2A · h′   with  h′(β) = (β/2)(−μ + iη)
    const cuscomplex hp          = cusfloat(0.5) * beta * (-mu + J * eta);
    const cuscomplex df2a        = f2a * hp;

    const cuscomplex hbz1        = cusfloat(0.5) * beta * z1;
    const cuscomplex inv_hbz1_sq = cusfloat(1) / (hbz1 * hbz1);
    const cuscomplex j_over_eta  = J / eta;

    // Accumulate S0 = f2B and S1 = β·∂f2B/∂β
    cuscomplex r0_pow     { cusfloat(1), cusfloat(0) };
    cuscomplex hbz1_n_pow = hbz1;
    cuscomplex S0 { cusfloat(0), cusfloat(0) };
    cuscomplex S1 { cusfloat(0), cusfloat(0) };

    for (int n = 0; n < nmax; ++n) {
        cuscomplex hbz1_p = hbz1_n_pow;
        cuscomplex z2_pow { cusfloat(1), cusfloat(0) };
        cuscomplex r1     { cusfloat(0), cusfloat(0) };
        cuscomplex s1     { cusfloat(0), cusfloat(0) };
        cusfloat   p      = cusfloat(1) - cusfloat(2 * n);   // p = 1−2n at m=0

        for (int m = 0; m <= n; ++m) {
            const cuscomplex base = D(n, m) * hbz1_p * z2_pow;
            r1     += base;
            s1     += p * base;
            hbz1_p *= inv_hbz1_sq;
            z2_pow *= z2;
            p      -= cusfloat(2);
        }

        S0         += r0_pow * r1;
        S1         += r0_pow * s1;
        r0_pow     *= j_over_eta;
        hbz1_n_pow *= inv_hbz1_sq;
    }

    // d(f2A·f2B)/dβ = df2A·S0 + f2A·(S1/β)
    cusfloat result = (df2a * S0 + f2a * (S1 / beta)).real();

    // df0/dβ = 4 Σ α[n]·(2n+3)·β^{−2n−4}·P_n(μ)
    if (cfg.use_f0) {
        cusfloat P[TD_ASYMP_NMAX]{};
        cusfloat dP[TD_ASYMP_NMAX]{};
        compute_legendre(nmax, mu, eta * eta, P, dP);

        const cusfloat inv_beta2 = cusfloat(1) / (beta * beta);
        cusfloat beta_pow        = inv_beta2 * inv_beta2;   // β^{−4}
        cusfloat df0             = cusfloat(0);

        for (int n = 0; n < nmax; ++n) {
            df0      += A(n) * cusfloat(2*n + 3) * beta_pow * P[n];
            beta_pow *= inv_beta2;
        }
        result += cusfloat(4) * df0;
    }

    return result;
}


// ----------------------------------------------------------------------------
//  dGdttt  (∂²/∂β² of dGdt)
// ----------------------------------------------------------------------------

cusfloat dGdttt_asymptotic(cusfloat beta, cusfloat mu, TDAsympConfig cfg) {

    const int      nmax = cfg.nmax;
    const cusfloat eta  = std::sqrt(cusfloat(1) - mu * mu);
    const cuscomplex J  { cusfloat(0), cusfloat(1) };
    const cuscomplex z1 = J * mu + eta;
    const cuscomplex z2 = mu - J * eta;

    const cuscomplex f2a         = f2a_envelope(beta, mu, eta);
    // h′ = (β/2)(−μ+iη),  h″ = (1/2)(−μ+iη)
    const cuscomplex neg_mu_jeta = -mu + J * eta;
    const cuscomplex hp          = cusfloat(0.5) * beta * neg_mu_jeta;
    const cuscomplex hpp         = cusfloat(0.5)        * neg_mu_jeta;
    const cuscomplex df2a        = f2a * hp;
    const cuscomplex d2f2a       = f2a * (hp * hp + hpp);

    const cuscomplex hbz1        = cusfloat(0.5) * beta * z1;
    const cuscomplex inv_hbz1_sq = cusfloat(1) / (hbz1 * hbz1);
    const cuscomplex j_over_eta  = J / eta;

    // Accumulate S0, S1 = β·∂f2B/∂β, S2 = β²·∂²f2B/∂β²
    cuscomplex r0_pow     { cusfloat(1), cusfloat(0) };
    cuscomplex hbz1_n_pow = hbz1;
    cuscomplex S0 { cusfloat(0), cusfloat(0) };
    cuscomplex S1 { cusfloat(0), cusfloat(0) };
    cuscomplex S2 { cusfloat(0), cusfloat(0) };

    for (int n = 0; n < nmax; ++n) {
        cuscomplex hbz1_p = hbz1_n_pow;
        cuscomplex z2_pow { cusfloat(1), cusfloat(0) };
        cuscomplex s0     { cusfloat(0), cusfloat(0) };
        cuscomplex s1     { cusfloat(0), cusfloat(0) };
        cuscomplex s2     { cusfloat(0), cusfloat(0) };
        cusfloat   p      = cusfloat(1) - cusfloat(2 * n);

        for (int m = 0; m <= n; ++m) {
            const cuscomplex base = D(n, m) * hbz1_p * z2_pow;
            s0 += base;
            s1 += p               * base;
            s2 += p * (p - cusfloat(1)) * base;
            hbz1_p *= inv_hbz1_sq;
            z2_pow *= z2;
            p      -= cusfloat(2);
        }

        S0         += r0_pow * s0;
        S1         += r0_pow * s1;
        S2         += r0_pow * s2;
        r0_pow     *= j_over_eta;
        hbz1_n_pow *= inv_hbz1_sq;
    }

    // d²(f2A·f2B)/dβ² = d²f2A·S0 + 2·df2A·(S1/β) + f2A·(S2/β²)
    const cusfloat inv_beta  = cusfloat(1) / beta;
    const cusfloat inv_beta2 = inv_beta * inv_beta;
    const cuscomplex d2f2    = d2f2a * S0
                              + cusfloat(2) * df2a * (S1 * inv_beta)
                              + f2a               * (S2 * inv_beta2);

    cusfloat result = d2f2.real();

    // d²f0/dβ² = −4 Σ α[n]·(2n+3)(2n+4)·β^{−2n−5}·P_n(μ)
    if (cfg.use_f0) {
        cusfloat P[TD_ASYMP_NMAX]{};
        cusfloat dP[TD_ASYMP_NMAX]{};
        compute_legendre(nmax, mu, eta * eta, P, dP);

        const cusfloat inv_b2    = cusfloat(1) / (beta * beta);
        cusfloat       beta_pow  = inv_b2 * inv_b2 / beta;   // β^{−5}
        cusfloat       d2f0      = cusfloat(0);

        for (int n = 0; n < nmax; ++n) {
            d2f0     += A(n) * cusfloat((2*n+3)*(2*n+4)) * beta_pow * P[n];
            beta_pow *= inv_b2;
        }
        result += cusfloat(-4) * d2f0;
    }

    return result;
}


// ----------------------------------------------------------------------------
//  dGdtx  (∂/∂μ of dGdt)
// ----------------------------------------------------------------------------

cusfloat dGdtx_asymptotic(cusfloat beta, cusfloat mu, TDAsympConfig cfg) {

    const int      nmax = cfg.nmax;
    const cusfloat eta  = std::sqrt(cusfloat(1) - mu * mu);
    const cusfloat eta2 = eta * eta;
    const cuscomplex J  { cusfloat(0), cusfloat(1) };
    const cuscomplex z1 = J * mu + eta;
    const cuscomplex z2 = mu - J * eta;

    const cuscomplex f2a = f2a_envelope(beta, mu, eta);

    // ∂f2A/∂μ = f2A · c_μ
    //   c_μ = μ/(2η²) − β²/4 + i·(2 − β²μ)/(4η)
    const cuscomplex c_mu    = cusfloat(mu) / (cusfloat(2) * eta2)
                              - cusfloat(0.25) * beta * beta
                              + J * ((cusfloat(2) - beta * beta * mu) / (cusfloat(4) * eta));
    const cuscomplex df2a_mu = f2a * c_mu;

    const cuscomplex hbz1        = cusfloat(0.5) * beta * z1;
    const cuscomplex inv_hbz1_sq = cusfloat(1) / (hbz1 * hbz1);
    const cuscomplex j_over_eta  = J / eta;

    // Precompute  1/(η·z1)  and  1/(η·z2)  for the c_mn log-derivative
    const cuscomplex inv_eta_z1  = cusfloat(1) / (eta * z1);
    const cuscomplex inv_eta_z2  = cusfloat(1) / (eta * z2);
    const cusfloat   mu_eta2_inv = mu / eta2;    // μ/η²  (r0 log-deriv coefficient)

    // Accumulate S0 = f2B  and  dS0 = ∂f2B/∂μ
    cuscomplex r0_pow     { cusfloat(1), cusfloat(0) };
    cuscomplex hbz1_n_pow = hbz1;
    cuscomplex S0  { cusfloat(0), cusfloat(0) };
    cuscomplex dS0 { cusfloat(0), cusfloat(0) };   // ∂f2B/∂μ

    for (int n = 0; n < nmax; ++n) {
        cuscomplex hbz1_p = hbz1_n_pow;
        cuscomplex z2_pow { cusfloat(1), cusfloat(0) };
        cuscomplex r1     { cusfloat(0), cusfloat(0) };
        cuscomplex dr1    { cusfloat(0), cusfloat(0) };   // Σ_m base·c_mn (log-part only)
        cusfloat   p      = cusfloat(1) - cusfloat(2 * n);
        cusfloat   mf     = cusfloat(0);   // m as cusfloat

        for (int m = 0; m <= n; ++m) {
            const cuscomplex base = D(n, m) * hbz1_p * z2_pow;
            // c_mn = −p·z2·inv_eta_z1  +  m·z1·inv_eta_z2
            // (second term vanishes at m=0 because mf=0)
            const cuscomplex c_mn = -p * z2 * inv_eta_z1 + mf * z1 * inv_eta_z2;
            r1   += base;
            dr1  += base * c_mn;
            hbz1_p *= inv_hbz1_sq;
            z2_pow *= z2;
            p  -= cusfloat(2);
            mf += cusfloat(1);
        }

        // r0 log-deriv: d/dμ[(j/η)^n] contributes n·μ/η² · r1
        S0  += r0_pow * r1;
        dS0 += r0_pow * (cusfloat(n) * mu_eta2_inv * r1 + dr1);
        r0_pow     *= j_over_eta;
        hbz1_n_pow *= inv_hbz1_sq;
    }

    // ∂(f2A·f2B)/∂μ = ∂f2A/∂μ · S0 + f2A · ∂f2B/∂μ
    cusfloat result = (df2a_mu * S0 + f2a * dS0).real();

    // ∂f0/∂μ = −4 Σ α[n]·β^{−2n−3}·P_n′(μ)
    if (cfg.use_f0) {
        cusfloat P[TD_ASYMP_NMAX]{};
        cusfloat dP[TD_ASYMP_NMAX]{};
        compute_legendre(nmax, mu, eta2, P, dP);

        const cusfloat inv_beta2 = cusfloat(1) / (beta * beta);
        cusfloat beta_pow        = cusfloat(1) / (beta * beta * beta);  // β^{−3}
        cusfloat df0_dmu         = cusfloat(0);

        for (int n = 0; n < nmax; ++n) {
            df0_dmu  += A(n) * beta_pow * dP[n];
            beta_pow *= inv_beta2;
        }
        result += cusfloat(-4) * df0_dmu;
    }

    return result;
}


// ----------------------------------------------------------------------------
//  dGdttx  (∂²(dGdt) / ∂β ∂μ)
// ----------------------------------------------------------------------------

cusfloat dGdttx_asymptotic(cusfloat beta, cusfloat mu, TDAsympConfig cfg) {

    const int      nmax = cfg.nmax;
    const cusfloat eta  = std::sqrt(cusfloat(1) - mu * mu);
    const cusfloat eta2 = eta * eta;
    const cuscomplex J  { cusfloat(0), cusfloat(1) };
    const cuscomplex z1 = J * mu + eta;
    const cuscomplex z2 = mu - J * eta;

    const cuscomplex f2a = f2a_envelope(beta, mu, eta);

    // β-derivative of f2A:   f2A_β = f2A · h′,   h′ = (β/2)(−μ+iη)
    const cuscomplex hp      = cusfloat(0.5) * beta * (-mu + J * eta);
    const cuscomplex f2a_b   = f2a * hp;

    // μ-derivative of f2A:   f2A_μ = f2A · c_μ
    const cuscomplex c_mu    = cusfloat(mu) / (cusfloat(2) * eta2)
                              - cusfloat(0.25) * beta * beta
                              + J * ((cusfloat(2) - beta * beta * mu) / (cusfloat(4) * eta));
    const cuscomplex f2a_mu  = f2a * c_mu;

    // Mixed (β,μ)-derivative of f2A:
    //   h′_μ = ∂h′/∂μ = −β·z1/(2η)
    //   f2A_βμ = f2A · (c_μ·h′ + h′_μ)
    const cuscomplex hp_mu   = cusfloat(-0.5) * beta * z1 / eta;
    const cuscomplex f2a_bmu = f2a * (c_mu * hp + hp_mu);

    const cuscomplex hbz1        = cusfloat(0.5) * beta * z1;
    const cuscomplex inv_hbz1_sq = cusfloat(1) / (hbz1 * hbz1);
    const cuscomplex j_over_eta  = J / eta;

    const cuscomplex inv_eta_z1  = cusfloat(1) / (eta * z1);
    const cuscomplex inv_eta_z2  = cusfloat(1) / (eta * z2);
    const cusfloat   mu_eta2_inv = mu / eta2;

    // Accumulate all four inner sums in a single pass
    cuscomplex r0_pow     { cusfloat(1), cusfloat(0) };
    cuscomplex hbz1_n_pow = hbz1;
    cuscomplex S0  { cusfloat(0), cusfloat(0) };   // f2B
    cuscomplex dS0 { cusfloat(0), cusfloat(0) };   // ∂f2B/∂μ
    cuscomplex S1  { cusfloat(0), cusfloat(0) };   // β·∂f2B/∂β
    cuscomplex dS1 { cusfloat(0), cusfloat(0) };   // β·∂²f2B/(∂β∂μ)

    for (int n = 0; n < nmax; ++n) {
        cuscomplex hbz1_p = hbz1_n_pow;
        cuscomplex z2_pow { cusfloat(1), cusfloat(0) };
        cuscomplex r1     { cusfloat(0), cusfloat(0) };
        cuscomplex dr1    { cusfloat(0), cusfloat(0) };
        cuscomplex s1     { cusfloat(0), cusfloat(0) };
        cuscomplex ds1    { cusfloat(0), cusfloat(0) };
        cusfloat   p      = cusfloat(1) - cusfloat(2 * n);
        cusfloat   mf     = cusfloat(0);

        for (int m = 0; m <= n; ++m) {
            const cuscomplex base = D(n, m) * hbz1_p * z2_pow;
            const cuscomplex c_mn = -p * z2 * inv_eta_z1 + mf * z1 * inv_eta_z2;
            r1  += base;
            dr1 += base * c_mn;
            s1  += p * base;
            ds1 += p * base * c_mn;
            hbz1_p *= inv_hbz1_sq;
            z2_pow *= z2;
            p  -= cusfloat(2);
            mf += cusfloat(1);
        }

        const cusfloat n_mu_eta2 = cusfloat(n) * mu_eta2_inv;
        S0  += r0_pow * r1;
        dS0 += r0_pow * (n_mu_eta2 * r1 + dr1);
        S1  += r0_pow * s1;
        dS1 += r0_pow * (n_mu_eta2 * s1 + ds1);
        r0_pow     *= j_over_eta;
        hbz1_n_pow *= inv_hbz1_sq;
    }

    // ∂²(f2A·f2B)/(∂β∂μ) = f2A_βμ·S0 + f2A_β·dS0 + f2A_μ·(S1/β) + f2A·(dS1/β)
    const cusfloat   inv_beta    = cusfloat(1) / beta;
    const cuscomplex d2f2_mixed  = f2a_bmu * S0
                                  + f2a_b  * dS0
                                  + f2a_mu * (S1  * inv_beta)
                                  + f2a    * (dS1 * inv_beta);

    cusfloat result = d2f2_mixed.real();

    // ∂²f0/(∂β∂μ) = 4 Σ α[n]·(2n+3)·β^{−2n−4}·P_n′(μ)
    if (cfg.use_f0) {
        cusfloat P[TD_ASYMP_NMAX]{};
        cusfloat dP[TD_ASYMP_NMAX]{};
        compute_legendre(nmax, mu, eta2, P, dP);

        const cusfloat inv_b2  = cusfloat(1) / (beta * beta);
        cusfloat beta_pow      = inv_b2 * inv_b2;   // β^{−4}
        cusfloat d2f0_dxdt     = cusfloat(0);

        for (int n = 0; n < nmax; ++n) {
            d2f0_dxdt += A(n) * cusfloat(2*n + 3) * beta_pow * dP[n];
            beta_pow  *= inv_b2;
        }
        result += cusfloat(4) * d2f0_dxdt;
    }

    return result;
}
