
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

#pragma once

#include "../config.hpp"

// ---------------------------------------------------------------------------
//  Truncation-order constants
// ---------------------------------------------------------------------------

/// Maximum inner-series truncation order supported by the tables.
constexpr int TD_ASYMP_NMAX = 15;

/// Reduced truncation order validated for β ≥ 55
/// (see test_high_beta_simplifications in aux_tools/test_dGdt_asymptotic.py).
constexpr int TD_ASYMP_NMAX_FAST = 6;

// ---------------------------------------------------------------------------
//  Configuration struct
// ---------------------------------------------------------------------------

/**
 * @brief Configuration for the time-domain asymptotic Green function routines.
 *
 * Each function computes  g(β, μ) = f0(β, μ) + Re[ f2(β, μ) ]  where:
 *
 *   f0   — non-oscillatory power-series correction in 1/β.
 *   f2   — stationary-phase oscillatory contribution (dominant for large β).
 *
 * Two modes are provided:
 *
 *   Full (default):    use_f0 = true,  nmax = TD_ASYMP_NMAX (= 15).
 *
 *   Optimised (fast):  use_f0 = false, nmax = TD_ASYMP_NMAX_FAST (= 6).
 *     Valid for β ≥ 55.  For all five functions:
 *       |f0| / (|f0| + |f2|) < 1e-3,
 *       |f2(nmax=6) − f2(nmax=15)| / |f2(nmax=15)| < 1e-8.
 */
struct TDAsympConfig {
    bool use_f0 = true;           ///< Include f0 power-series correction.
    int  nmax   = TD_ASYMP_NMAX;  ///< Inner-series truncation order (1 … TD_ASYMP_NMAX).
};

/// Ready-made preset for the optimised high-β mode (β ≥ 55).
constexpr TDAsympConfig TD_ASYMP_FAST{ false, TD_ASYMP_NMAX_FAST };

// ---------------------------------------------------------------------------
//  Public interface
// ---------------------------------------------------------------------------

/**
 * @brief Asymptotic approximation of dG/dt.
 *
 * Asymptotic expansion of the time derivative of the transient free-surface
 * Green function, valid for large non-dimensional time β.
 *
 * @param beta  Non-dimensional time β = (t − τ) g / U.
 * @param mu    Angular integration parameter μ = cos θ ∈ (−1, 1).
 * @param cfg   Computation mode (default: full).
 * @return      Real-valued dG/dt.
 */
cusfloat dGdt_asymptotic  (cusfloat beta, cusfloat mu, TDAsympConfig cfg = {});

/**
 * @brief Asymptotic approximation of d²G/dt²  (∂/∂β of dG/dt).
 *
 * @param beta  Non-dimensional time β.
 * @param mu    Angular integration parameter μ ∈ (−1, 1).
 * @param cfg   Computation mode (default: full).
 * @return      Real-valued d²G/dt².
 */
cusfloat dGdtt_asymptotic (cusfloat beta, cusfloat mu, TDAsympConfig cfg = {});

/**
 * @brief Asymptotic approximation of d³G/dt³  (∂²/∂β² of dG/dt).
 *
 * @param beta  Non-dimensional time β.
 * @param mu    Angular integration parameter μ ∈ (−1, 1).
 * @param cfg   Computation mode (default: full).
 * @return      Real-valued d³G/dt³.
 */
cusfloat dGdttt_asymptotic(cusfloat beta, cusfloat mu, TDAsympConfig cfg = {});

/**
 * @brief Asymptotic approximation of ∂(dG/dt)/∂μ.
 *
 * @param beta  Non-dimensional time β.
 * @param mu    Angular integration parameter μ ∈ (−1, 1).
 * @param cfg   Computation mode (default: full).
 * @return      Real-valued ∂(dG/dt)/∂μ.
 */
cusfloat dGdtx_asymptotic (cusfloat beta, cusfloat mu, TDAsympConfig cfg = {});

/**
 * @brief Asymptotic approximation of ∂²(dG/dt)/(∂β ∂μ).
 *
 * @param beta  Non-dimensional time β.
 * @param mu    Angular integration parameter μ ∈ (−1, 1).
 * @param cfg   Computation mode (default: full).
 * @return      Real-valued ∂²(dG/dt)/(∂β ∂μ).
 */
cusfloat dGdttx_asymptotic(cusfloat beta, cusfloat mu, TDAsympConfig cfg = {});
