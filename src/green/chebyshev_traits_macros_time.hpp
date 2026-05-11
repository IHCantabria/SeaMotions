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

// Macro for time-domain 2D chebyshev coefficients (read-only 2D view).
// Carries G0 near-field correction data and full 2D coefficient arrays.
#define CHEBYSHEV_TIME_2D_TRAITS( NS )                                                                          \
template<>                                                                                                      \
struct ChebyshevTraits<NS>                                                                                      \
{                                                                                                               \
    static constexpr            int             max_ref_level               = NS::max_ref_level;                \
    static constexpr            int             intervals_np                = NS::intervals_np;                 \
    static constexpr            int             max_cheby_order             = NS::max_cheby_order;              \
    static constexpr const      std::size_t*    blocks_start                = NS::blocks_start;                 \
    static constexpr const      std::size_t*    blocks_coeffs_np            = NS::blocks_coeffs_np;             \
    static constexpr const      std::size_t*    blocks_max_cheby_order      = NS::blocks_max_cheby_order;       \
    static constexpr            bool            fcn_log_scale               = NS::fcn_log_scale;                \
    static constexpr            bool            x_log_scale                 = NS::x_log_scale;                  \
    static constexpr            cusfloat        x_max_global                = NS::x_max_global;                 \
    static constexpr            cusfloat        x_min_global                = NS::x_min_global;                 \
    static constexpr            cusfloat        dx_min_region               = NS::dx_min_region;                \
    static constexpr const      cusfloat*       x_max_region                = NS::x_max_region;                 \
    static constexpr const      cusfloat*       x_min_region                = NS::x_min_region;                 \
    static constexpr const      cusfloat*       dx_region                   = NS::dx_region;                    \
    static constexpr            bool            y_log_scale                 = NS::y_log_scale;                  \
    static constexpr            cusfloat        y_max_global                = NS::y_max_global;                 \
    static constexpr            cusfloat        y_min_global                = NS::y_min_global;                 \
    static constexpr            cusfloat        dy_min_region               = NS::dy_min_region;                \
    static constexpr const      cusfloat*       y_max_region                = NS::y_max_region;                 \
    static constexpr const      cusfloat*       y_min_region                = NS::y_min_region;                 \
    static constexpr const      cusfloat*       dy_region                   = NS::dy_region;                    \
    static constexpr            int             G0_alpha_shift_np           = NS::G0_alpha_shift_np;            \
    static constexpr const      cusfloat*       G0_alpha_shift              = NS::G0_alpha_shift;               \
    static constexpr const      cusfloat*       G0_cheby_intv               = NS::G0_cheby_intv;                \
    static constexpr            int             G0_cheby_np                 = NS::G0_cheby_np;                  \
    static constexpr const      cusfloat*       G0_cheby_mui_log            = NS::G0_cheby_mui_log;             \
    static constexpr const      cusfloat*       G0_cheby_mui                = NS::G0_cheby_mui;                 \
    static constexpr const      cusfloat*       G0_coeffs                   = NS::G0_coeffs;                    \
    static constexpr            cusfloat        G0_norm_f                   = NS::G0_norm_f;                    \
    static constexpr            std::size_t     num_c                       = NS::num_c;                        \
    static constexpr const      cusfloat*       coeffs                      = NS::c;                            \
    static constexpr const      std::size_t*    ncx                         = NS::ncx;                          \
    static constexpr const      std::size_t*    ncy                         = NS::ncy;                          \
};                                                                                                              \


// Macro for time-domain 2D Chebyshev coefficients that support fold-X (beta).
// Provides the same read-only 2D metadata as CHEBYSHEV_TIME_2D_TRAITS plus
// mutable inline-static storage for the 1D (Y = log_mu) folded result after
// calling fold_time_residual_x<NS>(beta).
//
// Size helpers:
//   nx_blocks_f  - number of X (beta) patches
//   ny_blocks_f  - number of Y (log_mu) patches
//   max_f_coeffs - conservative upper bound on folded coefficient count
//
// After fold_time_residual_x<NS>(beta) is called the following members are
// populated and can be consumed by ChebyshevEvaluatorBase1D /
// ChebyshevEvaluatorBase1DVector:
//   coeffs_f, ncy_f          - folded coefficient values and Y-orders
//   blocks_start_f           - start index in coeffs_f for each Y-patch
//   blocks_np_f              - number of folded terms per Y-patch
//   blocks_max_order_f       - maximum Y-order per Y-patch
//   y_min/max/dy_region_f    - Y-region bounds for mapping to [-1,1]
#define CHEBYSHEV_TIME_2DF_TRAITS( NS )                                                                                         \
template<>                                                                                                                      \
struct ChebyshevTraits<NS>                                                                                                      \
{                                                                                                                               \
    /* ---- size helpers ---- */                                                                                                \
    static constexpr std::size_t    nx_blocks_f                 = static_cast<std::size_t>( ( NS::x_max_global - NS::x_min_global ) / NS::dx_min_region + 0.5f );      \
    static constexpr std::size_t    ny_blocks_f                 = static_cast<std::size_t>( ( NS::y_max_global - NS::y_min_global ) / NS::dy_min_region + 0.5f );  \
    static constexpr std::size_t    max_f_coeffs                = ny_blocks_f * ( NS::max_cheby_order + 1 );                   \
    /* ---- read-only metadata (same as CHEBYSHEV_TIME_2D_TRAITS) ---- */                                                      \
    static constexpr            int             max_ref_level               = NS::max_ref_level;                               \
    static constexpr            int             intervals_np                = NS::intervals_np;                                \
    static constexpr            int             max_cheby_order             = NS::max_cheby_order;                             \
    static constexpr const      std::size_t*    blocks_start                = NS::blocks_start;                                \
    static constexpr const      std::size_t*    blocks_coeffs_np            = NS::blocks_coeffs_np;                            \
    static constexpr const      std::size_t*    blocks_max_cheby_order      = NS::blocks_max_cheby_order;                      \
    static constexpr            bool            fcn_log_scale               = NS::fcn_log_scale;                               \
    static constexpr            bool            x_log_scale                 = NS::x_log_scale;                                 \
    static constexpr            cusfloat        x_max_global                = NS::x_max_global;                                \
    static constexpr            cusfloat        x_min_global                = NS::x_min_global;                                \
    static constexpr            cusfloat        dx_min_region               = NS::dx_min_region;                               \
    static constexpr const      cusfloat*       x_max_region                = NS::x_max_region;                                \
    static constexpr const      cusfloat*       x_min_region                = NS::x_min_region;                                \
    static constexpr const      cusfloat*       dx_region                   = NS::dx_region;                                   \
    static constexpr            bool            y_log_scale                 = NS::y_log_scale;                                 \
    static constexpr            cusfloat        y_max_global                = NS::y_max_global;                                \
    static constexpr            cusfloat        y_min_global                = NS::y_min_global;                                \
    static constexpr            cusfloat        dy_min_region               = NS::dy_min_region;                               \
    static constexpr const      cusfloat*       y_max_region                = NS::y_max_region;                                \
    static constexpr const      cusfloat*       y_min_region                = NS::y_min_region;                                \
    static constexpr const      cusfloat*       dy_region                   = NS::dy_region;                                   \
    static constexpr            int             G0_alpha_shift_np           = NS::G0_alpha_shift_np;                           \
    static constexpr const      cusfloat*       G0_alpha_shift              = NS::G0_alpha_shift;                              \
    static constexpr const      cusfloat*       G0_cheby_intv               = NS::G0_cheby_intv;                               \
    static constexpr            int             G0_cheby_np                 = NS::G0_cheby_np;                                 \
    static constexpr const      cusfloat*       G0_cheby_mui_log            = NS::G0_cheby_mui_log;                            \
    static constexpr const      cusfloat*       G0_cheby_mui                = NS::G0_cheby_mui;                                \
    static constexpr const      cusfloat*       G0_coeffs                   = NS::G0_coeffs;                                   \
    static constexpr            cusfloat        G0_norm_f                   = NS::G0_norm_f;                                   \
    static constexpr            std::size_t     num_c                       = NS::num_c;                                       \
    static constexpr const      cusfloat*       coeffs                      = NS::c;                                           \
    static constexpr const      std::size_t*    ncx                         = NS::ncx;                                         \
    static constexpr const      std::size_t*    ncy                         = NS::ncy;                                         \
    /* ---- mutable fold-X storage: populated by fold_time_residual_x ---- */                                                  \
    inline static               cusfloat        coeffs_f[max_f_coeffs]          = {};                                          \
    inline static               std::size_t     ncy_f[max_f_coeffs]             = {};                                          \
    inline static               std::size_t     blocks_start_f[ny_blocks_f]     = {};                                          \
    inline static               std::size_t     blocks_np_f[ny_blocks_f]        = {};                                          \
    inline static               std::size_t     blocks_max_order_f[ny_blocks_f] = {};                                          \
    inline static               cusfloat        y_min_region_f[ny_blocks_f]     = {};                                          \
    inline static               cusfloat        y_max_region_f[ny_blocks_f]     = {};                                          \
    inline static               cusfloat        dy_region_f[ny_blocks_f]        = {};                                          \
};                                                                                                                              \
