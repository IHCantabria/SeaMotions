
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

#define CHEBYSHEV_1D_TRAITS( NS )                                                                               \
template<>                                                                                                      \
struct ChebyshevTraits<NS>                                                                                      \
{                                                                                                               \
    static constexpr            int             max_ref_level               = NS::max_ref_level;                \
    static constexpr            int             intervals_np                = NS::intervals_np;                 \
    static constexpr            int             max_cheby_order             = NS::max_cheby_order;              \
    static constexpr const      std::size_t*    blocks_start                = NS::blocks_start;                 \
    static constexpr const      std::size_t*    blocks_start_fall           = NS::blocks_start_fall;            \
    static constexpr const      std::size_t*    blocks_coeffs_np            = NS::blocks_coeffs_np;             \
    static constexpr const      std::size_t*    blocks_coeffs_np_fall       = NS::blocks_coeffs_np_fall;        \
    static constexpr const      std::size_t*    blocks_max_cheby_order      = NS::blocks_max_cheby_order;       \
    static constexpr const      std::size_t*    blocks_max_cheby_order_fall = NS::blocks_max_cheby_order_fall;  \
    static constexpr            bool            fcn_log_scale               = NS::fcn_log_scale;                \
    static constexpr            bool            x_log_scale                 = NS::x_log_scale;                  \
    static constexpr            cusfloat        x_max_global                = NS::x_max_global;                 \
    static constexpr            cusfloat        x_min_global                = NS::x_min_global;                 \
    static constexpr            cusfloat        dx_min_region               = NS::dx_min_region;                \
    static constexpr const      cusfloat*       x_max_region                = NS::x_max_region;                 \
    static constexpr const      cusfloat*       x_min_region                = NS::x_min_region;                 \
    static constexpr const      cusfloat*       dx_region                   = NS::dx_region;                    \
    static constexpr            std::size_t     num_c                       = NS::num_c;                        \
    static constexpr const      cusfloat*       coeffs                      = NS::c;                            \
    static constexpr const      std::size_t*    ncx                         = NS::ncx;                          \
};                                                                                                              \


#define CHEBYSHEV_1DF_TRAITS( NS )                                                                              \
template<>                                                                                                      \
struct ChebyshevTraits<NS>                                                                                      \
{                                                                                                               \
    inline static               cusfloat        coeffs                     = 0.0;                              \
};                                                                                                              \


#define CHEBYSHEV_2D_TRAITS( NS )                                                                               \
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
    static constexpr            std::size_t     num_c                       = NS::num_c;                        \
    static constexpr const      cusfloat*       coeffs                      = NS::c;                            \
    static constexpr const      std::size_t*    ncx                         = NS::ncx;                          \
    static constexpr const      std::size_t*    ncy                         = NS::ncy;                          \
};                                                                                                              \


#define CHEBYSHEV_3D_TRAITS( NS )                                                                               \
template<>                                                                                                      \
struct ChebyshevTraits<NS>                                                                                      \
{                                                                                                               \
    static constexpr            int             max_ref_level               = NS::max_ref_level;                \
    static constexpr            int             intervals_np                = NS::intervals_np;                 \
    static constexpr            int             max_cheby_order             = NS::max_cheby_order;              \
    static constexpr const      std::size_t*    blocks_start                = NS::blocks_start;                 \
    static constexpr const      std::size_t*    blocks_start_fall           = NS::blocks_start_fall;            \
    static constexpr const      std::size_t*    blocks_coeffs_np            = NS::blocks_coeffs_np;             \
    static constexpr const      std::size_t*    blocks_coeffs_np_fall       = NS::blocks_coeffs_np_fall;        \
    static constexpr const      std::size_t*    blocks_max_cheby_order      = NS::blocks_max_cheby_order;       \
    static constexpr const      std::size_t*    blocks_max_cheby_order_fall = NS::blocks_max_cheby_order_fall;  \
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
    static constexpr            bool            z_log_scale                 = NS::z_log_scale;                  \
    static constexpr            cusfloat        z_max_global                = NS::z_max_global;                 \
    static constexpr            cusfloat        z_min_global                = NS::z_min_global;                 \
    static constexpr            cusfloat        dz_min_region               = NS::dz_min_region;                \
    static constexpr const      cusfloat*       z_max_region                = NS::z_max_region;                 \
    static constexpr const      cusfloat*       z_min_region                = NS::z_min_region;                 \
    static constexpr const      cusfloat*       dz_region                   = NS::dz_region;                    \
    static constexpr            std::size_t     num_c                       = NS::num_c;                        \
    static constexpr const      cusfloat*       coeffs                      = NS::c;                            \
    static constexpr const      std::size_t*    ncx                         = NS::ncx;                          \
    static constexpr const      std::size_t*    ncy                         = NS::ncy;                          \
    static constexpr const      std::size_t*    ncz                         = NS::ncz;                          \
};                                                                                                              \


#define CHEBYSHEV_3DF_TRAITS( NS )                                                                                      \
template<>                                                                                                              \
struct ChebyshevTraits<NS>                                                                                              \
{                                                                                                                       \
            static constexpr            int             max_ref_level               = NS::max_ref_level;                \
            static constexpr            int             intervals_np                = NS::intervals_np;                 \
            static constexpr            int             max_cheby_order             = NS::max_cheby_order_f;            \
    inline  static                      std::size_t*    blocks_start                = NS::blocks_start_f;               \
    inline  static                      std::size_t*    blocks_coeffs_np            = NS::blocks_coeffs_np_f;           \
    inline  static                      std::size_t*    blocks_max_cheby_order      = NS::blocks_max_cheby_order_f;     \
            static constexpr            bool            fcn_log_scale               = NS::fcn_log_scale;                \
            static constexpr            bool            x_log_scale                 = NS::x_log_scale;                  \
            static constexpr            cusfloat        x_max_global                = NS::x_max_global;                 \
            static constexpr            cusfloat        x_min_global                = NS::x_min_global;                 \
            static constexpr            cusfloat        dx_min_region               = NS::dx_min_region;                \
    inline  static                      cusfloat*       x_max_region                = NS::x_max_region_f;               \
    inline  static                      cusfloat*       x_min_region                = NS::x_min_region_f;               \
    inline  static                      cusfloat*       dx_region                   = NS::dx_region_f;                  \
            static constexpr            bool            y_log_scale                 = NS::y_log_scale;                  \
            static constexpr            cusfloat        y_max_global                = NS::y_max_global;                 \
            static constexpr            cusfloat        y_min_global                = NS::y_min_global;                 \
            static constexpr            cusfloat        dy_min_region               = NS::dy_min_region;                \
    inline  static                      cusfloat*       y_max_region                = NS::y_max_region_f;               \
    inline  static                      cusfloat*       y_min_region                = NS::y_min_region_f;               \
    inline  static                      cusfloat*       dy_region                   = NS::dy_region_f;                  \
            static constexpr            std::size_t     num_c                       = NS::num_cf;                       \
    inline  static                      cusfloat*       coeffs                      = NS::cf;                           \
    inline  static                      std::size_t*    ncx                         = NS::ncxf;                         \
    inline  static                      std::size_t*    ncy                         = NS::ncyf;                         \
};                                                                                                                      \
