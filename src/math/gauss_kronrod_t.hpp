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

// Include local modules
#include "../config.hpp"

/**
 * @brief Gauss-Kronrod quadrature nodes and weights on the reference interval [-1, 1].
 *
 * GaussKronrodT<NGK> provides the NGK-point Kronrod rule that embeds the
 * standard Gauss-Legendre rule of order (NGK-1)/2.  Supported sizes:
 *
 *   NGK =  7  →  G3K7   (embeds 3-point GL; exact for polynomials of degree ≤ 13)
 *   NGK = 15  →  G7K15  (embeds 7-point GL; exact for polynomials of degree ≤ 29)
 *   NGK = 21  →  G10K21 (embeds 10-point GL; exact for polynomials of degree ≤ 41)
 *
 * Nodes are sorted in ascending order.  The GL nodes are a strict subset of
 * the Kronrod nodes (marked in comments), enabling both the high-order
 * Kronrod estimate and the embedded Gauss estimate from the same function
 * evaluations when needed for error control.
 *
 * Reference values are taken from QUADPACK (Piessens et al., 1983) and are
 * accurate to machine precision in double precision.
 */
template<int NGK>
struct GaussKronrodT
{
    static_assert( 
                    NGK == 15 || NGK == 21,
                   "GaussKronrodT is only specialised for NGK = 15 or 21." 
                );
};


// =============================================================================
// G7K15  —  15-point Kronrod rule   (embeds the 7-point Gauss-Legendre rule)
// Exact for polynomials of degree <= 29.
// GL7 nodes (marked with ←): ±0.9491, ±0.7415, ±0.4058, 0
// =============================================================================
template<>
struct GaussKronrodT<15>
{
public:
    static constexpr int        N           = 15;

    // Kronrod nodes on [-1, 1], sorted ascending
    static constexpr cusfloat roots_x[15] =
    {
        -9.914553711208126384E-01,   // new
        -9.491079123427585245E-01,   // ← GL7
        -8.648644233597498182E-01,   // new
        -7.415311855993944398E-01,   // ← GL7
        -5.860872354676911185E-01,   // new
        -4.058451513773971669E-01,   // ← GL7
        -2.077849550078984676E-01,   // new
        +0.000000000000000000E+00,   // ← GL7
        +2.077849550078984676E-01,   // new
        +4.058451513773971669E-01,   // ← GL7
        +5.860872354676911185E-01,   // new
        +7.415311855993944398E-01,   // ← GL7
        +8.648644233597498182E-01,   // new
        +9.491079123427585245E-01,   // ← GL7
        +9.914553711208126384E-01    // new
    };

    // Kronrod weights (corresponding to roots_x)
    static constexpr cusfloat weights_x[15] =
    {
        +2.293532201052922E-02,
        +6.309209262828563E-02,
        +1.047900103222502E-01,
        +1.406532597155259E-01,
        +1.690047266392679E-01,
        +1.903505780647854E-01,
        +2.044329400752989E-01,
        +2.094821410847278E-01,
        +2.044329400752989E-01,
        +1.903505780647854E-01,
        +1.690047266392679E-01,
        +1.406532597155259E-01,
        +1.047900103222502E-01,
        +6.309209262828563E-02,
        +2.293532201052922E-02
    };
};


// =============================================================================
// G10K21  —  21-point Kronrod rule   (embeds the 10-point Gauss-Legendre rule)
// Exact for polynomials of degree <= 41.
// GL10 nodes (marked with ←): ±0.9739, ±0.8651, ±0.6794, ±0.4334, ±0.1489
// =============================================================================
template<>
struct GaussKronrodT<21>
{
public:
    static constexpr int        N           = 21;

    // Kronrod nodes on [-1, 1], sorted ascending
    static constexpr cusfloat roots_x[21] =
    {
        -9.956571630258080807E-01,   // new
        -9.739065285171717201E-01,   // ← GL10
        -9.301574913557082261E-01,   // new
        -8.650633666889845107E-01,   // ← GL10
        -7.808177265864168970E-01,   // new
        -6.794095682990244062E-01,   // ← GL10
        -5.627571346686046605E-01,   // new
        -4.333953941292471907E-01,   // ← GL10
        -2.943928627014601981E-01,   // new
        -1.488743389816312109E-01,   // ← GL10
        +0.000000000000000000E+00,   // new (centre)
        +1.488743389816312109E-01,   // ← GL10
        +2.943928627014601981E-01,   // new
        +4.333953941292471907E-01,   // ← GL10
        +5.627571346686046605E-01,   // new
        +6.794095682990244062E-01,   // ← GL10
        +7.808177265864168970E-01,   // new
        +8.650633666889845107E-01,   // ← GL10
        +9.301574913557082261E-01,   // new
        +9.739065285171717201E-01,   // ← GL10
        +9.956571630258080807E-01    // new
    };

    // Kronrod weights (corresponding to roots_x)
    static constexpr cusfloat weights_x[21] =
    {
        +1.169463886737187E-02,
        +3.255816230796473E-02,
        +5.475589657435200E-02,
        +7.503967481091995E-02,
        +9.312545458369761E-02,
        +1.093871588022976E-01,
        +1.234919762620659E-01,
        +1.347092173181476E-01,
        +1.427759385770601E-01,
        +1.477391049013385E-01,
        +1.494455540029169E-01,
        +1.477391049013385E-01,
        +1.427759385770601E-01,
        +1.347092173181476E-01,
        +1.234919762620659E-01,
        +1.093871588022976E-01,
        +9.312545458369761E-02,
        +7.503967481091995E-02,
        +5.475589657435200E-02,
        +3.255816230796473E-02,
        +1.169463886737187E-02
    };
};
