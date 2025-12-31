
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

#include "../../config.hpp"

struct L2C
{
static constexpr int                                        max_ref_level = 1;
static constexpr int                                        intervals_np = 2;
static constexpr int                                        max_cheby_order = 11;
static constexpr int                                        max_cheby_order_f = 11;
static constexpr std::size_t                                blocks_start[2] = {0, 6};
static constexpr std::size_t                                blocks_start_fall[2] = {0, 0};
static constexpr std::size_t                                blocks_coeffs_np[2] = {6, 12};
static constexpr std::size_t                                blocks_coeffs_np_fall[2] = {1, 1};
static constexpr std::size_t                                blocks_max_cheby_order[2] = {5, 11};
static constexpr std::size_t                                blocks_max_cheby_order_fall[2] = {5, 11};
static constexpr cusfloat                                   x_min_region[2] = {-5.000000000000000E+00, -2.500000000000000E+00};
static constexpr cusfloat                                   x_max_region[2] = {-2.500000000000000E+00, 0.000000000000000E+00};
static constexpr cusfloat                                   dx_region[2] = {2.500000000000000E+00, 2.500000000000000E+00};
static constexpr bool                                       fcn_log_scale = 0;
static constexpr bool                                       x_log_scale = 1;
static constexpr cusfloat                                   x_max_global = 0.0000000000000000;
static constexpr cusfloat                                   x_min_global = -5.0000000000000000;
static constexpr cusfloat                                   dx_min_region = 2.5000000000000000;

static constexpr std::size_t                                num_c = 18;
alignas(FLOATING_PRECISION)  static constexpr cusfloat      c[18] = {
                                                                2.6230507433814698E+00,  // c[0]
                                                                -1.4408945282408883E+00,  // c[1]
                                                                -9.1645527450229247E-04,  // c[2]
                                                                -3.4373381262338409E-04,  // c[3]
                                                                -9.9580403236487758E-05,  // c[4]
                                                                -2.3144144640857567E-05,  // c[5]
                                                                -3.7275006307808456E-01,  // c[6]
                                                                -1.6077276829592893E+00,  // c[7]
                                                                -6.7211819034913578E-02,  // c[8]
                                                                -1.2086015630786526E-02,  // c[9]
                                                                3.5480756226971832E-03,  // c[10]
                                                                4.2837323006755246E-03,  // c[11]
                                                                2.4295209811072277E-03,  // c[12]
                                                                1.0962945315780015E-03,  // c[13]
                                                                4.3635068246852815E-04,  // c[14]
                                                                1.5628005268397693E-04,  // c[15]
                                                                4.9233646028018496E-05,  // c[16]
                                                                1.2497401632990957E-05,  // c[17]
                                                  };
alignas(FLOATING_PRECISION)  static constexpr std::size_t   ncx[18] = {
                                                                0,  // ncx[0]
                                                                1,  // ncx[1]
                                                                2,  // ncx[2]
                                                                3,  // ncx[3]
                                                                4,  // ncx[4]
                                                                5,  // ncx[5]
                                                                0,  // ncx[6]
                                                                1,  // ncx[7]
                                                                2,  // ncx[8]
                                                                3,  // ncx[9]
                                                                4,  // ncx[10]
                                                                5,  // ncx[11]
                                                                6,  // ncx[12]
                                                                7,  // ncx[13]
                                                                8,  // ncx[14]
                                                                9,  // ncx[15]
                                                                10,  // ncx[16]
                                                                11,  // ncx[17]
                                                  };
};

