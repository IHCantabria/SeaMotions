
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
static constexpr int                                        max_cheby_order = 15;
static constexpr int                                        max_cheby_order_f = 15;
static constexpr std::size_t                                blocks_start[2] = {0, 8};
static constexpr std::size_t                                blocks_start_fall[2] = {0, 0};
static constexpr std::size_t                                blocks_coeffs_np[2] = {8, 16};
static constexpr std::size_t                                blocks_coeffs_np_fall[2] = {1, 1};
static constexpr std::size_t                                blocks_max_cheby_order[2] = {7, 15};
static constexpr std::size_t                                blocks_max_cheby_order_fall[2] = {7, 15};
static constexpr cusfloat                                   x_min_region[2] = {-5.000000000000000E+00, -2.500000000000000E+00};
static constexpr cusfloat                                   x_max_region[2] = {-2.500000000000000E+00, 0.000000000000000E+00};
static constexpr cusfloat                                   dx_region[2] = {2.500000000000000E+00, 2.500000000000000E+00};
static constexpr bool                                       fcn_log_scale = 0;
static constexpr bool                                       x_log_scale = 1;
static constexpr cusfloat                                   x_max_global = 0.0000000000000000;
static constexpr cusfloat                                   x_min_global = -5.0000000000000000;
static constexpr cusfloat                                   dx_min_region = 2.5000000000000000;

static constexpr std::size_t                                num_c = 24;
alignas(FLOATING_PRECISION)  static constexpr cusfloat      c[24] = {
                                                                2.6230507433814698E+00,  // c[0]
                                                                -1.4408945282408883E+00,  // c[1]
                                                                -9.1645527450229247E-04,  // c[2]
                                                                -3.4373381262338409E-04,  // c[3]
                                                                -9.9580403236487758E-05,  // c[4]
                                                                -2.3144144640857567E-05,  // c[5]
                                                                -4.4150356670602520E-06,  // c[6]
                                                                -6.9792671425483732E-07,  // c[7]
                                                                -3.7275006307808456E-01,  // c[8]
                                                                -1.6077276829592893E+00,  // c[9]
                                                                -6.7211819034913578E-02,  // c[10]
                                                                -1.2086015630786526E-02,  // c[11]
                                                                3.5480756226971832E-03,  // c[12]
                                                                4.2837323006755246E-03,  // c[13]
                                                                2.4295209811072277E-03,  // c[14]
                                                                1.0962945315780015E-03,  // c[15]
                                                                4.3635068246852815E-04,  // c[16]
                                                                1.5628005268397693E-04,  // c[17]
                                                                4.9233646028018496E-05,  // c[18]
                                                                1.2497401632990957E-05,  // c[19]
                                                                1.6428402654367069E-06,  // c[20]
                                                                -7.4570808553109380E-07,  // c[21]
                                                                -8.2894289515911268E-07,  // c[22]
                                                                -5.1089898456257998E-07,  // c[23]
                                                  };
alignas(FLOATING_PRECISION)  static constexpr std::size_t   ncx[24] = {
                                                                0,  // ncx[0]
                                                                1,  // ncx[1]
                                                                2,  // ncx[2]
                                                                3,  // ncx[3]
                                                                4,  // ncx[4]
                                                                5,  // ncx[5]
                                                                6,  // ncx[6]
                                                                7,  // ncx[7]
                                                                0,  // ncx[8]
                                                                1,  // ncx[9]
                                                                2,  // ncx[10]
                                                                3,  // ncx[11]
                                                                4,  // ncx[12]
                                                                5,  // ncx[13]
                                                                6,  // ncx[14]
                                                                7,  // ncx[15]
                                                                8,  // ncx[16]
                                                                9,  // ncx[17]
                                                                10,  // ncx[18]
                                                                11,  // ncx[19]
                                                                12,  // ncx[20]
                                                                13,  // ncx[21]
                                                                14,  // ncx[22]
                                                                15,  // ncx[23]
                                                  };
};
