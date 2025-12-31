
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

#include "../../../config.hpp"

struct M2C
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
                                                                2.6684151694226572E+00,  // c[0]
                                                                -1.4524115697730069E+00,  // c[1]
                                                                -6.5621388989946516E-03,  // c[2]
                                                                -2.2979513772384552E-03,  // c[3]
                                                                -5.9873203163895994E-04,  // c[4]
                                                                -1.1685992090226804E-04,  // c[5]
                                                                -1.5906097340728786E-05,  // c[6]
                                                                -8.7512448787574471E-07,  // c[7]
                                                                -1.2586190156114321E-01,  // c[8]
                                                                -1.1084468343880072E+00,  // c[9]
                                                                3.9515333331942448E-01,  // c[10]
                                                                2.3858195190833995E-01,  // c[11]
                                                                7.5832673358674829E-02,  // c[12]
                                                                1.8843711086912393E-03,  // c[13]
                                                                -1.2102147310964866E-02,  // c[14]
                                                                -7.6129670969728863E-03,  // c[15]
                                                                -2.7364065530609445E-03,  // c[16]
                                                                -5.4882009948363508E-04,  // c[17]
                                                                4.1979700859465230E-05,  // c[18]
                                                                9.3811245582921254E-05,  // c[19]
                                                                4.9182916275029853E-05,  // c[20]
                                                                1.7512321032656097E-05,  // c[21]
                                                                4.6681789573899124E-06,  // c[22]
                                                                8.7230136432137062E-07,  // c[23]
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
