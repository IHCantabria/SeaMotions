
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

struct L0_dAC
{
static constexpr int                                        max_ref_level = 1;
static constexpr int                                        intervals_np = 2;
static constexpr int                                        max_cheby_order = 4;
static constexpr std::size_t                                blocks_start[4] = {0, 14, 32, 47};
static constexpr std::size_t                                blocks_coeffs_np[4] = {14, 18, 15, 17};
static constexpr std::size_t                                blocks_max_cheby_order[4] = {3, 4, 4, 4};
static constexpr cusfloat                                   x_min_region[4] = {0.000000000000000E+00, 0.000000000000000E+00, 5.000000000000000E-01, 5.000000000000000E-01};
static constexpr cusfloat                                   x_max_region[4] = {5.000000000000000E-01, 5.000000000000000E-01, 1.000000000000000E+00, 1.000000000000000E+00};
static constexpr cusfloat                                   dx_region[4] = {5.000000000000000E-01, 5.000000000000000E-01, 5.000000000000000E-01, 5.000000000000000E-01};
static constexpr cusfloat                                   y_min_region[4] = {0.000000000000000E+00, 5.000000000000000E-01, 0.000000000000000E+00, 5.000000000000000E-01};
static constexpr cusfloat                                   y_max_region[4] = {5.000000000000000E-01, 1.000000000000000E+00, 5.000000000000000E-01, 1.000000000000000E+00};
static constexpr cusfloat                                   dy_region[4] = {5.000000000000000E-01, 5.000000000000000E-01, 5.000000000000000E-01, 5.000000000000000E-01};
static constexpr bool                                       fcn_log_scale = 0;
static constexpr bool                                       x_log_scale = 0;
static constexpr cusfloat                                   x_max_global = 1.0000000000000000;
static constexpr cusfloat                                   x_min_global = 0.0000000000000000;
static constexpr cusfloat                                   dx_min_region = 0.5000000000000000;
static constexpr bool                                       y_log_scale = 0;
static constexpr cusfloat                                   y_max_global = 1.0000000000000000;
static constexpr cusfloat                                   y_min_global = 0.0000000000000000;
static constexpr cusfloat                                   dy_min_region = 0.5000000000000000;

static constexpr std::size_t                                num_c = 64;
alignas(FLOATING_PRECISION)  static constexpr cusfloat      c[64] = {
                                                                -1.2817289658354530E-02,  // c[0]
                                                                -4.3164145314454757E-04,  // c[1]
                                                                -1.1115903309811138E-04,  // c[2]
                                                                -1.8761297723288086E-06,  // c[3]
                                                                -1.2745455389584111E-02,  // c[4]
                                                                -4.2430074882661407E-04,  // c[5]
                                                                -1.0921386954131725E-04,  // c[6]
                                                                -1.8121692341893138E-06,  // c[7]
                                                                8.5874673188206349E-05,  // c[8]
                                                                8.7449582309583225E-06,  // c[9]
                                                                2.3166618298167902E-06,  // c[10]
                                                                1.3888840602223720E-05,  // c[11]
                                                                1.3747538467381295E-06,  // c[12]
                                                                3.6340673103467218E-07,  // c[13]
                                                                -1.4696732102782694E-02,  // c[14]
                                                                -1.4907187608924561E-03,  // c[15]
                                                                -1.6311420816020371E-04,  // c[16]
                                                                -7.3776840160155212E-06,  // c[17]
                                                                -4.9258922775848512E-07,  // c[18]
                                                                -1.4590069087749478E-02,  // c[19]
                                                                -1.4616320609882483E-03,  // c[20]
                                                                -1.5922212245275159E-04,  // c[21]
                                                                -7.0796533581470489E-06,  // c[22]
                                                                -4.6825145360706449E-07,  // c[23]
                                                                1.2733057707856033E-04,  // c[24]
                                                                3.4603870885589034E-05,  // c[25]
                                                                4.6217300090161360E-06,  // c[26]
                                                                3.5244788430729155E-07,  // c[27]
                                                                2.0359991600144402E-05,  // c[28]
                                                                5.3790254099582262E-06,  // c[29]
                                                                7.0727771138065759E-07,  // c[30]
                                                                -3.3907692866962276E-07,  // c[31]
                                                                -3.7139358272961775E-02,  // c[32]
                                                                -1.1666427302315339E-03,  // c[33]
                                                                -2.9959284202309838E-04,  // c[34]
                                                                -4.5748356539235550E-06,  // c[35]
                                                                -6.0433352539769874E-07,  // c[36]
                                                                -1.1458091575541254E-02,  // c[37]
                                                                -3.0072934536541336E-04,  // c[38]
                                                                -7.6616156118578622E-05,  // c[39]
                                                                -8.2241932233178064E-07,  // c[40]
                                                                2.2768029268995210E-04,  // c[41]
                                                                2.0650439024837413E-05,  // c[42]
                                                                5.4253726435456694E-06,  // c[43]
                                                                9.1993163372256111E-06,  // c[44]
                                                                5.3165882617852052E-07,  // c[45]
                                                                -3.8424414237881420E-07,  // c[46]
                                                                -4.2175291978250054E-02,  // c[47]
                                                                -3.9721252378074621E-03,  // c[48]
                                                                -4.2392883750144987E-04,  // c[49]
                                                                -1.7371333724297891E-05,  // c[50]
                                                                -1.1022011849899497E-06,  // c[51]
                                                                -1.2724864562036640E-02,  // c[52]
                                                                -9.8317826511477900E-04,  // c[53]
                                                                -9.7219321134630774E-05,  // c[54]
                                                                -2.6648403394001572E-06,  // c[55]
                                                                3.2298434034904583E-04,  // c[56]
                                                                7.8291851063562370E-05,  // c[57]
                                                                9.8417709505367167E-06,  // c[58]
                                                                6.4998273230148715E-07,  // c[59]
                                                                1.1342513859657313E-05,  // c[60]
                                                                1.6067455919331619E-06,  // c[61]
                                                                -6.9891647765584972E-07,  // c[62]
                                                                -2.6935044052296166E-07,  // c[63]
                                                  };
alignas(FLOATING_PRECISION)  static constexpr std::size_t   ncx[64] = {
                                                                0,  // ncx[0]
                                                                0,  // ncx[1]
                                                                0,  // ncx[2]
                                                                0,  // ncx[3]
                                                                1,  // ncx[4]
                                                                1,  // ncx[5]
                                                                1,  // ncx[6]
                                                                1,  // ncx[7]
                                                                2,  // ncx[8]
                                                                2,  // ncx[9]
                                                                2,  // ncx[10]
                                                                3,  // ncx[11]
                                                                3,  // ncx[12]
                                                                3,  // ncx[13]
                                                                0,  // ncx[14]
                                                                0,  // ncx[15]
                                                                0,  // ncx[16]
                                                                0,  // ncx[17]
                                                                0,  // ncx[18]
                                                                1,  // ncx[19]
                                                                1,  // ncx[20]
                                                                1,  // ncx[21]
                                                                1,  // ncx[22]
                                                                1,  // ncx[23]
                                                                2,  // ncx[24]
                                                                2,  // ncx[25]
                                                                2,  // ncx[26]
                                                                2,  // ncx[27]
                                                                3,  // ncx[28]
                                                                3,  // ncx[29]
                                                                3,  // ncx[30]
                                                                4,  // ncx[31]
                                                                0,  // ncx[32]
                                                                0,  // ncx[33]
                                                                0,  // ncx[34]
                                                                0,  // ncx[35]
                                                                0,  // ncx[36]
                                                                1,  // ncx[37]
                                                                1,  // ncx[38]
                                                                1,  // ncx[39]
                                                                1,  // ncx[40]
                                                                2,  // ncx[41]
                                                                2,  // ncx[42]
                                                                2,  // ncx[43]
                                                                3,  // ncx[44]
                                                                3,  // ncx[45]
                                                                4,  // ncx[46]
                                                                0,  // ncx[47]
                                                                0,  // ncx[48]
                                                                0,  // ncx[49]
                                                                0,  // ncx[50]
                                                                0,  // ncx[51]
                                                                1,  // ncx[52]
                                                                1,  // ncx[53]
                                                                1,  // ncx[54]
                                                                1,  // ncx[55]
                                                                2,  // ncx[56]
                                                                2,  // ncx[57]
                                                                2,  // ncx[58]
                                                                2,  // ncx[59]
                                                                3,  // ncx[60]
                                                                3,  // ncx[61]
                                                                4,  // ncx[62]
                                                                4,  // ncx[63]
                                                  };
alignas(FLOATING_PRECISION)  static constexpr std::size_t   ncy[64] = {
                                                                0,  // ncy[0]
                                                                1,  // ncy[1]
                                                                2,  // ncy[2]
                                                                3,  // ncy[3]
                                                                0,  // ncy[4]
                                                                1,  // ncy[5]
                                                                2,  // ncy[6]
                                                                3,  // ncy[7]
                                                                0,  // ncy[8]
                                                                1,  // ncy[9]
                                                                2,  // ncy[10]
                                                                0,  // ncy[11]
                                                                1,  // ncy[12]
                                                                2,  // ncy[13]
                                                                0,  // ncy[14]
                                                                1,  // ncy[15]
                                                                2,  // ncy[16]
                                                                3,  // ncy[17]
                                                                4,  // ncy[18]
                                                                0,  // ncy[19]
                                                                1,  // ncy[20]
                                                                2,  // ncy[21]
                                                                3,  // ncy[22]
                                                                4,  // ncy[23]
                                                                0,  // ncy[24]
                                                                1,  // ncy[25]
                                                                2,  // ncy[26]
                                                                3,  // ncy[27]
                                                                0,  // ncy[28]
                                                                1,  // ncy[29]
                                                                2,  // ncy[30]
                                                                0,  // ncy[31]
                                                                0,  // ncy[32]
                                                                1,  // ncy[33]
                                                                2,  // ncy[34]
                                                                3,  // ncy[35]
                                                                4,  // ncy[36]
                                                                0,  // ncy[37]
                                                                1,  // ncy[38]
                                                                2,  // ncy[39]
                                                                3,  // ncy[40]
                                                                0,  // ncy[41]
                                                                1,  // ncy[42]
                                                                2,  // ncy[43]
                                                                0,  // ncy[44]
                                                                1,  // ncy[45]
                                                                0,  // ncy[46]
                                                                0,  // ncy[47]
                                                                1,  // ncy[48]
                                                                2,  // ncy[49]
                                                                3,  // ncy[50]
                                                                4,  // ncy[51]
                                                                0,  // ncy[52]
                                                                1,  // ncy[53]
                                                                2,  // ncy[54]
                                                                3,  // ncy[55]
                                                                0,  // ncy[56]
                                                                1,  // ncy[57]
                                                                2,  // ncy[58]
                                                                3,  // ncy[59]
                                                                0,  // ncy[60]
                                                                1,  // ncy[61]
                                                                0,  // ncy[62]
                                                                1,  // ncy[63]
                                                  };
};

