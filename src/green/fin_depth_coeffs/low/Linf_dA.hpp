
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

struct Linf_dAC
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
                                                                -6.2959914138698368E-03,  // c[0]
                                                                -3.2573142608562620E-04,  // c[1]
                                                                -8.4308661720972742E-05,  // c[2]
                                                                -1.6620177515426932E-06,  // c[3]
                                                                -6.2411940913033332E-03,  // c[4]
                                                                -3.1920915290672947E-04,  // c[5]
                                                                -8.2573730358728422E-05,  // c[6]
                                                                -1.6012969920774429E-06,  // c[7]
                                                                6.5464693576109436E-05,  // c[8]
                                                                7.7660065412270122E-06,  // c[9]
                                                                2.0652099399090865E-06,  // c[10]
                                                                1.0531845343611203E-05,  // c[11]
                                                                1.2157075008860504E-06,  // c[12]
                                                                3.2257198640396178E-07,  // c[13]
                                                                -7.7350008869547550E-03,  // c[14]
                                                                -1.1516428669770470E-03,  // c[15]
                                                                -1.3077901767104586E-04,  // c[16]
                                                                -6.6498766521045163E-06,  // c[17]
                                                                -4.5380606152821109E-07,  // c[18]
                                                                -7.6489078555122079E-03,  // c[19]
                                                                -1.1253423997720450E-03,  // c[20]
                                                                -1.2718366784750213E-04,  // c[21]
                                                                -6.3637450926163556E-06,  // c[22]
                                                                -4.3018573631445215E-07,  // c[23]
                                                                1.0269529143377038E-04,  // c[24]
                                                                3.1272038167644766E-05,  // c[25]
                                                                4.2670471429907017E-06,  // c[26]
                                                                3.3825120241600119E-07,  // c[27]
                                                                1.6317402404165287E-05,  // c[28]
                                                                4.8388455277455938E-06,  // c[29]
                                                                6.4998174002393796E-07,  // c[30]
                                                                -3.1391305927894967E-07,  // c[31]
                                                                -1.7895330380916355E-02,  // c[32]
                                                                -8.6396862987175188E-04,  // c[33]
                                                                -2.2290649955104678E-04,  // c[34]
                                                                -3.9905342981945760E-06,  // c[35]
                                                                -5.2922909271653907E-07,  // c[36]
                                                                -5.2709310772034142E-03,  // c[37]
                                                                -2.1044397164662724E-04,  // c[38]
                                                                -5.3775631489094569E-05,  // c[39]
                                                                -6.6801108872867766E-07,  // c[40]
                                                                1.6975828815613041E-04,  // c[41]
                                                                1.8011671643874483E-05,  // c[42]
                                                                4.7487740715714379E-06,  // c[43]
                                                                6.3789457669828761E-06,  // c[44]
                                                                4.2018828395667604E-07,  // c[45]
                                                                -3.3674192814400761E-07,  // c[46]
                                                                -2.1674826335688638E-02,  // c[47]
                                                                -3.0059816433622710E-03,  // c[48]
                                                                -3.3233224422760110E-04,  // c[49]
                                                                -1.5399690058097692E-05,  // c[50]
                                                                -9.9849453880052844E-07,  // c[51]
                                                                -6.1645376640246420E-03,  // c[52]
                                                                -6.9705857868980557E-04,  // c[53]
                                                                -7.0480796261064614E-05,  // c[54]
                                                                -2.1544408983466745E-06,  // c[55]
                                                                2.5373182080468874E-04,  // c[56]
                                                                6.9386794247720276E-05,  // c[57]
                                                                8.9077062837754437E-06,  // c[58]
                                                                6.1485007877502265E-07,  // c[59]
                                                                8.0507270027235442E-06,  // c[60]
                                                                1.2398832380939835E-06,  // c[61]
                                                                -6.3322687772227612E-07,  // c[62]
                                                                -2.5468387513782192E-07,  // c[63]
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

