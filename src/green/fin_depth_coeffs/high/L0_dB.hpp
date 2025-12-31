
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

struct L0_dBC
{
static constexpr int                                        max_ref_level = 1;
static constexpr int                                        intervals_np = 2;
static constexpr int                                        max_cheby_order = 4;
static constexpr std::size_t                                blocks_start[4] = {0, 15, 33, 47};
static constexpr std::size_t                                blocks_coeffs_np[4] = {15, 18, 14, 16};
static constexpr std::size_t                                blocks_max_cheby_order[4] = {4, 4, 4, 4};
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

static constexpr std::size_t                                num_c = 63;
alignas(FLOATING_PRECISION)  static constexpr cusfloat      c[63] = {
                                                                2.5291192493034000E-02,  // c[0]
                                                                2.5468388543865975E-02,  // c[1]
                                                                2.1350916345448515E-04,  // c[2]
                                                                3.6726662582643401E-05,  // c[3]
                                                                4.6204288250262628E-07,  // c[4]
                                                                -4.4178116159213185E-04,  // c[5]
                                                                -4.5131490977328402E-04,  // c[6]
                                                                -1.1534458653982959E-05,  // c[7]
                                                                -2.0454537057269036E-06,  // c[8]
                                                                -1.0779237389784042E-04,  // c[9]
                                                                -1.1006347948499433E-04,  // c[10]
                                                                -2.7469964593528334E-06,  // c[11]
                                                                -4.8620320617558197E-07,  // c[12]
                                                                1.5013313674538736E-06,  // c[13]
                                                                1.5645227569383371E-06,  // c[14]
                                                                7.9475958554285117E-02,  // c[15]
                                                                2.9149911238277669E-02,  // c[16]
                                                                7.3629543035319111E-04,  // c[17]
                                                                5.3733147265341071E-05,  // c[18]
                                                                1.8125181569653755E-06,  // c[19]
                                                                -1.5307988306283457E-03,  // c[20]
                                                                -6.6577096477910772E-04,  // c[21]
                                                                -4.5556268586235266E-05,  // c[22]
                                                                -4.0706721202726391E-06,  // c[23]
                                                                -3.7212787968457853E-04,  // c[24]
                                                                -1.6087751825505153E-04,  // c[25]
                                                                -1.0750216170038159E-05,  // c[26]
                                                                -9.4811809093261340E-07,  // c[27]
                                                                5.9716780256087627E-06,  // c[28]
                                                                3.1363480870161796E-06,  // c[29]
                                                                3.5811216075659277E-07,  // c[30]
                                                                6.9379277977321338E-07,  // c[31]
                                                                3.5885230096916604E-07,  // c[32]
                                                                2.3631708817943169E-02,  // c[33]
                                                                2.3775146347769845E-02,  // c[34]
                                                                1.7269057987098920E-04,  // c[35]
                                                                2.9519965119163068E-05,  // c[36]
                                                                2.9775102871254164E-07,  // c[37]
                                                                -1.1909738093214182E-03,  // c[38]
                                                                -1.2141510806542586E-03,  // c[39]
                                                                -2.8011719155069791E-05,  // c[40]
                                                                -4.9289673067753870E-06,  // c[41]
                                                                -7.5933438579959926E-05,  // c[42]
                                                                -7.6960862807686674E-05,  // c[43]
                                                                -1.2363750640025876E-06,  // c[44]
                                                                3.5292416403078604E-06,  // c[45]
                                                                3.6561047099484961E-06,  // c[46]
                                                                7.3776479193214381E-02,  // c[47]
                                                                2.6706480204476506E-02,  // c[48]
                                                                5.7840112105234543E-04,  // c[49]
                                                                4.0034751002359465E-05,  // c[50]
                                                                1.0792512145504311E-06,  // c[51]
                                                                -4.0645885812131497E-03,  // c[52]
                                                                -1.7244349203966646E-03,  // c[53]
                                                                -1.0663483574784504E-04,  // c[54]
                                                                -9.0360284898000640E-06,  // c[55]
                                                                -4.5688520505093996E-07,  // c[56]
                                                                -2.4819227234885804E-04,  // c[57]
                                                                -9.7594358098614890E-05,  // c[58]
                                                                -3.9920393443443770E-06,  // c[59]
                                                                1.3422337935407869E-05,  // c[60]
                                                                6.6515360315957256E-06,  // c[61]
                                                                6.5760870278873916E-07,  // c[62]
                                                  };
alignas(FLOATING_PRECISION)  static constexpr std::size_t   ncx[63] = {
                                                                0,  // ncx[0]
                                                                0,  // ncx[1]
                                                                0,  // ncx[2]
                                                                0,  // ncx[3]
                                                                0,  // ncx[4]
                                                                1,  // ncx[5]
                                                                1,  // ncx[6]
                                                                1,  // ncx[7]
                                                                1,  // ncx[8]
                                                                2,  // ncx[9]
                                                                2,  // ncx[10]
                                                                2,  // ncx[11]
                                                                2,  // ncx[12]
                                                                3,  // ncx[13]
                                                                3,  // ncx[14]
                                                                0,  // ncx[15]
                                                                0,  // ncx[16]
                                                                0,  // ncx[17]
                                                                0,  // ncx[18]
                                                                0,  // ncx[19]
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
                                                                4,  // ncx[32]
                                                                0,  // ncx[33]
                                                                0,  // ncx[34]
                                                                0,  // ncx[35]
                                                                0,  // ncx[36]
                                                                0,  // ncx[37]
                                                                1,  // ncx[38]
                                                                1,  // ncx[39]
                                                                1,  // ncx[40]
                                                                1,  // ncx[41]
                                                                2,  // ncx[42]
                                                                2,  // ncx[43]
                                                                2,  // ncx[44]
                                                                3,  // ncx[45]
                                                                3,  // ncx[46]
                                                                0,  // ncx[47]
                                                                0,  // ncx[48]
                                                                0,  // ncx[49]
                                                                0,  // ncx[50]
                                                                0,  // ncx[51]
                                                                1,  // ncx[52]
                                                                1,  // ncx[53]
                                                                1,  // ncx[54]
                                                                1,  // ncx[55]
                                                                1,  // ncx[56]
                                                                2,  // ncx[57]
                                                                2,  // ncx[58]
                                                                2,  // ncx[59]
                                                                3,  // ncx[60]
                                                                3,  // ncx[61]
                                                                3,  // ncx[62]
                                                  };
alignas(FLOATING_PRECISION)  static constexpr std::size_t   ncy[63] = {
                                                                0,  // ncy[0]
                                                                1,  // ncy[1]
                                                                2,  // ncy[2]
                                                                3,  // ncy[3]
                                                                4,  // ncy[4]
                                                                0,  // ncy[5]
                                                                1,  // ncy[6]
                                                                2,  // ncy[7]
                                                                3,  // ncy[8]
                                                                0,  // ncy[9]
                                                                1,  // ncy[10]
                                                                2,  // ncy[11]
                                                                3,  // ncy[12]
                                                                0,  // ncy[13]
                                                                1,  // ncy[14]
                                                                0,  // ncy[15]
                                                                1,  // ncy[16]
                                                                2,  // ncy[17]
                                                                3,  // ncy[18]
                                                                4,  // ncy[19]
                                                                0,  // ncy[20]
                                                                1,  // ncy[21]
                                                                2,  // ncy[22]
                                                                3,  // ncy[23]
                                                                0,  // ncy[24]
                                                                1,  // ncy[25]
                                                                2,  // ncy[26]
                                                                3,  // ncy[27]
                                                                0,  // ncy[28]
                                                                1,  // ncy[29]
                                                                2,  // ncy[30]
                                                                0,  // ncy[31]
                                                                1,  // ncy[32]
                                                                0,  // ncy[33]
                                                                1,  // ncy[34]
                                                                2,  // ncy[35]
                                                                3,  // ncy[36]
                                                                4,  // ncy[37]
                                                                0,  // ncy[38]
                                                                1,  // ncy[39]
                                                                2,  // ncy[40]
                                                                3,  // ncy[41]
                                                                0,  // ncy[42]
                                                                1,  // ncy[43]
                                                                2,  // ncy[44]
                                                                0,  // ncy[45]
                                                                1,  // ncy[46]
                                                                0,  // ncy[47]
                                                                1,  // ncy[48]
                                                                2,  // ncy[49]
                                                                3,  // ncy[50]
                                                                4,  // ncy[51]
                                                                0,  // ncy[52]
                                                                1,  // ncy[53]
                                                                2,  // ncy[54]
                                                                3,  // ncy[55]
                                                                4,  // ncy[56]
                                                                0,  // ncy[57]
                                                                1,  // ncy[58]
                                                                2,  // ncy[59]
                                                                0,  // ncy[60]
                                                                1,  // ncy[61]
                                                                2,  // ncy[62]
                                                  };
};

