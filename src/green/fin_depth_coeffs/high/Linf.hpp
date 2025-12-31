
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

struct LinfC
{
static constexpr int                                        max_ref_level = 1;
static constexpr int                                        intervals_np = 2;
static constexpr int                                        max_cheby_order = 4;
static constexpr std::size_t                                blocks_start[4] = {0, 15, 31, 45};
static constexpr std::size_t                                blocks_coeffs_np[4] = {15, 16, 14, 17};
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

static constexpr std::size_t                                num_c = 62;
alignas(FLOATING_PRECISION)  static constexpr cusfloat      c[62] = {
                                                                3.0800155877736268E-01,  // c[0]
                                                                3.0630567800966871E-03,  // c[1]
                                                                7.7738732407384636E-04,  // c[2]
                                                                6.6838018989185278E-06,  // c[3]
                                                                8.6778875794122148E-07,  // c[4]
                                                                -1.5821809401644167E-03,  // c[5]
                                                                -8.2403607339100876E-05,  // c[6]
                                                                -2.1335316672838056E-05,  // c[7]
                                                                -4.2450056366707609E-07,  // c[8]
                                                                -3.9073287104041295E-04,  // c[9]
                                                                -2.0026553775441547E-05,  // c[10]
                                                                -5.1810188966081722E-06,  // c[11]
                                                                2.7339326412485071E-06,  // c[12]
                                                                3.2486943656907796E-07,  // c[13]
                                                                3.2955509013388173E-07,  // c[14]
                                                                3.2078299938234400E-01,  // c[15]
                                                                9.8631429787881548E-03,  // c[16]
                                                                9.5199151214065292E-04,  // c[17]
                                                                2.3584550622076034E-05,  // c[18]
                                                                1.3411571864723609E-06,  // c[19]
                                                                -1.9465871331678557E-03,  // c[20]
                                                                -2.9181972151523169E-04,  // c[21]
                                                                -3.3228135310712018E-05,  // c[22]
                                                                -1.7047505633522968E-06,  // c[23]
                                                                -4.7907657861975990E-04,  // c[24]
                                                                -7.0636327831167384E-05,  // c[25]
                                                                -7.9896030992371303E-06,  // c[26]
                                                                -4.0084606962127386E-07,  // c[27]
                                                                4.2920501871568720E-06,  // c[28]
                                                                1.3090792409936602E-06,  // c[29]
                                                                5.1080461809915184E-07,  // c[30]
                                                                3.0187336340288523E-01,  // c[31]
                                                                2.7566624726674377E-03,  // c[32]
                                                                6.9820945988715915E-04,  // c[33]
                                                                5.1941472592827078E-06,  // c[34]
                                                                6.6942272511230985E-07,  // c[35]
                                                                -4.4950523812485911E-03,  // c[36]
                                                                -2.1824361642347467E-04,  // c[37]
                                                                -5.6320221646824487E-05,  // c[38]
                                                                -1.0154535167242185E-06,  // c[39]
                                                                -3.2983187643566367E-04,  // c[40]
                                                                -1.3179009995640048E-05,  // c[41]
                                                                -3.3675844457225004E-06,  // c[42]
                                                                7.0872929201467108E-06,  // c[43]
                                                                7.5295614350955054E-07,  // c[44]
                                                                3.1330782489146991E-01,  // c[45]
                                                                8.7884752322347492E-03,  // c[46]
                                                                8.3155541009932851E-04,  // c[47]
                                                                1.7728025581955666E-05,  // c[48]
                                                                9.5114117975789170E-07,  // c[49]
                                                                -5.4504230615226076E-03,  // c[50]
                                                                -7.6016876012140536E-04,  // c[51]
                                                                -8.4196524342373263E-05,  // c[52]
                                                                -3.9267787743627221E-06,  // c[53]
                                                                -2.5533791360855644E-07,  // c[54]
                                                                -3.8578677443921756E-04,  // c[55]
                                                                -4.3643653870439891E-05,  // c[56]
                                                                -4.4106614996508933E-06,  // c[57]
                                                                1.0598543653398892E-05,  // c[58]
                                                                2.9017282550895113E-06,  // c[59]
                                                                3.7270829793614735E-07,  // c[60]
                                                                2.5165257904347788E-07,  // c[61]
                                                  };
alignas(FLOATING_PRECISION)  static constexpr std::size_t   ncx[62] = {
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
                                                                3,  // ncx[12]
                                                                3,  // ncx[13]
                                                                4,  // ncx[14]
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
                                                                4,  // ncx[30]
                                                                0,  // ncx[31]
                                                                0,  // ncx[32]
                                                                0,  // ncx[33]
                                                                0,  // ncx[34]
                                                                0,  // ncx[35]
                                                                1,  // ncx[36]
                                                                1,  // ncx[37]
                                                                1,  // ncx[38]
                                                                1,  // ncx[39]
                                                                2,  // ncx[40]
                                                                2,  // ncx[41]
                                                                2,  // ncx[42]
                                                                3,  // ncx[43]
                                                                3,  // ncx[44]
                                                                0,  // ncx[45]
                                                                0,  // ncx[46]
                                                                0,  // ncx[47]
                                                                0,  // ncx[48]
                                                                0,  // ncx[49]
                                                                1,  // ncx[50]
                                                                1,  // ncx[51]
                                                                1,  // ncx[52]
                                                                1,  // ncx[53]
                                                                1,  // ncx[54]
                                                                2,  // ncx[55]
                                                                2,  // ncx[56]
                                                                2,  // ncx[57]
                                                                3,  // ncx[58]
                                                                3,  // ncx[59]
                                                                3,  // ncx[60]
                                                                4,  // ncx[61]
                                                  };
alignas(FLOATING_PRECISION)  static constexpr std::size_t   ncy[62] = {
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
                                                                0,  // ncy[12]
                                                                1,  // ncy[13]
                                                                0,  // ncy[14]
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
                                                                0,  // ncy[30]
                                                                0,  // ncy[31]
                                                                1,  // ncy[32]
                                                                2,  // ncy[33]
                                                                3,  // ncy[34]
                                                                4,  // ncy[35]
                                                                0,  // ncy[36]
                                                                1,  // ncy[37]
                                                                2,  // ncy[38]
                                                                3,  // ncy[39]
                                                                0,  // ncy[40]
                                                                1,  // ncy[41]
                                                                2,  // ncy[42]
                                                                0,  // ncy[43]
                                                                1,  // ncy[44]
                                                                0,  // ncy[45]
                                                                1,  // ncy[46]
                                                                2,  // ncy[47]
                                                                3,  // ncy[48]
                                                                4,  // ncy[49]
                                                                0,  // ncy[50]
                                                                1,  // ncy[51]
                                                                2,  // ncy[52]
                                                                3,  // ncy[53]
                                                                4,  // ncy[54]
                                                                0,  // ncy[55]
                                                                1,  // ncy[56]
                                                                2,  // ncy[57]
                                                                0,  // ncy[58]
                                                                1,  // ncy[59]
                                                                2,  // ncy[60]
                                                                0,  // ncy[61]
                                                  };
};

