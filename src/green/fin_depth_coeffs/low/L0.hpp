
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

struct L0C
{
static constexpr int                                        max_ref_level = 1;
static constexpr int                                        intervals_np = 2;
static constexpr int                                        max_cheby_order = 4;
static constexpr std::size_t                                blocks_start[4] = {0, 15, 31, 46};
static constexpr std::size_t                                blocks_coeffs_np[4] = {15, 16, 15, 17};
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
                                                                -6.9315715694827843E+04,  // c[0]
                                                                6.2961094749880431E-03,  // c[1]
                                                                1.5894788666628301E-03,  // c[2]
                                                                8.8769697867974173E-06,  // c[3]
                                                                1.1461829672043677E-06,  // c[4]
                                                                -3.2150567494682036E-03,  // c[5]
                                                                -1.0900347479037009E-04,  // c[6]
                                                                -2.8079316052753711E-05,  // c[7]
                                                                -4.7850608098087832E-07,  // c[8]
                                                                -7.9745902075956110E-04,  // c[9]
                                                                -2.6604726144796587E-05,  // c[10]
                                                                -6.8485755946312565E-06,  // c[11]
                                                                3.5850966924044769E-06,  // c[12]
                                                                3.6573760553437751E-07,  // c[13]
                                                                4.3450688735902077E-07,  // c[14]
                                                                -6.9315689811791613E+04,  // c[15]
                                                                1.9776952704660289E-02,  // c[16]
                                                                1.8185111307502666E-03,  // c[17]
                                                                3.0603461027567391E-05,  // c[18]
                                                                1.6761546248744708E-06,  // c[19]
                                                                -3.6900993552535510E-03,  // c[20]
                                                                -3.7700516804761719E-04,  // c[21]
                                                                -4.1356244082635385E-05,  // c[22]
                                                                -1.8884688870457467E-06,  // c[23]
                                                                -9.1315182066864509E-04,  // c[24]
                                                                -9.1688200200223946E-05,  // c[25]
                                                                -9.9955839232279686E-06,  // c[26]
                                                                -4.4573471313924529E-07,  // c[27]
                                                                5.3195838063402334E-06,  // c[28]
                                                                1.4481656762654893E-06,  // c[29]
                                                                6.3720608522999100E-07,  // c[30]
                                                                -6.9315728290716899E+04,  // c[31]
                                                                5.8863408769411762E-03,  // c[32]
                                                                1.4841016487707748E-03,  // c[33]
                                                                7.1830405659056851E-06,  // c[34]
                                                                9.2153561581653776E-07,  // c[35]
                                                                -9.3132996116764843E-03,  // c[36]
                                                                -2.9424198010019609E-04,  // c[37]
                                                                -7.5576356721285265E-05,  // c[38]
                                                                -1.1627372487055254E-06,  // c[39]
                                                                -7.1670568536319479E-04,  // c[40]
                                                                -1.8828821794159012E-05,  // c[41]
                                                                -4.7968903800210683E-06,  // c[42]
                                                                9.5027016868698411E-06,  // c[43]
                                                                8.6307954916264862E-07,  // c[44]
                                                                2.8763440695911413E-07,  // c[45]
                                                                -6.9315704184163391E+04,  // c[46]
                                                                1.8371819652656995E-02,  // c[47]
                                                                1.6666528399582603E-03,  // c[48]
                                                                2.4055085305008106E-05,  // c[49]
                                                                1.2495073633544962E-06,  // c[50]
                                                                -1.0584196034187698E-02,  // c[51]
                                                                -1.0028177803178551E-03,  // c[52]
                                                                -1.0721240619204764E-04,  // c[53]
                                                                -4.4240733814149280E-06,  // c[54]
                                                                -2.8152101094747195E-07,  // c[55]
                                                                -7.9601295010434114E-04,  // c[56]
                                                                -6.1549070778710302E-05,  // c[57]
                                                                -6.0841116464871448E-06,  // c[58]
                                                                1.3486813486451865E-05,  // c[59]
                                                                3.2733930765971309E-06,  // c[60]
                                                                4.1171097109327093E-07,  // c[61]
                                                                3.5456059777061455E-07,  // c[62]
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
                                                                4,  // ncx[45]
                                                                0,  // ncx[46]
                                                                0,  // ncx[47]
                                                                0,  // ncx[48]
                                                                0,  // ncx[49]
                                                                0,  // ncx[50]
                                                                1,  // ncx[51]
                                                                1,  // ncx[52]
                                                                1,  // ncx[53]
                                                                1,  // ncx[54]
                                                                1,  // ncx[55]
                                                                2,  // ncx[56]
                                                                2,  // ncx[57]
                                                                2,  // ncx[58]
                                                                3,  // ncx[59]
                                                                3,  // ncx[60]
                                                                3,  // ncx[61]
                                                                4,  // ncx[62]
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
                                                                0,  // ncy[46]
                                                                1,  // ncy[47]
                                                                2,  // ncy[48]
                                                                3,  // ncy[49]
                                                                4,  // ncy[50]
                                                                0,  // ncy[51]
                                                                1,  // ncy[52]
                                                                2,  // ncy[53]
                                                                3,  // ncy[54]
                                                                4,  // ncy[55]
                                                                0,  // ncy[56]
                                                                1,  // ncy[57]
                                                                2,  // ncy[58]
                                                                0,  // ncy[59]
                                                                1,  // ncy[60]
                                                                2,  // ncy[61]
                                                                0,  // ncy[62]
                                                  };
};

