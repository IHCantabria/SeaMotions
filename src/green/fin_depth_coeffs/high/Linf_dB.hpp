
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

struct Linf_dBC
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
                                                                1.2332637187249185E-02,  // c[0]
                                                                1.2466010169700167E-02,  // c[1]
                                                                1.6082013372506314E-04,  // c[2]
                                                                2.7812984518737962E-05,  // c[3]
                                                                4.0888815004783705E-07,  // c[4]
                                                                -3.3473223706630942E-04,  // c[5]
                                                                -3.4318767995041868E-04,  // c[6]
                                                                -1.0235615420151186E-05,  // c[7]
                                                                -1.8226131867180462E-06,  // c[8]
                                                                -8.1320785543886138E-05,  // c[9]
                                                                -8.3327992637186759E-05,  // c[10]
                                                                -2.4291408839558415E-06,  // c[11]
                                                                -4.3169029206001353E-07,  // c[12]
                                                                1.3359556297438892E-06,  // c[13]
                                                                1.3959248036478963E-06,  // c[14]
                                                                3.9736402488564054E-02,  // c[15]
                                                                1.5274870127756465E-02,  // c[16]
                                                                5.6766114682868995E-04,  // c[17]
                                                                4.3005933514836815E-05,  // c[18]
                                                                1.6319319042542701E-06,  // c[19]
                                                                -1.1878475891177521E-03,  // c[20]
                                                                -5.3540669595105134E-04,  // c[21]
                                                                -4.1137406113579298E-05,  // c[22]
                                                                -3.7565309806165135E-06,  // c[23]
                                                                -2.8738065965874262E-04,  // c[24]
                                                                -1.2870525696419841E-04,  // c[25]
                                                                -9.6706966676239061E-06,  // c[26]
                                                                -8.7160737669109579E-07,  // c[27]
                                                                5.4082284714163055E-06,  // c[28]
                                                                2.8982658970017755E-06,  // c[29]
                                                                3.4382301470617035E-07,  // c[30]
                                                                6.2534735635600119E-07,  // c[31]
                                                                3.3003449786245774E-07,  // c[32]
                                                                1.1089106646735597E-02,  // c[33]
                                                                1.1192799562504142E-02,  // c[34]
                                                                1.2491351213200330E-04,  // c[35]
                                                                2.1448204309640361E-05,  // c[36]
                                                                2.5397790819829246E-07,  // c[37]
                                                                -8.8520952752959420E-04,  // c[38]
                                                                -9.0544590403802006E-04,  // c[39]
                                                                -2.4470123671817647E-05,  // c[40]
                                                                -4.3223576906965507E-06,  // c[41]
                                                                -5.3217998746977254E-05,  // c[42]
                                                                -5.4052274328949927E-05,  // c[43]
                                                                -1.0039175286583198E-06,  // c[44]
                                                                3.0838775193823783E-06,  // c[45]
                                                                3.2026321564338628E-06,  // c[46]
                                                                3.5367103638361694E-02,  // c[47]
                                                                1.3335367932399842E-02,  // c[48]
                                                                4.2640541885303980E-04,  // c[49]
                                                                3.0481370822053573E-05,  // c[50]
                                                                9.3280489347477092E-07,  // c[51]
                                                                -3.0880129112609948E-03,  // c[52]
                                                                -1.3553413929229300E-03,  // c[53]
                                                                -9.4675741549585745E-05,  // c[54]
                                                                -8.1970034446736265E-06,  // c[55]
                                                                -4.3305096417712023E-07,  // c[56]
                                                                -1.7618639719238049E-04,  // c[57]
                                                                -7.0767205735561408E-05,  // c[58]
                                                                -3.2235634207970027E-06,  // c[59]
                                                                1.1918048329286621E-05,  // c[60]
                                                                6.0249205855174119E-06,  // c[61]
                                                                6.2227061762350645E-07,  // c[62]
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

