#ifndef __dGdttxA1_coeffs_hpp
#define __dGdttxA1_coeffs_hpp

#include "../../config.hpp"

struct dGdttxA1C
{
static constexpr int                                        max_ref_level = 2;
static constexpr int                                        intervals_np = 4;
static constexpr int                                        max_cheby_order = 4;
static constexpr int                                        max_cheby_order_f = 4;
static constexpr std::size_t                                blocks_start[4] = {0, 5, 10, 15};
static constexpr std::size_t                                blocks_start_fall[4] = {0, 0, 0, 0};
static constexpr std::size_t                                blocks_coeffs_np[4] = {5, 5, 5, 5};
static constexpr std::size_t                                blocks_coeffs_np_fall[4] = {1, 1, 1, 1};
static constexpr std::size_t                                blocks_max_cheby_order[4] = {4, 4, 4, 4};
static constexpr std::size_t                                blocks_max_cheby_order_fall[4] = {4, 4, 4, 4};
static constexpr cusfloat                                   x_min_region[4] = {-4.000E+00, -3.000E+00, -2.000E+00, -1.000E+00};
static constexpr cusfloat                                   x_max_region[4] = {-3.000E+00, -2.000E+00, -1.000E+00, -4.343E-05};
static constexpr cusfloat                                   dx_region[4] = {1.000E+00, 1.000E+00, 1.000E+00, 1.000E+00};
static constexpr bool                                       fcn_log_scale = 0;
static constexpr bool                                       x_log_scale = 0;
static constexpr cusfloat                                   x_max_global = -0.0000434316198075;
static constexpr cusfloat                                   x_min_global = -4.0000000000000000;
static constexpr cusfloat                                   dx_min_region = 0.9999891420950481;

static constexpr std::size_t                                num_c = 20;
alignas(FLOATING_PRECISION)  static constexpr cusfloat      c[20] = {
                                                                -755.4073355608804832,  // c[0]
                                                                17.7353964297403124,  // c[1]
                                                                8.0918854747870910,  // c[2]
                                                                2.5773794183444068,  // c[3]
                                                                0.6014030113704507,  // c[4]
                                                                -508.9819374026350829,  // c[5]
                                                                238.2524171555323278,  // c[6]
                                                                -4.4365700976640596,  // c[7]
                                                                -32.6902735114172600,  // c[8]
                                                                -4.3801003134288123,  // c[9]
                                                                -296.5301182037686090,  // c[10]
                                                                13.8390073962552975,  // c[11]
                                                                5.2278167020401156,  // c[12]
                                                                1.3975764164682971,  // c[13]
                                                                0.1835342344840285,  // c[14]
                                                                -237.5661859347375753,  // c[15]
                                                                -5.4873327098973803,  // c[16]
                                                                -91.0562255759949153,  // c[17]
                                                                -75.0975818274898472,  // c[18]
                                                                -40.4291219305139649,  // c[19]
                                                  };
alignas(FLOATING_PRECISION)  static constexpr std::size_t   ncx[20] = {
                                                                0,  // ncx[0]
                                                                1,  // ncx[1]
                                                                2,  // ncx[2]
                                                                3,  // ncx[3]
                                                                4,  // ncx[4]
                                                                0,  // ncx[5]
                                                                1,  // ncx[6]
                                                                2,  // ncx[7]
                                                                3,  // ncx[8]
                                                                4,  // ncx[9]
                                                                0,  // ncx[10]
                                                                1,  // ncx[11]
                                                                2,  // ncx[12]
                                                                3,  // ncx[13]
                                                                4,  // ncx[14]
                                                                0,  // ncx[15]
                                                                1,  // ncx[16]
                                                                2,  // ncx[17]
                                                                3,  // ncx[18]
                                                                4,  // ncx[19]
                                                  };
};
#endif
