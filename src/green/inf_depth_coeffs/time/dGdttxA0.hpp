#ifndef __dGdttxA0_coeffs_hpp
#define __dGdttxA0_coeffs_hpp

#include "../../config.hpp"

struct dGdttxA0C
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
                                                                688.5202510178934290,  // c[0]
                                                                -16.1455776336903867,  // c[1]
                                                                -7.3655277303765843,  // c[2]
                                                                -2.3458225683534408,  // c[3]
                                                                -0.5473425046283751,  // c[4]
                                                                464.2307539431501482,  // c[5]
                                                                -216.8455856371702737,  // c[6]
                                                                4.0373256815533409,  // c[7]
                                                                29.7506674751392381,  // c[8]
                                                                3.9856597608514335,  // c[9]
                                                                270.7868860249420777,  // c[10]
                                                                -12.6936811270449770,  // c[11]
                                                                -4.7822156042090569,  // c[12]
                                                                -1.2753131705873813,  // c[13]
                                                                -0.1671407319647216,  // c[14]
                                                                216.2746422016124654,  // c[15]
                                                                3.8901138078391977,  // c[16]
                                                                82.2313700146786459,  // c[17]
                                                                67.9376233877178777,  // c[18]
                                                                36.5771599760560306,  // c[19]
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
