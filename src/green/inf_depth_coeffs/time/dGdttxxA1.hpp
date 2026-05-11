#ifndef __dGdttxxA1_coeffs_hpp
#define __dGdttxxA1_coeffs_hpp

#include "../../../config.hpp"

struct dGdttxxA1C
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
                                                                -1838.1498184782133194,  // c[0]
                                                                1.6083451655054830,  // c[1]
                                                                0.7890971399722275,  // c[2]
                                                                0.2761799265953186,  // c[3]
                                                                0.0741600161180713,  // c[4]
                                                                -1746.3191601333842300,  // c[5]
                                                                133.4795658521625512,  // c[6]
                                                                58.6107024337480880,  // c[7]
                                                                17.1786315579398092,  // c[8]
                                                                3.3176773810082523,  // c[9]
                                                                -963.6119341695094818,  // c[10]
                                                                442.8437791063237228,  // c[11]
                                                                -169.5074025167886020,  // c[12]
                                                                -4.7274785533811041,  // c[13]
                                                                37.9119789305141808,  // c[14]
                                                                -404.9363545020349306,  // c[15]
                                                                363.2171815903428183,  // c[16]
                                                                117.6308571589736971,  // c[17]
                                                                16.0418344177450081,  // c[18]
                                                                -4.5031978958005539,  // c[19]
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
