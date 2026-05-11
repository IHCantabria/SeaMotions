#ifndef __dGdttxxA2_coeffs_hpp
#define __dGdttxxA2_coeffs_hpp

#include "../../../config.hpp"

struct dGdttxxA2C
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
                                                                168.1881937999417289,  // c[0]
                                                                -0.1398988315737473,  // c[1]
                                                                -0.0702917808814902,  // c[2]
                                                                -0.0249408972757124,  // c[3]
                                                                -0.0067465221034624,  // c[4]
                                                                159.8327742876945763,  // c[5]
                                                                -12.1694917209143032,  // c[6]
                                                                -5.3570145642206768,  // c[7]
                                                                -1.5728253671017924,  // c[8]
                                                                -0.3041926641251820,  // c[9]
                                                                88.2979533251867537,  // c[10]
                                                                -40.4637358159742817,  // c[11]
                                                                15.5471338214211485,  // c[12]
                                                                0.4450856161491892,  // c[13]
                                                                -3.4740850475291865,  // c[14]
                                                                37.2484054673747238,  // c[15]
                                                                -33.4889060032852157,  // c[16]
                                                                -11.0829828621952302,  // c[17]
                                                                -1.5689821797093835,  // c[18]
                                                                0.4240238530865685,  // c[19]
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
