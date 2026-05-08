#ifndef __dGdttxA2_coeffs_hpp
#define __dGdttxA2_coeffs_hpp

#include "../../config.hpp"

struct dGdttxA2C
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
                                                                69.1148119766277915,  // c[0]
                                                                -1.6185680866057810,  // c[1]
                                                                -0.7394942866124339,  // c[2]
                                                                -0.2357449871409081,  // c[3]
                                                                -0.0550383392068561,  // c[4]
                                                                46.5784985860591547,  // c[5]
                                                                -21.7940994295912205,  // c[6]
                                                                0.4064158183434117,  // c[7]
                                                                2.9927566530785743,  // c[8]
                                                                0.4015705531296820,  // c[9]
                                                                27.2242027705335445,  // c[10]
                                                                -1.1694355940854266,  // c[11]
                                                                -0.4549559084238517,  // c[12]
                                                                -0.1248353629794785,  // c[13]
                                                                -0.0167663066475860,  // c[14]
                                                                22.6221278218618380,  // c[15]
                                                                1.5268560493969898,  // c[16]
                                                                8.9272576971509050,  // c[17]
                                                                7.2501033343335557,  // c[18]
                                                                3.8952532136661695,  // c[19]
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
