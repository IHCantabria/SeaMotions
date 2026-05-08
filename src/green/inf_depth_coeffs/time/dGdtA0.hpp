#ifndef __dGdtA0_coeffs_hpp
#define __dGdtA0_coeffs_hpp

#include "../../config.hpp"

struct dGdtA0C
{
static constexpr int                                        max_ref_level = 1;
static constexpr int                                        intervals_np = 2;
static constexpr int                                        max_cheby_order = 0;
static constexpr int                                        max_cheby_order_f = 0;
static constexpr std::size_t                                blocks_start[2] = {0, 1};
static constexpr std::size_t                                blocks_start_fall[2] = {0, 0};
static constexpr std::size_t                                blocks_coeffs_np[2] = {1, 1};
static constexpr std::size_t                                blocks_coeffs_np_fall[2] = {1, 1};
static constexpr std::size_t                                blocks_max_cheby_order[2] = {0, 0};
static constexpr std::size_t                                blocks_max_cheby_order_fall[2] = {0, 0};
static constexpr cusfloat                                   x_min_region[2] = {-4.000E+00, -2.000E+00};
static constexpr cusfloat                                   x_max_region[2] = {-2.000E+00, -8.687E-05};
static constexpr cusfloat                                   dx_region[2] = {2.000E+00, 2.000E+00};
static constexpr bool                                       fcn_log_scale = 0;
static constexpr bool                                       x_log_scale = 0;
static constexpr cusfloat                                   x_max_global = -0.0000868675834286;
static constexpr cusfloat                                   x_min_global = -4.0000000000000000;
static constexpr cusfloat                                   dx_min_region = 1.9999565662082857;

static constexpr std::size_t                                num_c = 2;
alignas(FLOATING_PRECISION)  static constexpr cusfloat      c[2] = {
                                                                1.0000002261025218,  // c[0]
                                                                1.0000002261025218,  // c[1]
                                                  };
alignas(FLOATING_PRECISION)  static constexpr std::size_t   ncx[2] = {
                                                                0,  // ncx[0]
                                                                0,  // ncx[1]
                                                  };
};
#endif
