#ifndef __M2_coeffs_hpp
#define __M2_coeffs_hpp

#include "../../config.hpp"

namespace M2C
{
    extern int          dims;
    extern int          num_intervals;
    extern cusfloat     interval_bounds[];
    extern int          num_points[];
    extern int          num_points_cum[];
    extern int          max_size_fold;
    extern bool         x_log_scale[];
    extern cusfloat     x_map_scale[];
    extern cusfloat     x_map_scale_log[];
    extern cusfloat     x_max[];
    extern cusfloat     x_min[];
    extern cusfloat     x_min_l10[];

    extern int          num_c;
    extern cusfloat     c[];
    extern cusfloat     cf[];
    extern cusfloat     cf2[];
    extern int          ncx[];
    extern int          ncxf[];
    extern int          ncxf2[];
}
#endif
