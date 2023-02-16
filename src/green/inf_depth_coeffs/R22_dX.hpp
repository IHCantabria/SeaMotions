#ifndef __R22_dX_coeffs_hpp
#define __R22_dX_coeffs_hpp

#include "../../config.hpp"

namespace R22_dXC
{
    extern int          dims;
    extern int          num_intervals_x;
    extern int          num_intervals_y;
    extern cusfloat     interval_bounds_x[];
    extern cusfloat     interval_bounds_y[];
    extern int          num_points[];
    extern int          num_points_cum[];
    extern bool         x_log_scale[];
    extern cusfloat     x_map_scale[];
    extern cusfloat     x_map_scale_log[];
    extern cusfloat     x_max[];
    extern cusfloat     x_min[];
    extern cusfloat     x_min_l10[];
    extern bool         y_log_scale[];
    extern cusfloat     y_map_scale[];
    extern cusfloat     y_map_scale_log[];
    extern cusfloat     y_max[];
    extern cusfloat     y_min[];
    extern cusfloat     y_min_l10[];

    extern int          num_c;
    extern cusfloat     c[];
    extern int          ncx[];
    extern int          ncy[];
}
#endif
