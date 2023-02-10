#ifndef __L1_dA_coeffs_hpp
#define __L1_dA_coeffs_hpp

#include "../../config.hpp"

namespace L1_dAC
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
    extern bool         y_log_scale[];
    extern cusfloat     y_map_scale[];
    extern cusfloat     y_map_scale_log[];
    extern cusfloat     y_max[];
    extern cusfloat     y_min[];
    extern cusfloat     y_min_l10[];
    extern bool         z_log_scale[];
    extern cusfloat     z_map_scale[];
    extern cusfloat     z_map_scale_log[];
    extern cusfloat     z_max[];
    extern cusfloat     z_min[];
    extern cusfloat     z_min_l10[];

    extern int          num_c;
    extern cusfloat     c[];
    extern cusfloat     cf[];
    extern int          ncx[];
    extern int          ncxf[];
    extern int          ncy[];
    extern int          ncyf[];
    extern int          ncz[];
}
#endif
