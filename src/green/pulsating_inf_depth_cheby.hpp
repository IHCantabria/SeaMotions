
#ifndef __pulsating_inf_depth_cheby.hpp
#define __pulsating_inf_depth_cheby.hpp

// Include local modules
#include "../config.hpp"


class P2
{
public:
    cusfloat*       c;
    int             current_inter = -1;
    cusfloat*       intervals_bounds;
    bool            is_build = false;
    int*            nx;
    int*            ny;
    int             num_intervals = 0;
    int*            num_points;
    int*            num_points_cum;
    bool*           x_log_scale;
    cusfloat*       x_max;
    cusfloat*       x_map_scale;
    cusfloat*       x_map_scale_log;
    cusfloat*       x_min;
    cusfloat*       x_min_l10;
    bool*           y_log_scale;
    cusfloat*       y_max;
    cusfloat*       y_map_scale;
    cusfloat*       y_map_scale_log;
    cusfloat*       y_min;
    cusfloat*       y_min_l10;

    void initialize(void);
    cusfloat x_map(cusfloat x);
    cusfloat y_map(cusfloat y);

};

#endif