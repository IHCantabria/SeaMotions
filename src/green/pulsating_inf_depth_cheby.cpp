
// Include general usage libraries
#include <cassert>

// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "pulsating_inf_depth_cheby.hpp"


void P2::initialize(void)
{
    assert(this->is_build != -1 && "Database not loaded into P2!\n");

    // Calculate map scale coefficients
    for (int i=0; i<this->num_intervals; i++)
    {
        // Calculate for x coordinate
        this->x_map_scale[i] = 2.0/(this->x_max[i]-this->x_min[i]);
        this->x_map_scale_log[i] = 2.0/(
                                        std::log10(this->x_max[i])
                                        -
                                        std::log10(this->x_min[i])
                                        );
        this->x_min_l10[i] = std::log10(this->x_min[i]);
        
        // Calculate for y coordinate
        this->y_map_scale[i] = 2.0/(this->y_max[i]-this->y_min[i]);
        this->y_map_scale_log[i] = 2.0/(
                                        std::log10(this->y_max[i])
                                        -
                                        std::log10(this->y_min[i])
                                        );
        this->y_min_l10[i] = std::log10(this->y_min[i]);
    }

};


cusfloat P2::x_map(cusfloat x)
{
    cusfloat xi = 0.0;
    if (this->x_log_scale[this->current_inter])
    {
        xi = (
            this->x_map_scale_log[this->current_inter]
            *
            (std::log10(x)-this->x_min_l10[this->current_inter])
            -
            1.0
            );
    }
    else
    {
        xi = (
            this->x_map_scale[this->current_inter]
            *
            (x-this->x_min[this->current_inter])
            -
            1.0
            );
    }

    return xi;
}


cusfloat P2::y_map(cusfloat y)
{
    cusfloat yi = 0.0;
    if (this->y_log_scale[this->current_inter])
    {
        yi = (
            this->y_map_scale_log[this->current_inter]
            *
            (std::log10(y)-this->y_min_l10[this->current_inter])
            -
            1.0
            );
    }
    else
    {
        yi = (
            this->y_map_scale[this->current_inter]
            *
            (y-this->y_min[this->current_inter])
            -
            1.0
            );
    }

    return yi;
}