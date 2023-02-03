
#ifndef __pulsating_fin_depth_cheby_hpp
#define __pulsating_fin_depth_cheby_hpp


// Include general usage libraries
#include <cassert>
#include <iostream>

// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "../config.hpp"
#include "../math/chebyshev.hpp"
#include "./fin_depth_coeffs/L1.hpp"


class P3
{
public:
    const cusfloat*     c;
    cusfloat*           cf;
    int                 current_inter = -1;
    const cusfloat*     intervals_bounds;
    const int*          nx;
    int*                nxf;
    const int*          ny;
    int*                nyf;
    const int*          nz;
    const int           num_intervals = 0;
    const int*          num_points;
    int                 num_points_f;
    const int*          num_points_cum;
    const bool*         x_log_scale;
    const cusfloat*     x_max;
    cusfloat*           x_map_scale;
    cusfloat*           x_map_scale_log;
    const cusfloat*     x_min;
    cusfloat*            x_min_l10;
    const bool*         y_log_scale;
    const cusfloat*     y_max;
    cusfloat*           y_map_scale;
    cusfloat*           y_map_scale_log;
    const cusfloat*     y_min;
    cusfloat*           y_min_l10;
    const bool*         z_log_scale;
    const cusfloat*     z_max;
    cusfloat*           z_map_scale;
    cusfloat*            z_map_scale_log;
    const cusfloat*     z_min;
    cusfloat*           z_min_l10;

    void initialize(void)
    {
        // std::cout << "Constructor working..." << std::endl;
        // Calculate map scale coefficients
        for (int i=0; i<L1C::num_intervals; i++)
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

            // Calculate for z coordinate
            this->z_map_scale[i] = 2.0/(this->z_max[i]-this->z_min[i]);
            this->z_map_scale_log[i] = 2.0/(
                                            std::log10(this->z_max[i])
                                            -
                                            std::log10(this->z_min[i])
                                            );
            this->z_min_l10[i] = std::log10(this->z_min[i]);
            
        }


    };

    void fold_h(cusfloat H)
    {
        // Get H working interval
        this->current_inter = this->get_interval_h(H);

        // Fold coefficients
        cusfloat zeta = this->z_map(H);
        int count = 0;
        int first_index = this->num_points_cum[this->current_inter];
        cusfloat coeff_i = this->c[first_index]*chebyshev_poly(this->nz[first_index], zeta);
        for (int j=this->num_points_cum[this->current_inter]+1; j<this->num_points_cum[this->current_inter+1]; j++)
        {
            // Check if there is a change in z branch
            if (this->ny[j] != this->ny[j-1])
            {
                // Storage data of the previous coefficient
                this->cf[count] = coeff_i;
                this->nxf[count] = this->nx[j-1];
                this->nyf[count] = this->ny[j-1];
                count++;

                // Restart coefficient for the current index
                coeff_i = this->c[j]*chebyshev_poly(this->nz[j], zeta);
            }
            else
            {
                coeff_i += this->c[j]*chebyshev_poly(this->nz[j], zeta);
            }
        }
        // Save last coefficient and update count to be the total
        // number of points (not the position in the vector) and also 
        // the first index for the next interval
        this->cf[count] = coeff_i;
        this->nxf[count] = this->nx[this->num_points_cum[this->current_inter+1]-1];
        this->nyf[count] = this->ny[this->num_points_cum[this->current_inter+1]-1];
        count++;

        // Storage intervals length
        this->num_points_f = count;

    }
    
    int get_interval_h(cusfloat H)
    {
        cusfloat d0=0.0, d1=0.0;
        int inter = -1;
        for (int i=0; i<this->num_intervals; i++)
        {
            // Calculate distances to bounds
            d0 = H-this->intervals_bounds[i];
            d1 = H-this->intervals_bounds[i+1];

            // Check if the value is in the interval
            if ((d0*d1) <= 0)
            {
                inter = i;
                break;
            }
        }

        // Check if the interval could be found
        assert(inter != -1 && "Could not be found an interval for H during folding");

        return inter;
    }

    cusfloat get_value_ab(cusfloat a, cusfloat b)
    {
        cusfloat xi = this->x_map(a);
        cusfloat eta = this->y_map(b);
        // std::cout << "xi: " << xi << std::endl;
        // std::cout << "eta: " << eta << std::endl;
        cusfloat sol = 0.0;
        for (int i=0; i<this->num_points_f; i++)
        {
            sol += this->cf[i]*chebyshev_poly_raw(this->nxf[i], xi)*chebyshev_poly_raw(this->nyf[i], eta);
        }

        return sol;
    }

    cusfloat get_value_abh(cusfloat a, cusfloat b, cusfloat h)
    {
        cusfloat sol = 0.0;
        for (int i=this->num_points_cum[3]; i<this->num_points_cum[4]; i++)
        {
            sol += (
                    this->c[i]
                    *
                    chebyshev_poly(this->nx[i], a)
                    *
                    chebyshev_poly(this->ny[i], b)
                    *
                    chebyshev_poly(this->nz[i], h)
                    );
        }

        return sol;
    }

    cusfloat x_map(cusfloat x)
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

    cusfloat y_map(cusfloat y)
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

    cusfloat z_map(cusfloat z)
    {
        cusfloat zi = 0.0;
        if (this->z_log_scale[this->current_inter])
        {
            zi = (
                this->z_map_scale_log[this->current_inter]
                *
                (std::log10(z)-this->z_min_l10[this->current_inter])
                -
                1.0
                );
        }
        else
        {
            zi = (
                this->z_map_scale[this->current_inter]
                *
                (z-this->z_min[this->current_inter])
                -
                1.0
            );
        }

        return zi;
    }

};

#endif