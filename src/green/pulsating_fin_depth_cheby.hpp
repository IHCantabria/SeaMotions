
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
#include "./fin_depth_coeffs/L1_dA.hpp"
#include "./fin_depth_coeffs/L1_dB.hpp"
#include "./fin_depth_coeffs/L2.hpp"
#include "./fin_depth_coeffs/L3.hpp"
#include "./fin_depth_coeffs/M1.hpp"
#include "./fin_depth_coeffs/M2.hpp"
#include "./fin_depth_coeffs/M3.hpp"


class P3
{
public:
    cusfloat*       c;
    cusfloat*       cf;
    int             current_inter = -1;
    int             dims = -1;
    cusfloat        int_1d = 0.0;
    cusfloat*       intervals_bounds;
    int*            nx;
    int*            nxf;
    int*            ny;
    int*            nyf;
    int*            nz;
    int             num_intervals = 0;
    int*            num_points;
    int             num_points_f = 0;
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
    bool*           z_log_scale;
    cusfloat*       z_max;
    cusfloat*       z_map_scale;
    cusfloat*       z_map_scale_log;
    cusfloat*       z_min;
    cusfloat*       z_min_l10;

    void calculate_h_1D(cusfloat H)
    {
        // Get working interval
        this->current_inter = this->get_interval_h(H);
        std::cout << "current_interval: " << this->current_inter << std::endl;
        // Calcualte Integral value
        cusfloat zeta = this->x_map(H);
        std::cout << "zeta: " << zeta << std::endl;
        cusfloat coeff_i = 0.0;
        for (int j=this->num_points_cum[this->current_inter]; j<this->num_points_cum[this->current_inter+1]; j++)
        {
            coeff_i += this->c[j]*chebyshev_poly(this->nx[j], zeta);
        }

        // Storage intergral value
        this->int_1d = coeff_i;

    }

    void initialize(void)
    {
        assert(this->dims == 1 && "Dims value has not been loaded!\n");

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
            
            if (this->dims >= 2)
            {
                // Calculate for y coordinate
                this->y_map_scale[i] = 2.0/(this->y_max[i]-this->y_min[i]);
                this->y_map_scale_log[i] = 2.0/(
                                                std::log10(this->y_max[i])
                                                -
                                                std::log10(this->y_min[i])
                                                );
                this->y_min_l10[i] = std::log10(this->y_min[i]);

                if (this->dims == 3)
                {
                    // Calculate for z coordinate
                    this->z_map_scale[i] = 2.0/(this->z_max[i]-this->z_min[i]);
                    this->z_map_scale_log[i] = 2.0/(
                                                    std::log10(this->z_max[i])
                                                    -
                                                    std::log10(this->z_min[i])
                                                    );
                    this->z_min_l10[i] = std::log10(this->z_min[i]);
                }

            }
            
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


void set_data_l1(P3 *l1)
{
    // Storage data into target object
    l1->c                =   L1C::c;
    l1->cf               =   L1C::cf;
    l1->current_inter    =   -1;
    l1->dims             =   L1C::dims;
    l1->intervals_bounds =   L1C::interval_bounds;
    l1->nx               =   L1C::ncx;
    l1->nxf              =   L1C::ncxf;
    l1->ny               =   L1C::ncy;
    l1->nyf              =   L1C::ncyf;
    l1->nz               =   L1C::ncz;
    l1->num_intervals    =   L1C::num_intervals;
    l1->num_points       =   L1C::num_points;
    l1->num_points_f     =   0;
    l1->num_points_cum   =   L1C::num_points_cum;
    l1->x_log_scale      =   L1C::x_log_scale;
    l1->x_max            =   L1C::x_max;
    l1->x_map_scale      =   L1C::x_map_scale;
    l1->x_map_scale_log  =   L1C::x_map_scale_log;
    l1->x_min            =   L1C::x_min;
    l1->x_min_l10        =   L1C::x_min_l10;
    l1->y_log_scale      =   L1C::y_log_scale;
    l1->y_max            =   L1C::y_max;
    l1->y_map_scale      =   L1C::y_map_scale;
    l1->y_map_scale_log  =   L1C::y_map_scale_log;
    l1->y_min            =   L1C::y_min;
    l1->y_min_l10        =   L1C::y_min_l10;
    l1->z_log_scale      =   L1C::z_log_scale;
    l1->z_max            =   L1C::z_max;
    l1->z_map_scale      =   L1C::z_map_scale;
    l1->z_map_scale_log  =   L1C::z_map_scale_log;
    l1->z_min            =   L1C::z_min;
    l1->z_min_l10        =   L1C::z_min_l10;

    // Initialize object
    l1->initialize();
}


void set_data_l1_da(P3 *l1_da)
{
    // Storage data into target object
    l1_da->c                =   L1_dAC::c;
    l1_da->cf               =   L1_dAC::cf;
    l1_da->current_inter    =   -1;
    l1_da->dims             =   L1_dAC::dims;
    l1_da->intervals_bounds =   L1_dAC::interval_bounds;
    l1_da->nx               =   L1_dAC::ncx;
    l1_da->nxf              =   L1_dAC::ncxf;
    l1_da->ny               =   L1_dAC::ncy;
    l1_da->nyf              =   L1_dAC::ncyf;
    l1_da->nz               =   L1_dAC::ncz;
    l1_da->num_intervals    =   L1_dAC::num_intervals;
    l1_da->num_points       =   L1_dAC::num_points;
    l1_da->num_points_f     =   0;
    l1_da->num_points_cum   =   L1_dAC::num_points_cum;
    l1_da->x_log_scale      =   L1_dAC::x_log_scale;
    l1_da->x_max            =   L1_dAC::x_max;
    l1_da->x_map_scale      =   L1_dAC::x_map_scale;
    l1_da->x_map_scale_log  =   L1_dAC::x_map_scale_log;
    l1_da->x_min            =   L1_dAC::x_min;
    l1_da->x_min_l10        =   L1_dAC::x_min_l10;
    l1_da->y_log_scale      =   L1_dAC::y_log_scale;
    l1_da->y_max            =   L1_dAC::y_max;
    l1_da->y_map_scale      =   L1_dAC::y_map_scale;
    l1_da->y_map_scale_log  =   L1_dAC::y_map_scale_log;
    l1_da->y_min            =   L1_dAC::y_min;
    l1_da->y_min_l10        =   L1_dAC::y_min_l10;
    l1_da->z_log_scale      =   L1_dAC::z_log_scale;
    l1_da->z_max            =   L1_dAC::z_max;
    l1_da->z_map_scale      =   L1_dAC::z_map_scale;
    l1_da->z_map_scale_log  =   L1_dAC::z_map_scale_log;
    l1_da->z_min            =   L1_dAC::z_min;
    l1_da->z_min_l10        =   L1_dAC::z_min_l10;

    // Initialize object
    l1_da->initialize();
}


void set_data_l1_db(P3 *l1_db)
{
    // Storage data into target object
    l1_db->c                =   L1_dBC::c;
    l1_db->cf               =   L1_dBC::cf;
    l1_db->current_inter    =   -1;
    l1_db->dims             =   L1_dBC::dims;
    l1_db->intervals_bounds =   L1_dBC::interval_bounds;
    l1_db->nx               =   L1_dBC::ncx;
    l1_db->nxf              =   L1_dBC::ncxf;
    l1_db->ny               =   L1_dBC::ncy;
    l1_db->nyf              =   L1_dBC::ncyf;
    l1_db->nz               =   L1_dBC::ncz;
    l1_db->num_intervals    =   L1_dBC::num_intervals;
    l1_db->num_points       =   L1_dBC::num_points;
    l1_db->num_points_f     =   0;
    l1_db->num_points_cum   =   L1_dBC::num_points_cum;
    l1_db->x_log_scale      =   L1_dBC::x_log_scale;
    l1_db->x_max            =   L1_dBC::x_max;
    l1_db->x_map_scale      =   L1_dBC::x_map_scale;
    l1_db->x_map_scale_log  =   L1_dBC::x_map_scale_log;
    l1_db->x_min            =   L1_dBC::x_min;
    l1_db->x_min_l10        =   L1_dBC::x_min_l10;
    l1_db->y_log_scale      =   L1_dBC::y_log_scale;
    l1_db->y_max            =   L1_dBC::y_max;
    l1_db->y_map_scale      =   L1_dBC::y_map_scale;
    l1_db->y_map_scale_log  =   L1_dBC::y_map_scale_log;
    l1_db->y_min            =   L1_dBC::y_min;
    l1_db->y_min_l10        =   L1_dBC::y_min_l10;
    l1_db->z_log_scale      =   L1_dBC::z_log_scale;
    l1_db->z_max            =   L1_dBC::z_max;
    l1_db->z_map_scale      =   L1_dBC::z_map_scale;
    l1_db->z_map_scale_log  =   L1_dBC::z_map_scale_log;
    l1_db->z_min            =   L1_dBC::z_min;
    l1_db->z_min_l10        =   L1_dBC::z_min_l10;

    // Initialize object
    l1_db->initialize();
}


void set_data_l2(P3 *l2)
{
    // Storage data into target object
    l2->c                =   L2C::c;
    l2->cf               =   L2C::cf;
    l2->current_inter    =   -1;
    l2->dims             =   L2C::dims;
    l2->intervals_bounds =   L2C::interval_bounds;
    l2->nx               =   L2C::ncx;
    l2->nxf              =   L2C::ncxf;
    l2->num_intervals    =   L2C::num_intervals;
    l2->num_points       =   L2C::num_points;
    l2->num_points_f     =   0;
    l2->num_points_cum   =   L2C::num_points_cum;
    l2->x_log_scale      =   L2C::x_log_scale;
    l2->x_max            =   L2C::x_max;
    l2->x_map_scale      =   L2C::x_map_scale;
    l2->x_map_scale_log  =   L2C::x_map_scale_log;
    l2->x_min            =   L2C::x_min;
    l2->x_min_l10        =   L2C::x_min_l10;

    // Initialize object
    l2->initialize();
}


void set_data_l3(P3 *l3)
{
    // Storage data into target object
    l3->c                =   L3C::c;
    l3->cf               =   L3C::cf;
    l3->current_inter    =   -1;
    l3->dims             =   L3C::dims;
    l3->intervals_bounds =   L3C::interval_bounds;
    l3->nx               =   L3C::ncx;
    l3->nxf              =   L3C::ncxf;
    l3->ny               =   L3C::ncy;
    l3->nyf              =   L3C::ncyf;
    l3->nz               =   L3C::ncz;
    l3->num_intervals    =   L3C::num_intervals;
    l3->num_points       =   L3C::num_points;
    l3->num_points_f     =   0;
    l3->num_points_cum   =   L3C::num_points_cum;
    l3->x_log_scale      =   L3C::x_log_scale;
    l3->x_max            =   L3C::x_max;
    l3->x_map_scale      =   L3C::x_map_scale;
    l3->x_map_scale_log  =   L3C::x_map_scale_log;
    l3->x_min            =   L3C::x_min;
    l3->x_min_l10        =   L3C::x_min_l10;
    l3->y_log_scale      =   L3C::y_log_scale;
    l3->y_max            =   L3C::y_max;
    l3->y_map_scale      =   L3C::y_map_scale;
    l3->y_map_scale_log  =   L3C::y_map_scale_log;
    l3->y_min            =   L3C::y_min;
    l3->y_min_l10        =   L3C::y_min_l10;
    l3->z_log_scale      =   L3C::z_log_scale;
    l3->z_max            =   L3C::z_max;
    l3->z_map_scale      =   L3C::z_map_scale;
    l3->z_map_scale_log  =   L3C::z_map_scale_log;
    l3->z_min            =   L3C::z_min;
    l3->z_min_l10        =   L3C::z_min_l10;

    // Initialize object
    l3->initialize();
}


void set_data_m1(P3 *m1)
{
    // Storage data into target object
    m1->c                =   M1C::c;
    m1->cf               =   M1C::cf;
    m1->current_inter    =   -1;
    m1->dims             =   M1C::dims;
    m1->intervals_bounds =   M1C::interval_bounds;
    m1->nx               =   M1C::ncx;
    m1->nxf              =   M1C::ncxf;
    m1->ny               =   M1C::ncy;
    m1->nyf              =   M1C::ncyf;
    m1->nz               =   M1C::ncz;
    m1->num_intervals    =   M1C::num_intervals;
    m1->num_points       =   M1C::num_points;
    m1->num_points_f     =   0;
    m1->num_points_cum   =   M1C::num_points_cum;
    m1->x_log_scale      =   M1C::x_log_scale;
    m1->x_max            =   M1C::x_max;
    m1->x_map_scale      =   M1C::x_map_scale;
    m1->x_map_scale_log  =   M1C::x_map_scale_log;
    m1->x_min            =   M1C::x_min;
    m1->x_min_l10        =   M1C::x_min_l10;
    m1->y_log_scale      =   M1C::y_log_scale;
    m1->y_max            =   M1C::y_max;
    m1->y_map_scale      =   M1C::y_map_scale;
    m1->y_map_scale_log  =   M1C::y_map_scale_log;
    m1->y_min            =   M1C::y_min;
    m1->y_min_l10        =   M1C::y_min_l10;
    m1->z_log_scale      =   M1C::z_log_scale;
    m1->z_max            =   M1C::z_max;
    m1->z_map_scale      =   M1C::z_map_scale;
    m1->z_map_scale_log  =   M1C::z_map_scale_log;
    m1->z_min            =   M1C::z_min;
    m1->z_min_l10        =   M1C::z_min_l10;

    // Initialize object
    m1->initialize();
}


void set_data_m2(P3 *m2)
{
    // Storage data into target object
    m2->c                =   M2C::c;
    m2->cf               =   M2C::cf;
    m2->current_inter    =   -1;
    m2->dims             =   M2C::dims;
    m2->intervals_bounds =   M2C::interval_bounds;
    m2->nx               =   M2C::ncx;
    m2->nxf              =   M2C::ncxf;
    m2->num_intervals    =   M2C::num_intervals;
    m2->num_points       =   M2C::num_points;
    m2->num_points_f     =   0;
    m2->num_points_cum   =   M2C::num_points_cum;
    m2->x_log_scale      =   M2C::x_log_scale;
    m2->x_max            =   M2C::x_max;
    m2->x_map_scale      =   M2C::x_map_scale;
    m2->x_map_scale_log  =   M2C::x_map_scale_log;
    m2->x_min            =   M2C::x_min;
    m2->x_min_l10        =   M2C::x_min_l10;

    // Initialize object
    m2->initialize();
}


void set_data_m3(P3 *m3)
{
    // Storage data into target object
    m3->c                =   M3C::c;
    m3->cf               =   M3C::cf;
    m3->current_inter    =   -1;
    m3->dims             =   M3C::dims;
    m3->intervals_bounds =   M3C::interval_bounds;
    m3->nx               =   M3C::ncx;
    m3->nxf              =   M3C::ncxf;
    m3->ny               =   M3C::ncy;
    m3->nyf              =   M3C::ncyf;
    m3->nz               =   M3C::ncz;
    m3->num_intervals    =   M3C::num_intervals;
    m3->num_points       =   M3C::num_points;
    m3->num_points_f     =   0;
    m3->num_points_cum   =   M3C::num_points_cum;
    m3->x_log_scale      =   M3C::x_log_scale;
    m3->x_max            =   M3C::x_max;
    m3->x_map_scale      =   M3C::x_map_scale;
    m3->x_map_scale_log  =   M3C::x_map_scale_log;
    m3->x_min            =   M3C::x_min;
    m3->x_min_l10        =   M3C::x_min_l10;
    m3->y_log_scale      =   M3C::y_log_scale;
    m3->y_max            =   M3C::y_max;
    m3->y_map_scale      =   M3C::y_map_scale;
    m3->y_map_scale_log  =   M3C::y_map_scale_log;
    m3->y_min            =   M3C::y_min;
    m3->y_min_l10        =   M3C::y_min_l10;
    m3->z_log_scale      =   M3C::z_log_scale;
    m3->z_max            =   M3C::z_max;
    m3->z_map_scale      =   M3C::z_map_scale;
    m3->z_map_scale_log  =   M3C::z_map_scale_log;
    m3->z_min            =   M3C::z_min;
    m3->z_min_l10        =   M3C::z_min_l10;

    // Initialize object
    m3->initialize();
}

#endif
