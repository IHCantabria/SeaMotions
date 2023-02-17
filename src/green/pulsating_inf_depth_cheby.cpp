
// Include general usage libraries
#include <cassert>

// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "./inf_depth_coeffs/R11.hpp"
#include "./inf_depth_coeffs/R11A_dX.hpp"
#include "./inf_depth_coeffs/R11B_dX.hpp"
#include "./inf_depth_coeffs/R12.hpp"
#include "./inf_depth_coeffs/R12_dX.hpp"
#include "./inf_depth_coeffs/R21.hpp"
#include "./inf_depth_coeffs/R21_dX.hpp"
#include "./inf_depth_coeffs/R22.hpp"
#include "./inf_depth_coeffs/R22_dX.hpp"
#include "../math/chebyshev.hpp"
#include "../math/math_tools.hpp"
#include "../math/special_math.hpp"
#include "pulsating_inf_depth_cheby.hpp"


// Include some namespaces into the local one
using namespace std;


int P2::calculate_interval(cusfloat x, cusfloat y)
{
    // Define auxiliar local variables
    cusfloat d0=0.0, d1=0.0;

    // Calculate x interval
    int inter_x = -1;
    if (x > this->intervals_bounds_x[this->num_intervals_x])
    {
        inter_x = this->num_intervals_x-1;
    }
    else if (x < this->intervals_bounds_x[0])
    {
        inter_x = 0;
    }
    else
    {
        for (int i=0; i<this->num_intervals_x; i++)
        {
            // Calculate distances to bounds
            d0 = x-this->intervals_bounds_x[i];
            d1 = x-this->intervals_bounds_x[i+1];

            // Check if the value is in the interval
            if ((d0*d1) <= 0)
            {
                inter_x = i;
                break;
            }
        }

        // Check if the interval could be found
        assert(inter_x != -1 && "Could not be found an interval for H during folding");

    }

    // Calculate y interval
    int inter_y = -1;
    if (y > this->intervals_bounds_y[this->num_intervals_y])
    {
        inter_y = this->num_intervals_y-1;
    }
    else if (y < this->intervals_bounds_y[0])
    {
        inter_y = 0;
    }
    else
    {
        for (int i=0; i<this->num_intervals_y; i++)
        {
            // Calculate distances to bounds
            d0 = y-this->intervals_bounds_y[i];
            d1 = y-this->intervals_bounds_y[i+1];

            // Check if the value is in the interval
            if ((d0*d1) <= 0)
            {
                inter_y = i;
                break;
            }
        }

        // Check if the interval could be found
        // std::cout << "H: " << std::endl;
        assert(inter_y != -1 && "Could not be found an interval for H during folding");

    }

    // Calculate global interval
    int inter = inter_x*this->num_intervals_y + inter_y;

    return inter;
}


cusfloat P2::calculate_xy_cheby(cusfloat x, cusfloat y)
{
    // Get current interval
    int inter = this->calculate_interval(x, y);

    // Map variables from real parametric space to 
    // approximation space
    cusfloat xi = this->x_map(x);
    cusfloat eta = this->y_map(y);

    // Loop over coefficients to calculate integral value
    cusfloat sol = 0.0;
    for (int i=this->num_points_cum[inter]; i<this->num_points_cum[inter+1]; i++)
    {
        sol += this->c[i]*chebyshev_poly_raw(this->nx[i], xi)*chebyshev_poly_raw(this->ny[i], eta);
    }

    return sol;
}


cusfloat P2::fit_boundary_x(cusfloat x)
{
    cusfloat xf = 0.0;
    if (x < this->intervals_bounds_x[0])
    {
        xf = this->intervals_bounds_x[0];
    }
    else if (x > this->intervals_bounds_x[this->num_intervals_x-1])
    {
        xf = this->intervals_bounds_x[this->num_intervals_x-1];
    }
    else
    {
        xf = x;
    }
    
    return xf;
}


cusfloat P2::fit_boundary_y(cusfloat y)
{
    cusfloat yf = 0.0;
    if (y < this->intervals_bounds_y[0])
    {
        yf = this->intervals_bounds_y[0];
    }
    else if (y > this->intervals_bounds_y[this->num_intervals_y-1])
    {
        yf = this->intervals_bounds_y[this->num_intervals_y-1];
    }
    else
    {
        yf = y;
    }
    
    return yf;
}


void P2::initialize(void)
{
    assert(this->is_build && "Database not loaded into P2!\n");

    // Calculate map scale coefficients
    for (int i=0; i<this->num_intervals_x; i++)
    {
        // Calculate for x coordinate
        this->x_map_scale[i] = 2.0/(this->x_max[i]-this->x_min[i]);
        this->x_map_scale_log[i] = 2.0/(
                                        std::log10(this->x_max[i])
                                        -
                                        std::log10(this->x_min[i])
                                        );
        this->x_min_l10[i] = std::log10(this->x_min[i]);
        
    }

    for (int i=0; i<this->num_intervals_y; i++)
    {
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


cusfloat R11::calculate_xy(cusfloat x, cusfloat y)
{
    // Fit input values to boundaries
    cusfloat xf = this->fit_boundary_x(x);
    cusfloat yf = this->fit_boundary_y(y);

    // Calculate chebyshev part
    cusfloat cheby_sol = this->calculate_xy_cheby(xf, yf);

    // Calculate total solution
    cusfloat R = sqrt(pow2s(xf) + pow2s(yf));
    cusfloat c0 = (
                    -2*besselj0(xf)*log(R+yf)
                    -(PI*bessely0(xf)-2*besselj0(xf)*log(xf))
                    )*exp(-yf);
    cusfloat sol = exp(-yf)*R*cheby_sol + c0;

    return sol;
}


cusfloat R11A_dX::calculate_xy(cusfloat x, cusfloat y)
{
    // Fit input values to boundaries
    cusfloat xf = this->fit_boundary_x(x);
    cusfloat yf = this->fit_boundary_y(y);

    // Calculate chebyshev part
    cusfloat cheby_sol = this->calculate_xy_cheby(xf, yf);

    // Calculate total solution
    cusfloat sol = pow(10.0, cheby_sol);
    sol += PI*exp(-yf)*(bessely1(xf)+struve1(xf)-2.0/PI);

    return sol;
}


cusfloat R11B_dX::calculate_xy(cusfloat x, cusfloat y)
{
    // Fit input values to boundaries
    cusfloat xf = this->fit_boundary_x(x);
    cusfloat yf = this->fit_boundary_y(y);

    // Calculate chebyshev part
    cusfloat cheby_sol = this->calculate_xy_cheby(xf, yf);

    // Calculate total solution
    cusfloat R = sqrt(pow2s(xf) + pow2s(yf));
    cusfloat c0 = (
                    +2*besselj1(xf)*log(R+yf)
                    -2*besselj0(xf)*xf/R/(R+yf)
                    +PI*bessely1(xf)
                    -2*besselj1(xf)*log(xf)
                    +2*besselj0(xf)/xf
                    )*exp(-yf);
    cusfloat sol = exp(-yf)*R*cheby_sol + c0;

    return sol;
}


cusfloat R12::calculate_xy(cusfloat x, cusfloat y)
{
    // Fit input values to boundaries
    cusfloat xf = this->fit_boundary_x(x);
    cusfloat yf = this->fit_boundary_y(y);

    // Calculate chebyshev part
    cusfloat cheby_sol = this->calculate_xy_cheby(xf, yf);

    return cheby_sol;
}


cusfloat R12_dX::calculate_xy(cusfloat x, cusfloat y)
{
    // Fit input values to boundaries
    cusfloat xf = this->fit_boundary_x(x);
    cusfloat yf = this->fit_boundary_y(y);

    // Calculate chebyshev part
    cusfloat cheby_sol = this->calculate_xy_cheby(xf, yf);

    return cheby_sol;
}


cusfloat R21::calculate_xy(cusfloat x, cusfloat y)
{
    // Fit input values to boundaries
    cusfloat xf = this->fit_boundary_x(x);
    cusfloat yf = this->fit_boundary_y(y);

    // Calculate chebyshev part
    cusfloat cheby_sol = this->calculate_xy_cheby(xf, yf);

    // Calculate total solution
    cusfloat sol = cheby_sol - 2*PI*exp(-yf)*bessely0(xf);

    return sol;
}


cusfloat R21_dX::calculate_xy(cusfloat x, cusfloat y)
{
    // Fit input values to boundaries
    cusfloat xf = this->fit_boundary_x(x);
    cusfloat yf = this->fit_boundary_y(y);

    // Calculate chebyshev part
    cusfloat cheby_sol = this->calculate_xy_cheby(xf, yf);

    // Calculate total solution
    cusfloat sol = cheby_sol + 2*PI*exp(-yf)*bessely1(xf);

    return sol;
}


cusfloat R22::calculate_xy(cusfloat x, cusfloat y)
{
    // Fit input values to boundaries
    cusfloat xf = this->fit_boundary_x(x);
    cusfloat yf = this->fit_boundary_y(y);

    // Calculate chebyshev part
    cusfloat cheby_sol = this->calculate_xy_cheby(xf, yf);

    // Calculate total solution
    cusfloat sol = cheby_sol - 2*PI*exp(-yf)*bessely0(xf);
    cusfloat R = sqrt(pow2s(xf) + pow2s(yf));
    cusfloat pn = 1.0;
    cusfloat rf = 1/R;
    cusfloat rf2 = rf;
    cusfloat costh = yf*rf2;
    for (int i=0; i<3; i++)
    {
        // Add new series term
        sol += pn*legendre_poly_raw(i, costh)*rf;

        // Update series coefficients
        pn *= static_cast<cusfloat>(i+1);
        rf *= rf2;
    }

    return sol;
}


cusfloat R22_dX::calculate_xy(cusfloat x, cusfloat y)
{
    // Fit input values to boundaries
    cusfloat xf = this->fit_boundary_x(x);
    cusfloat yf = this->fit_boundary_y(y);

    // Calculate chebyshev part
    cusfloat cheby_sol = this->calculate_xy_cheby(xf, yf);

    // Calculate total solution
    cusfloat sol = cheby_sol + 2*PI*exp(-yf)*bessely1(xf);
    cusfloat R = sqrt(pow2s(xf) + pow2s(yf));
    cusfloat pn = 1.0;
    cusfloat rf = 1/R;
    cusfloat rf2 = rf;
    cusfloat costh = yf*rf2;
    cusfloat costhd = -xf*costh*pow2s(rf2);
    for (int i=0; i<3; i++)
    {
        // Add new series term
        sol += pn* (
                    legendre_poly_der_raw(i, xf)*costhd*rf
                    +
                    legendre_poly_raw(i, xf)*2*(i+1)*xf*rf*rf2
                    );

        // Update series coefficients
        pn *= static_cast<cusfloat>(i+1);
        rf *= rf2;
    }

    return sol;
}


void set_data_r11(R11* r11)
{
    // Load data into target object
    r11->c                  =   R11C::c;
    r11->intervals_bounds_x =   R11C::interval_bounds_x;
    r11->intervals_bounds_y =   R11C::interval_bounds_y;
    r11->num_c              =   R11C::num_c;
    r11->num_intervals_x    =   R11C::num_intervals_x;
    r11->num_intervals_y    =   R11C::num_intervals_y;
    r11->num_points         =   R11C::num_points;
    r11->num_points_cum     =   R11C::num_points_cum;
    r11->nx                 =   R11C::ncx;
    r11->ny                 =   R11C::ncy;
    r11->x_log_scale        =   R11C::x_log_scale;
    r11->x_map_scale        =   R11C::x_map_scale;
    r11->x_map_scale_log    =   R11C::x_map_scale_log;
    r11->x_max              =   R11C::x_max;
    r11->x_min              =   R11C::x_min;
    r11->x_min_l10          =   R11C::x_min_l10;
    r11->y_log_scale        =   R11C::y_log_scale;
    r11->y_map_scale        =   R11C::y_map_scale;
    r11->y_map_scale_log    =   R11C::y_map_scale_log;
    r11->y_max              =   R11C::y_max;
    r11->y_min              =   R11C::y_min;
    r11->y_min_l10          =   R11C::y_min_l10;
    r11->is_build           =   true;

    // Intialize object
    r11->initialize();

}


void set_data_r11a_dx(R11A_dX* r11a_dx)
{
    // Load data into target object
    r11a_dx->c                  =   R11A_dXC::c;
    r11a_dx->intervals_bounds_x =   R11A_dXC::interval_bounds_x;
    r11a_dx->intervals_bounds_y =   R11A_dXC::interval_bounds_y;
    r11a_dx->num_c              =   R11A_dXC::num_c;
    r11a_dx->num_intervals_x    =   R11A_dXC::num_intervals_x;
    r11a_dx->num_intervals_y    =   R11A_dXC::num_intervals_y;
    r11a_dx->num_points         =   R11A_dXC::num_points;
    r11a_dx->num_points_cum     =   R11A_dXC::num_points_cum;
    r11a_dx->nx                 =   R11A_dXC::ncx;
    r11a_dx->ny                 =   R11A_dXC::ncy;
    r11a_dx->x_log_scale        =   R11A_dXC::x_log_scale;
    r11a_dx->x_map_scale        =   R11A_dXC::x_map_scale;
    r11a_dx->x_map_scale_log    =   R11A_dXC::x_map_scale_log;
    r11a_dx->x_max              =   R11A_dXC::x_max;
    r11a_dx->x_min              =   R11A_dXC::x_min;
    r11a_dx->x_min_l10          =   R11A_dXC::x_min_l10;
    r11a_dx->y_log_scale        =   R11A_dXC::y_log_scale;
    r11a_dx->y_map_scale        =   R11A_dXC::y_map_scale;
    r11a_dx->y_map_scale_log    =   R11A_dXC::y_map_scale_log;
    r11a_dx->y_max              =   R11A_dXC::y_max;
    r11a_dx->y_min              =   R11A_dXC::y_min;
    r11a_dx->y_min_l10          =   R11A_dXC::y_min_l10;
    r11a_dx->is_build           =   true;

    // Intialize object
    r11a_dx->initialize();

}


void set_data_r11b_dx(R11B_dX* r11b_dx)
{
    // Load data into target object
    r11b_dx->c                  =   R11B_dXC::c;
    r11b_dx->intervals_bounds_x =   R11B_dXC::interval_bounds_x;
    r11b_dx->intervals_bounds_y =   R11B_dXC::interval_bounds_y;
    r11b_dx->num_c              =   R11B_dXC::num_c;
    r11b_dx->num_intervals_x    =   R11B_dXC::num_intervals_x;
    r11b_dx->num_intervals_y    =   R11B_dXC::num_intervals_y;
    r11b_dx->num_points         =   R11B_dXC::num_points;
    r11b_dx->num_points_cum     =   R11B_dXC::num_points_cum;
    r11b_dx->nx                 =   R11B_dXC::ncx;
    r11b_dx->ny                 =   R11B_dXC::ncy;
    r11b_dx->x_log_scale        =   R11B_dXC::x_log_scale;
    r11b_dx->x_map_scale        =   R11B_dXC::x_map_scale;
    r11b_dx->x_map_scale_log    =   R11B_dXC::x_map_scale_log;
    r11b_dx->x_max              =   R11B_dXC::x_max;
    r11b_dx->x_min              =   R11B_dXC::x_min;
    r11b_dx->x_min_l10          =   R11B_dXC::x_min_l10;
    r11b_dx->y_log_scale        =   R11B_dXC::y_log_scale;
    r11b_dx->y_map_scale        =   R11B_dXC::y_map_scale;
    r11b_dx->y_map_scale_log    =   R11B_dXC::y_map_scale_log;
    r11b_dx->y_max              =   R11B_dXC::y_max;
    r11b_dx->y_min              =   R11B_dXC::y_min;
    r11b_dx->y_min_l10          =   R11B_dXC::y_min_l10;
    r11b_dx->is_build           =   true;

    // Intialize object
    r11b_dx->initialize();

}


void set_data_r12(R12* r12)
{
    // Load data into target object
    r12->c                  =   R12C::c;
    r12->intervals_bounds_x =   R12C::interval_bounds_x;
    r12->intervals_bounds_y =   R12C::interval_bounds_y;
    r12->num_c              =   R12C::num_c;
    r12->num_intervals_x    =   R12C::num_intervals_x;
    r12->num_intervals_y    =   R12C::num_intervals_y;
    r12->num_points         =   R12C::num_points;
    r12->num_points_cum     =   R12C::num_points_cum;
    r12->nx                 =   R12C::ncx;
    r12->ny                 =   R12C::ncy;
    r12->x_log_scale        =   R12C::x_log_scale;
    r12->x_map_scale        =   R12C::x_map_scale;
    r12->x_map_scale_log    =   R12C::x_map_scale_log;
    r12->x_max              =   R12C::x_max;
    r12->x_min              =   R12C::x_min;
    r12->x_min_l10          =   R12C::x_min_l10;
    r12->y_log_scale        =   R12C::y_log_scale;
    r12->y_map_scale        =   R12C::y_map_scale;
    r12->y_map_scale_log    =   R12C::y_map_scale_log;
    r12->y_max              =   R12C::y_max;
    r12->y_min              =   R12C::y_min;
    r12->y_min_l10          =   R12C::y_min_l10;
    r12->is_build           =   true;

    // Intialize object
    r12->initialize();

}


void set_data_r12_dx(R12_dX* r12_dx)
{
    // Load data into target object
    r12_dx->c                  =   R12_dXC::c;
    r12_dx->intervals_bounds_x =   R12_dXC::interval_bounds_x;
    r12_dx->intervals_bounds_y =   R12_dXC::interval_bounds_y;
    r12_dx->num_c              =   R12_dXC::num_c;
    r12_dx->num_intervals_x    =   R12_dXC::num_intervals_x;
    r12_dx->num_intervals_y    =   R12_dXC::num_intervals_y;
    r12_dx->num_points         =   R12_dXC::num_points;
    r12_dx->num_points_cum     =   R12_dXC::num_points_cum;
    r12_dx->nx                 =   R12_dXC::ncx;
    r12_dx->ny                 =   R12_dXC::ncy;
    r12_dx->x_log_scale        =   R12_dXC::x_log_scale;
    r12_dx->x_map_scale        =   R12_dXC::x_map_scale;
    r12_dx->x_map_scale_log    =   R12_dXC::x_map_scale_log;
    r12_dx->x_max              =   R12_dXC::x_max;
    r12_dx->x_min              =   R12_dXC::x_min;
    r12_dx->x_min_l10          =   R12_dXC::x_min_l10;
    r12_dx->y_log_scale        =   R12_dXC::y_log_scale;
    r12_dx->y_map_scale        =   R12_dXC::y_map_scale;
    r12_dx->y_map_scale_log    =   R12_dXC::y_map_scale_log;
    r12_dx->y_max              =   R12_dXC::y_max;
    r12_dx->y_min              =   R12_dXC::y_min;
    r12_dx->y_min_l10          =   R12_dXC::y_min_l10;
    r12_dx->is_build           =   true;

    // Intialize object
    r12_dx->initialize();

}


void set_data_r21(R21* r21)
{
    // Load data into target object
    r21->c                  =   R21C::c;
    r21->intervals_bounds_x =   R21C::interval_bounds_x;
    r21->intervals_bounds_y =   R21C::interval_bounds_y;
    r21->num_c              =   R21C::num_c;
    r21->num_intervals_x    =   R21C::num_intervals_x;
    r21->num_intervals_y    =   R21C::num_intervals_y;
    r21->num_points         =   R21C::num_points;
    r21->num_points_cum     =   R21C::num_points_cum;
    r21->nx                 =   R21C::ncx;
    r21->ny                 =   R21C::ncy;
    r21->x_log_scale        =   R21C::x_log_scale;
    r21->x_map_scale        =   R21C::x_map_scale;
    r21->x_map_scale_log    =   R21C::x_map_scale_log;
    r21->x_max              =   R21C::x_max;
    r21->x_min              =   R21C::x_min;
    r21->x_min_l10          =   R21C::x_min_l10;
    r21->y_log_scale        =   R21C::y_log_scale;
    r21->y_map_scale        =   R21C::y_map_scale;
    r21->y_map_scale_log    =   R21C::y_map_scale_log;
    r21->y_max              =   R21C::y_max;
    r21->y_min              =   R21C::y_min;
    r21->y_min_l10          =   R21C::y_min_l10;
    r21->is_build           =   true;

    // Intialize object
    r21->initialize();

}


void set_data_r21_dx(R21_dX* r21_dx)
{
    // Load data into target object
    r21_dx->c                  =   R21_dXC::c;
    r21_dx->intervals_bounds_x =   R21_dXC::interval_bounds_x;
    r21_dx->intervals_bounds_y =   R21_dXC::interval_bounds_y;
    r21_dx->num_c              =   R21_dXC::num_c;
    r21_dx->num_intervals_x    =   R21_dXC::num_intervals_x;
    r21_dx->num_intervals_y    =   R21_dXC::num_intervals_y;
    r21_dx->num_points         =   R21_dXC::num_points;
    r21_dx->num_points_cum     =   R21_dXC::num_points_cum;
    r21_dx->nx                 =   R21_dXC::ncx;
    r21_dx->ny                 =   R21_dXC::ncy;
    r21_dx->x_log_scale        =   R21_dXC::x_log_scale;
    r21_dx->x_map_scale        =   R21_dXC::x_map_scale;
    r21_dx->x_map_scale_log    =   R21_dXC::x_map_scale_log;
    r21_dx->x_max              =   R21_dXC::x_max;
    r21_dx->x_min              =   R21_dXC::x_min;
    r21_dx->x_min_l10          =   R21_dXC::x_min_l10;
    r21_dx->y_log_scale        =   R21_dXC::y_log_scale;
    r21_dx->y_map_scale        =   R21_dXC::y_map_scale;
    r21_dx->y_map_scale_log    =   R21_dXC::y_map_scale_log;
    r21_dx->y_max              =   R21_dXC::y_max;
    r21_dx->y_min              =   R21_dXC::y_min;
    r21_dx->y_min_l10          =   R21_dXC::y_min_l10;
    r21_dx->is_build           =   true;

    // Intialize object
    r21_dx->initialize();

}


void set_data_r22(R22* r22)
{
    // Load data into target object
    r22->c                  =   R22C::c;
    r22->intervals_bounds_x =   R22C::interval_bounds_x;
    r22->intervals_bounds_y =   R22C::interval_bounds_y;
    r22->num_c              =   R22C::num_c;
    r22->num_intervals_x    =   R22C::num_intervals_x;
    r22->num_intervals_y    =   R22C::num_intervals_y;
    r22->num_points         =   R22C::num_points;
    r22->num_points_cum     =   R22C::num_points_cum;
    r22->nx                 =   R22C::ncx;
    r22->ny                 =   R22C::ncy;
    r22->x_log_scale        =   R22C::x_log_scale;
    r22->x_map_scale        =   R22C::x_map_scale;
    r22->x_map_scale_log    =   R22C::x_map_scale_log;
    r22->x_max              =   R22C::x_max;
    r22->x_min              =   R22C::x_min;
    r22->x_min_l10          =   R22C::x_min_l10;
    r22->y_log_scale        =   R22C::y_log_scale;
    r22->y_map_scale        =   R22C::y_map_scale;
    r22->y_map_scale_log    =   R22C::y_map_scale_log;
    r22->y_max              =   R22C::y_max;
    r22->y_min              =   R22C::y_min;
    r22->y_min_l10          =   R22C::y_min_l10;
    r22->is_build           =   true;

    // Intialize object
    r22->initialize();

}


void set_data_r22_dx(R22_dX* r22_dx)
{
    // Load data into target object
    r22_dx->c                  =   R22_dXC::c;
    r22_dx->intervals_bounds_x =   R22_dXC::interval_bounds_x;
    r22_dx->intervals_bounds_y =   R22_dXC::interval_bounds_y;
    r22_dx->num_c              =   R22_dXC::num_c;
    r22_dx->num_intervals_x    =   R22_dXC::num_intervals_x;
    r22_dx->num_intervals_y    =   R22_dXC::num_intervals_y;
    r22_dx->num_points         =   R22_dXC::num_points;
    r22_dx->num_points_cum     =   R22_dXC::num_points_cum;
    r22_dx->nx                 =   R22_dXC::ncx;
    r22_dx->ny                 =   R22_dXC::ncy;
    r22_dx->x_log_scale        =   R22_dXC::x_log_scale;
    r22_dx->x_map_scale        =   R22_dXC::x_map_scale;
    r22_dx->x_map_scale_log    =   R22_dXC::x_map_scale_log;
    r22_dx->x_max              =   R22_dXC::x_max;
    r22_dx->x_min              =   R22_dXC::x_min;
    r22_dx->x_min_l10          =   R22_dXC::x_min_l10;
    r22_dx->y_log_scale        =   R22_dXC::y_log_scale;
    r22_dx->y_map_scale        =   R22_dXC::y_map_scale;
    r22_dx->y_map_scale_log    =   R22_dXC::y_map_scale_log;
    r22_dx->y_max              =   R22_dXC::y_max;
    r22_dx->y_min              =   R22_dXC::y_min;
    r22_dx->y_min_l10          =   R22_dXC::y_min_l10;
    r22_dx->is_build           =   true;

    // Intialize object
    r22_dx->initialize();

}