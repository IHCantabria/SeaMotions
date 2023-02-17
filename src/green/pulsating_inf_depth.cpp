
// Include general usage libraries
#include <tuple>

// Include local modules
#include "../config.hpp"
#include "integrals_db.hpp"
#include "../math/math_tools.hpp"
#include "pulsating_inf_depth.hpp"


cusfloat calculate_dr_dx(cusfloat R, cusfloat dx)
{
    /**
     * @brief Calcualte derivative of R with respect to X
     * 
     * dr_dx = -dx/R
     * where R:sqrt((x-xi)**2.0+(y-yi)**2.0)
     * 
     * \param R Euclidean distance between the field and the source points in the horizontal plane.
     * \param dx X distance in between the field and the source points in the horizontal plane
     */

    return -dx/R;
}


cusfloat calculate_r(
                    cusfloat x,
                    cusfloat y,
                    cusfloat xi,
                    cusfloat eta
                    )
{
    /**
     * @brief Calculate euclidean distance in between field and source point in the
     *        horizontal plane
     * 
     * \param x X Coordinate of the field point
     * \param y Y Coordinate of the field point
     * \param xi X Coordinate of the source point
     * \param eta Y Coordinate of the source point
    * 
    */

    return sqrt( pow2s(x-xi) + pow2s(y-eta) );
}


cusfloat wave_term_inf_depth(
                            cusfloat x_ndim,
                            cusfloat y_ndim,
                            IntegralsDb &idb
                            )
{
    cusfloat sol = 0.0;
    if (x_ndim < 3.0 && y_ndim < 4.0)
    {
        sol = idb.r11->calculate_xy(x_ndim, y_ndim);
    }
    else if (x_ndim < 3.0)
    {
        sol = idb.r12->calculate_xy(x_ndim, y_ndim);
    }
    else if (y_ndim < 4.0)
    {
        sol = idb.r21->calculate_xy(x_ndim, y_ndim);
    }
    else
    {
        sol = idb.r22->calculate_xy(x_ndim, y_ndim);
    }

    return sol;
}


std::tuple<cusfloat, cusfloat> wave_term_inf_depth_dhoriz(
                                                            cusfloat x,
                                                            cusfloat y,
                                                            cusfloat z,
                                                            cusfloat xi,
                                                            cusfloat eta,
                                                            cusfloat zeta,
                                                            cusfloat nu,
                                                            IntegralsDb &idb
                                                            )
{
    // Calculate derivative properties
    cusfloat x_ndim = nu*calculate_r(x, y, xi, eta);
    cusfloat y_ndim = nu*abs(z+zeta);

    // Calculate dG_dX
    cusfloat dg_dxndim = wave_term_inf_depth_dxndim(x_ndim, y_ndim, idb);

    // Calculate dX_dx
    cusfloat dxndim_dx = nu*calculate_dr_dx(x_ndim/nu, x-xi);

    // Calculate dX_dy
    cusfloat dxndim_dy = nu*calculate_dr_dx(x_ndim/nu, y-eta);

    return std::make_tuple(
                            dg_dxndim*dxndim_dx,
                            dg_dxndim*dxndim_dy
                            );
}


cusfloat wave_term_inf_depth_dx(
                                cusfloat x,
                                cusfloat y,
                                cusfloat z,
                                cusfloat xi,
                                cusfloat eta,
                                cusfloat zeta,
                                cusfloat nu,
                                IntegralsDb &idb
                                )
{
    // Calculate derivative properties
    cusfloat x_ndim = nu*calculate_r(x, y, xi, eta);
    cusfloat y_ndim = nu*abs(z+zeta);

    // Calculate dG_dX
    cusfloat dg_dxndim = wave_term_inf_depth_dxndim(x_ndim, y_ndim, idb);

    // Calculate dX_dx
    cusfloat dxndim_dx = nu*calculate_dr_dx(x_ndim/nu, x-xi);

    return dg_dxndim*dxndim_dx;
}


cusfloat wave_term_inf_depth_dxndim(
                                cusfloat x_ndim,
                                cusfloat y_ndim,
                                IntegralsDb &idb
                                )
{
    // Get integral value depending on the region
    cusfloat sol = 0.0;
    if (x_ndim < 3.0 && y_ndim < 4.0)
    {
        if (x_ndim < 1.0)
        {
            sol = idb.r11a_dx->calculate_xy(x_ndim, y_ndim);
        }
        else
        {
            sol = idb.r11b_dx->calculate_xy(x_ndim, y_ndim);
        }
    }
    else if (x_ndim < 3.0)
    {
        sol = idb.r12_dx->calculate_xy(x_ndim, y_ndim);
    }
    else if (y_ndim < 4.0)
    {
        sol = idb.r21_dx->calculate_xy(x_ndim, y_ndim);
    }
    else
    {
        sol = idb.r22_dx->calculate_xy(x_ndim, y_ndim);
    }

    return sol;
}


cusfloat wave_term_inf_depth_dy(
                                cusfloat x,
                                cusfloat y,
                                cusfloat z,
                                cusfloat xi,
                                cusfloat eta,
                                cusfloat zeta,
                                cusfloat nu,
                                IntegralsDb &idb
                                )
{
    // Calculate derivative properties
    cusfloat x_ndim = nu*calculate_r(x, y, xi, eta);
    cusfloat y_ndim = nu*abs(z+zeta);

    // Calculate dG_dX
    cusfloat dg_dxndim = wave_term_inf_depth_dxndim(x_ndim, y_ndim, idb);

    // Calculate dX_dy
    cusfloat dxndim_dy = nu*calculate_dr_dx(x_ndim/nu, y-eta);

    return dg_dxndim*dxndim_dy;
}



cusfloat wave_term_inf_depth_dyndim(
                                    cusfloat x_ndim,
                                    cusfloat y_ndim,
                                    IntegralsDb &idb
                                    )
{
    return  (
                -
                2/sqrt(pow2s(x_ndim) + pow2s(y_ndim))
                -
                wave_term_inf_depth(x_ndim, y_ndim, idb)
            );

}


cusfloat wave_term_inf_depth_dz(
                                cusfloat x_ndim,
                                cusfloat z,
                                cusfloat zeta,
                                cusfloat nu,
                                IntegralsDb &idb
                                )
{
    // Calculate derivative properties
    cusfloat y_ndim = nu*abs(z+zeta);

    // Calculate derivative with respet to Y non dimensional variable
    cusfloat dG_dY = wave_term_inf_depth_dyndim(x_ndim, y_ndim, idb);

    // Calculate derivative with respecto to y variable
    cusfloat dY_dz = nu*sign(z+zeta);

    return dG_dY*dY_dz;
}


cusfloat wave_term_inf_depth_dz(
                                cusfloat x,
                                cusfloat y,
                                cusfloat z,
                                cusfloat xi,
                                cusfloat eta,
                                cusfloat zeta,
                                cusfloat nu,
                                IntegralsDb &idb
                                )
{
    // Calculate derivative properties
    cusfloat x_ndim = nu*calculate_r(x, y, xi, eta);

    return wave_term_inf_depth_dz(x_ndim, z, zeta, nu, idb);
}