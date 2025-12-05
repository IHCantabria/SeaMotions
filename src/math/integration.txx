
/*
 * Copyright (c) 2025 Sergio Fernández Ruano / IHCantabria
 *
 * This file is part of SeaMotions Software.
 *
 * SeaMotions is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SeaMotions is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

// Include local modules
#include "../config.hpp"
#include "../containers/panel_geom_list.hpp"
#include "gauss.hpp"
#include "gauss_t.hpp"
#include "topology.hpp"
#include "../mesh/tools.hpp"

// Include namespaces
using namespace std::literals::complex_literals;


template<typename T, typename U>   
inline cuscomplex   _adaptive_quadrature_panel(
                                                    T*              panel,
                                                    U               target_fcn,
                                                    cuscomplex      prev_int,
                                                    cusfloat        abs_tol,
                                                    cusfloat        rel_tol,
                                                    int             prev_gpo,
                                                    int             adapt_level
                                            )
{
    // Get current gauss points
    int gpo = prev_gpo + 2;
    GaussPoints gp( gpo );

    // Re-mesh current panel
    PanelGeomList* panel_list;
    refine_element( 
                        panel,
                        panel_list
                    );

    // Loop over child panels to integrate the target function
    // over them
    cuscomplex    int_sol     = 0.0 + 0.0i;
    cuscomplex*   int_values  = generate_empty_vector<cuscomplex>( panel_list->panel_np );
    for ( int i=0; i<panel_list->panel_np; i++ )
    {
        int_values[i]   = quadrature_panel(
                                              panel_list->panels[i],
                                              target_fcn,
                                              &gp
                                          );
        int_sol         += int_values[i];
    }

    // Compare the cumulative integral value with the previous
    // solution
    bool is_equal_abs = assert_complex_equality( 
                                                    prev_int, 
                                                    int_sol,
                                                    abs_tol 
                                                );
    bool is_equal_rel = false;
    if ( int_sol != 0.0 )
    {
        assert_complex_equality(
                                    ( prev_int - int_sol ) / int_sol,
                                    cuscomplex( 0.0, 0.0 ),
                                    rel_tol
                                );
    }

    bool is_equal = !( is_equal_abs || is_equal_rel );

    if ( 
            is_equal
            &&
            ( adapt_level < 6 )
        )
    {
        // Increase adaption level
        adapt_level += 1;

        // Loop over child panels
        int_sol = 0.0;
        for ( int i=0; i<panel_list->panel_np; i++ )
        {
            int_sol +=  _adaptive_quadrature_panel(
                                                        panel_list->panels[i],
                                                        target_fcn,
                                                        int_values[i],
                                                        abs_tol,
                                                        rel_tol,
                                                        gpo,
                                                        adapt_level
                                                    );
        }
    }
    else if ( is_equal )
    {
        std::cerr << std::endl;
        std::cerr << "ERROR: Adaptive quadrature could not find the solution ";
        std::cerr << "with the required accuracy. Maximum adaption levels reached";
        std::cerr << std::endl;
        throw std::runtime_error( "" );
    }

    // Delete local heap memory
    mkl_free( int_values );

    return int_sol;
}


template<typename T, typename U>   
inline cuscomplex adaptive_quadrature_panel(
                                                    T*              panel,
                                                    U               target_fcn,
                                                    cusfloat        abs_tol,
                                                    cusfloat        rel_tol,
                                                    bool            block_adaption,
                                                    bool            even_order,
                                                    int             go_fixed
                                            )
{
    // Define local variables
    bool is_equal_abs   = false;
    bool is_equal_rel   = false;

    // Integrate parent panel
    cuscomplex  int_sol_0, int_sol_00, int_sol_1;
    
    if ( !block_adaption )
    {
        // Define polynomial adaption parameters
        int first_order     = 2;
        int last_order      = 10;
        int step_order      = 2;
        if ( !even_order )
        {
            first_order     = 1;
            step_order      = 1;
        }
        int_sol_0   =  quadrature_panel(
                                        panel,
                                        target_fcn,
                                        first_order
                                    );
        int_sol_00  = int_sol_0;
        for ( int igo=first_order+step_order; igo<last_order; igo+=step_order )
        {
            // Integrate function with the new
            int_sol_1   =  quadrature_panel(
                                                panel,
                                                target_fcn,
                                                igo
                                            );

            // Check for convergence
            is_equal_abs = assert_complex_equality( 
                                                        int_sol_0, 
                                                        int_sol_1, 
                                                        abs_tol 
                                                    );
            if ( int_sol_1 != 0.0 )
            {
                is_equal_rel = assert_complex_equality( 
                                                            ( int_sol_0 - int_sol_1 ) / int_sol_1, 
                                                            cuscomplex( 0.0, 0.0 ), 
                                                            rel_tol 
                                                        );
            }

            if ( is_equal_abs || is_equal_rel )
            {
                break;
            }

            // Save last integration value to compare with the
            // next integration value
            int_sol_0 = int_sol_1;
        }

    }
    else
    {
        int_sol_1   =  quadrature_panel(
                                            panel,
                                            target_fcn,
                                            go_fixed
                                        );
    }

    // Define adaption level
    int adapt_level = 0;
                                        
    // Launch adaptive interation
    if ( 
            !block_adaption 
            && 
            !( is_equal_abs || is_equal_rel ) 
        )
    {
        int_sol_1   =   _adaptive_quadrature_panel(
                                                        panel,
                                                        target_fcn,
                                                        int_sol_1,
                                                        abs_tol,
                                                        rel_tol,
                                                        4,
                                                        adapt_level
                                                    );
    }

    return int_sol_1;
}


template<int NGP, typename T, typename F, typename P>
inline void gauss2d_loop(T& storage, F&& f, const P* panel)
{
    for (int i = 0; i < NGP*NGP; i++)
    {
        storage += GaussPointsT<2,NGP>::weights_x[i] *
                   GaussPointsT<2,NGP>::weights_y[i] *
                   f(i) *
                   panel->jac_det_gauss_points[i];
    }
}


template<typename T, typename U>
cuscomplex  quadrature_panel(
                                T*          panel,
                                U           target_fcn,
                                int         gp_order
                            )
{
    // Define Gauss points
    GaussPoints* gp = new GaussPoints( gp_order );

    // Launch quadrature calculation
    cuscomplex  int_value   = quadrature_panel(
                                                    panel,
                                                    target_fcn,
                                                    gp
                                                );

    // Delete Gauss Points
    delete gp;

    return int_value;
}


template<typename T, typename U>
cuscomplex  quadrature_panel(
                                T*              panel,
                                U               target_fcn,
                                GaussPoints*    gp
                            )
{
    // Define variable to hold the integration value
    cuscomplex int_value( 0.0, 0.0 );

    // Loop over gauss points to perform the integration
    cuscomplex  fcn_val          = 0.0;
    cusfloat    gp_global[3]     = { 0.0, 0.0, 0.0 };
    for ( int i=0; i<gp->np; i++ )
    {
        for ( int j=0; j<gp->np; j++ )
        {
            // Get global coordinates for the gauss points
            panel->local_to_global( 
                                        gp->roots[i],
                                        gp->roots[j],
                                        gp_global
                                    );

            
            // Calculate target function value
            fcn_val = target_fcn( 
                                    gp->roots[i],
                                    gp->roots[j],
                                    gp_global[0], 
                                    gp_global[1], 
                                    gp_global[2] 
                                );

            // Calculate function integral function value
            int_value += gp->weights[i]*gp->weights[j]*fcn_val*jacobi_det_2d( 
                                                                                panel->num_nodes,
                                                                                panel->xl,
                                                                                panel->yl,
                                                                                gp->roots[i],
                                                                                gp->roots[j]
                                                                            );
        }
    }

    return int_value;
}


template<typename T, typename U, auto Kernel, int NGP>
void        quadrature_panel_t(
                                        T*                  panel,
                                        U                   target_fcn,
                                        cuscomplex&         result_G,
                                        cuscomplex&         result_G_dn_sf,
                                        cuscomplex&         result_G_dn_pf,
                                        bool                verbose
                                    )
{
    // Calculate function value
    target_fcn.template operator()<Kernel>( 
                                                panel->xl, 
                                                panel->yl, 
                                                panel->gauss_points_global_x, 
                                                panel->gauss_points_global_y, 
                                                panel->gauss_points_global_z,
                                                verbose
                                            );

    result_G        = 0.0;
    result_G_dn_sf  = 0.0;
    result_G_dn_pf  = 0.0;

    gauss2d_loop<NGP>(
                            result_G,
                            [&](int i){ return target_fcn.G[i]; },
                            panel
                        );

    gauss2d_loop<NGP>(
                            result_G_dn_sf,
                            [&](int i){ return target_fcn.dG_dn_sf[i]; },
                            panel
                        );

    gauss2d_loop<NGP>(
                            result_G_dn_pf,
                            [&](int i){ return target_fcn.dG_dn_pf[i]; },
                            panel
                        );

}


template<typename Functor>
cusfloat romberg_quadrature(
                                Functor f, 
                                cusfloat a, 
                                cusfloat b, 
                                cusfloat precision
                            )
{
    // Define buffers
    const int max_steps = 50;
    cusfloat* rp = new cusfloat [max_steps];
    cusfloat* rc = new cusfloat [max_steps];

    // Calculate maximum step size
    cusfloat h = b - a;

    // Calculate first trapezoidal point
    rp[0] = (f(a)+f(b))*h/2.0;

    // Loop until the maximum step limit
    cusfloat* aux_ptr = nullptr;
    cusfloat local_trap = 0.0;
    int trap_steps = 0;
    for (int i=1; i<max_steps; i++)
    {
        // Calculate i trapezoidal rule value
        local_trap = 0.0;
        trap_steps = 1 << (i-1);
        h /= 2.0;
        for (int j=1; j<=trap_steps; j++)
        {
            local_trap += f(a+(2*j-1)*h);
        }
        local_trap *= h;

        // Calculate current buffer values
        rc[0] = local_trap + rp[0]/2.0;
        for (int j=1; j<=i; j++)
        {
            rc[j] = rc[j-1] + 1/(pow(4, j)-1)*(rc[j-1]-rp[j-1]);
        }

        // Check for convergence
        if (std::abs(rc[i]-rp[i-1])<precision)
        {
            return rc[i];
        }

        // Change buffers to have the current as the previous
        aux_ptr = rc;
        rc = rp;
        rp = aux_ptr;

    }

    // Get value to return
    cusfloat last_value = rp[max_steps-1];

    // Delete heap memory
    delete [] rp;
    delete [] rc;

    // Return the last value
    std::cout << "WARNING: Romberg quadrature could not find the integral ";
    std::cout << "value with the requested accuracy." << std::endl;
    return last_value;
}