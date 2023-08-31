
// Include local modules
#include "../config.hpp"
#include "../containers/panel_geom_list.hpp"
#include "gauss.hpp"
#include "topology.hpp"
#include "../mesh/tools.hpp"

// Include namespaces
using namespace std::literals::complex_literals;


template<typename T>   
inline cuscomplex   _adaptive_quadrature_panel(
                                                    PanelGeom*      panel,
                                                    T               target_fcn,
                                                    cuscomplex      prev_int,
                                                    cusfloat        tol,
                                                    GaussPoints*    gp,
                                                    int             adapt_level
                                            )
{
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
                                              gp
                                          );
        int_sol         += int_values[i];
    }

    // Compare the cumulative integral value with the previous
    // solution
    bool is_equal = !assert_complex_equality( prev_int, int_sol, tol );
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
                                                        tol,
                                                        gp,
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
        std::runtime_error( "" );
    }

    // Delete local heap memory
    mkl_free( int_values );

    return int_sol;
}


template<typename T>   
inline cuscomplex   adaptive_quadrature_panel(
                                                    PanelGeom*  panel,
                                                    T           target_fcn,
                                                    cusfloat    tol,
                                                    int         gp_order
                                            )
{
    // Start gauss points
    GaussPoints* gp = new GaussPoints( gp_order );

    // Launch adaptive interation
    cuscomplex  int_sol =  adaptive_quadrature_panel(
                                                        panel,
                                                        target_fcn,
                                                        tol,
                                                        gp
                                                    );
    
    // Delete local heap memory
    delete gp;

    return int_sol;
}


template<typename T>   
inline cuscomplex adaptive_quadrature_panel(
                                                    PanelGeom*      panel,
                                                    T               target_fcn,
                                                    cusfloat        tol,
                                                    GaussPoints*    gp
                                            )
{
    // Integrate parent panel
    cuscomplex int_sol  =  quadrature_panel(
                                                panel,
                                                target_fcn,
                                                gp
                                            );

    // Define adaption level
    int adapt_level = 0;
                                        
    // Launch adaptive interation
    int_sol     =   _adaptive_quadrature_panel(
                                                    panel,
                                                    target_fcn,
                                                    int_sol,
                                                    tol,
                                                    gp,
                                                    adapt_level
                                                );

    return int_sol;
}


template<typename T>
cuscomplex  quadrature_panel(
                                PanelGeom*  panel,
                                T           target_fcn,
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

    return int_value;
}


template<typename T>
cuscomplex  quadrature_panel(
                                PanelGeom*      panel,
                                T               target_fcn,
                                GaussPoints*    gp
                            )
{
    // Define variable to hold the integration value
    cuscomplex int_value        = 0.0;

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
            fcn_val = target_fcn( gp_global[0], gp_global[1], gp_global[2] );

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
            // std::cout << " --> Root: " << a+(2*j-1)*h << " - Value: " << f(a+(2*j-1)*h) << std::endl;
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