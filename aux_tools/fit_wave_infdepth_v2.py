
# Import general usage libraries
import copy
from typing import Callable

# Import general usage scientific libraries
import matplotlib.pyplot as plt
import numpy as np
from numpy import array, exp, linspace, log, log10, pi, polyval, ndarray, sin, sqrt, zeros
from numpy import abs as np_abs
from scipy.special import eval_legendre, jv, legendre, roots_chebyt, yv, erf

# Import local modules
from base_integrals_v2 import fxy, fxy_polar, fxy_dx, fxy_dx_polar
from fit_cheby_v2 import FitProperties, fit_integral_1d, fit_integral_2d, RefLevel, fit_residual_2D_adaptive_interface
from fit_functions import FcnParameters, step_function_1d, exp_squared, sin_wave, sin_wave_exp, besselj0


FIT_TOL = 1E-6


def factorial(n: int)->int:
    if n < 2:
        return 1
    else:
        f = 1
        for i in range(2, n+1):
            f *= i
        
        return f
    

def fcn_1d( x ):
    return np.log10( 10**x )


def fcn_2d( x, y ):
    xi = np.log10( x )
    yi = np.log10( y )
    return np.log10( 10**xi * 10**yi )
    

def fit_residual_region_1d_test( )->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.x_log_scale = False
    fit_props.x_max = -3
    fit_props.x_min = 3
    fit_props.y_log_scale = False
    fit_props.y_max = 4.0
    fit_props.y_min = 1e-3
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_tol = FIT_TOL
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_1d(fcn_1d, fit_props, "region_test", show_figs=True)


def fit_residual_region_2d_test( )->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.x_log_scale = True
    fit_props.x_max =  1e-3
    fit_props.x_min = 1e3
    fit_props.y_log_scale = True
    fit_props.y_max = 1e-3
    fit_props.y_min = 1e3
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_tol = FIT_TOL
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(fcn_2d, fit_props, "region_test", stop_to_show=True)


def fit_residual_region_test()->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.x_log_scale = False
    fit_props.x_max = 1e-12
    fit_props.x_min = 1e-3
    fit_props.y_log_scale = False
    fit_props.y_max = 4.0
    fit_props.y_min = 1e-3
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.cheby_tol = FIT_TOL
    fit_props.x_map_fcn = lambda x: x
    fit_props.y_map_fcn = lambda y: y

    # Launch fit
    fit_integral_2d(residual_region_11, fit_props, "region_test", show_figs=True)


def fit_residual_region_11(show_figs=False, show_summary_fig=False):
    # Define boundary values
    # x = 10**array( [ np.log10( 1.0 ), np.log10( 3.16e2 ) ] )
    # y = 10**array( [ np.log10( 1.0 ), np.log10( 3.16e2 ) ] )
    # y = 10**array( [ -5.0, np.log10( 1.5 ) ] )
    # y = 10**array( [ -5.0, np.log10( 0.005 ) ] )
    # y = 10**array( [ -5.0, 3.4 ] )
    # x = array( [ 1e-12, 1e-3, 3 ] )
    # y = array( [ 1e-12, 1e-3, 4 ] )
    # x = array( [ 1e-8, 1e-3, 1e-2, 1e-1, 5e-1, 1.0, 1.5, 2.0, 2.5, 3 ] ) # There is room to lower the number of coefficients by using 0.25 spacing instead of 0.5
    # y = array( [ 1e-8, 1e-3, 1e-2, 1e-1, 5e-1, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4 ] )
    # x = np.linspace( 2.48, 2.5, 2 ) # There is room to lower the number of coefficients by using 0.25 spacing instead of 0.5
    # y = np.linspace( 4.10125, 4.12125, 2 )  
    x = np.linspace( 2.46, 2.5, 2 ) # There is room to lower the number of coefficients by using 0.25 spacing instead of 0.5
    y = np.linspace( 4.10125, 4.14, 2 ) 
    # x = 10**np.linspace( -5, 3.4, 10 ) # There is room to lower the number of coefficients by using 0.25 spacing instead of 0.5
    # y = 10**np.linspace( -5, np.log10( 1.2 ), 10 )
    # y = 10**np.linspace( -5, np.log10( np.pi/2.0-1e-5 ), 20 )
    # x = array( [ 1e-8, 1e-3 ] )
    # y = array( [ 1e-2, 1e-1 ] )
    # x = array([1e-12, 5e-12, 1e-11, 5e-11, 1e-10, 5e-10, 1e-9, 5e-9, 1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1, 1, 2.0, 3.0])
    # y = array([1e-12, 5e-12, 1e-11, 5e-11, 1e-10, 5e-10, 1e-9, 5e-9, 1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1, 1, 2.0, 4.0])
    # x = 10**linspace(-8.0, log10(3.0), 100)
    # y = 10**linspace(-5.0, log10(4.0), 100)
    results = []
    for i in range(x.shape[0]-1):
        for j in range(y.shape[0]-1):
            region_name = f"fit_residual_region_11 -> X: {x[i]:0.3E} - {x[i+1]:0.3E} | Y: {y[j]:0.3E} - {y[j+1]:0.3E}"

            # Set up domain extension and fit properties
            fit_props = FitProperties()
            fit_props.region_name = "residual_region_11"
            fit_props.cheby_order_x = 20
            fit_props.cheby_order_y = 20
            fit_props.fcn_log_scale = True
            fit_props.x_log_scale = True
            fit_props.x_max = x[i+1]
            fit_props.x_min = x[i]
            fit_props.y_log_scale = True
            fit_props.y_max = y[j+1]
            fit_props.y_min = y[j]
            fit_props.cheby_abs_tol = 1e-6
            fit_props.cheby_rel_tol = 1e-14
            fit_props.x_map_fcn = lambda x: x
            fit_props.y_map_fcn = lambda y: y

            # Launch fit
            results.append(fit_integral_2d(residual_region_11, fit_props, show_figs=show_figs))

            # raise ValueError( "Stop by user" )

    # Plot results summay
    plot_results_summary(x, y, results, "region_11", log_scale=True, show_figs=show_summary_fig)

    return x, y, results


def fit_residual_region_11_adaptive( folder_path: str, show_figs=False, is_square_ref=True, show_summary_fig=False ) -> None:
    # Define boundary values
    fit_props               = FitProperties( )
    fit_props.dims          = 2
    fit_props.region_name   = "region_11"
    fit_props.cheby_order_x = 5
    fit_props.cheby_order_y = 5
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = True
    fit_props.x_max         = 5.0
    fit_props.x_min         = -5.0
    fit_props.y_log_scale   = True
    fit_props.y_max         = 5.0
    fit_props.y_min         = -5.0
    fit_props.cheby_abs_tol = 1E-6
    fit_props.cheby_rel_tol = 1E-14
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.max_ref_level = 10
    
    fit_props.generate_fitting_matrix( )

    fit_function            = residual_region_11
    


    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_2d( fit_function, fit_props, show_figs=show_figs ) )

    # Start adaptive refinement loop
    fit_residual_2D_adaptive_interface(fit_function, ref_level, is_square_ref=is_square_ref, show_figs=show_figs )

    # Set starting position
    ref_level.set_start_index( 0 )

    # Plot results summary
    # ref_level.show_summary( folder_path )

    return ref_level


def fit_residual_region_11_dx_adaptive( folder_path: str, show_figs=False, is_square_ref=True, show_summary_fig=False ) -> None:
    # Define boundary values
    fit_props               = FitProperties( )
    fit_props.dims          = 2
    fit_props.region_name   = "region_11_dx"
    fit_props.cheby_order_x = 7
    fit_props.cheby_order_y = 7
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = True
    fit_props.x_max         = 5.0
    fit_props.x_min         = -5.0
    fit_props.y_log_scale   = True
    fit_props.y_max         = 5.0
    fit_props.y_min         = -5.0
    fit_props.cheby_abs_tol = 1E-6
    fit_props.cheby_rel_tol = 1E-14
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.max_ref_level = 10

    fit_props.generate_fitting_matrix( )

    fit_function            = residual_region_11_dx

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_2d( fit_function, fit_props, show_figs=show_figs ) )

    # Start adaptive refinement loop
    fit_residual_2D_adaptive_interface( fit_function, ref_level, is_square_ref=is_square_ref, show_figs=show_figs )

    # Set starting position
    ref_level.set_start_index( 0 )

    # Plot results summary
    # ref_level.show_summary( folder_path )

    return ref_level


def fit_residual_region_AA_dx_adaptive( folder_path: str, show_figs=False, is_square_ref=True, show_summary_fig=False ) -> None:
    # Define boundary values
    fit_props               = FitProperties( )
    fit_props.region_name   = "region_11"
    fit_props.cheby_order_x = 20
    fit_props.cheby_order_y = 20
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = False
    fit_props.x_max         = 30.0
    fit_props.x_min         = 0.0
    fit_props.y_log_scale   = False
    fit_props.y_max         = 30.0
    fit_props.y_min         = 0.0
    fit_props.cheby_abs_tol = 1E-6
    fit_props.cheby_rel_tol = 1E-14
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.max_ref_level = 8

    fcn_params              = FcnParameters( )
    # fit_function            = lambda x, y: step_function_1d( fcn_params, x )
    # fit_function            = lambda x, y: exp_squared( x, y )
    fit_function            = lambda x, y: sin_wave( x )
    # fit_function            = lambda x, y: sin_wave_exp( x, y )
    # fit_function            = lambda x, y: besselj0( x )

    # Create root FitRegion
    ref_level               = RefLevel( copy.deepcopy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_2d( fit_function , fit_props, show_figs=show_figs ) )

    # Start adaptive refinement loop
    fit_residual_2D_adaptive_interface( fit_function, ref_level, is_square_ref=is_square_ref, show_figs=show_figs )

    # Set starting position
    ref_level.set_start_index( 0 )

    # Plot results summary
    ref_level.show_summary( folder_path )

    return ref_level


# def fit_residual_2D_adaptive_interface(fit_function, ref_level: RefLevel, is_square_ref=True, show_figs=False) -> None:
#     # Define boundary values
#     x_new = np.array( [ 
#                             ref_level.fit_props.x_min, 
#                             ( ref_level.fit_props.x_min + ref_level.fit_props.x_max ) / 2.0,
#                             ref_level.fit_props.x_max
#                         ] )
#     y_new = np.array( [ 
#                             ref_level.fit_props.y_min, 
#                             ( ref_level.fit_props.y_min + ref_level.fit_props.y_max ) / 2.0,
#                             ref_level.fit_props.y_max
#                         ] )
    
#     if is_square_ref:
#         # Loop over refinements to create new refinement levels
#         for i in range( x_new.shape[0] - 1 ):
#             for j in range( y_new.shape[0] - 1 ):
#                 # Apply new sub-region boundaries
#                 fit_props       = copy.deepcopy( ref_level.fit_props )
#                 fit_props.x_max = x_new[i+1]
#                 fit_props.x_min = x_new[i]
#                 fit_props.y_max = y_new[j+1]
#                 fit_props.y_min = y_new[j]

#                 # Create new refinement level
#                 ref_level_i                 = RefLevel( fit_props, parent=ref_level )

#                 # Fit function over the interval
#                 ref_level_i.add_data( *fit_integral_2d( fit_function, fit_props, show_figs=show_figs ) )

#                 print( f" -> Ref.Level:     ", ref_level_i.level, flush=True )
#                 print( f" -> Passed?:       ", ref_level_i.check_tolerances( ), flush=True )
#                 print( "\n" )

#                 # Check fit tolerances and refine if any
#                 if not ref_level_i.check_tolerances( ) and ref_level_i.level < ref_level_i.fit_props.max_ref_level:
#                     fit_residual_2D_adaptive_interface(fit_function, ref_level_i, show_figs=show_figs )
                
#                 # Add ith refinement level region to the children list
#                 ref_level.child.append( ref_level_i )

#     else:
#         # Loop over refinements to create new refinement levels
#         ref_levels_x = [ ]
#         for i in range( x_new.shape[0] - 1 ):
#             # Apply new sub-region boundaries
#             fit_props       = copy.deepcopy( ref_level.fit_props )
#             fit_props.x_max = x_new[i+1]
#             fit_props.x_min = x_new[i]

#             # Create new refinement level
#             ref_level_i     = RefLevel( fit_props, parent=ref_level )

#             # Fit function over the interval
#             coeffs, fit_stats           = fit_integral_2d( fit_function, fit_props, show_figs=show_figs )
#             ref_level_i.cheby_coeffs    = coeffs
#             ref_level_i.fit_stats       = fit_stats

#             # # Check fit tolerances and refine if any
#             # if not ref_level_i.check_tolerances( ):
#             #     fit_residual_2D_adaptive_interface( fit_function, ref_level_i, show_figs=show_figs )
            
#             # Add ith refinement level region to the children list
#             ref_levels_x.append( ref_level_i )

#         # Loop over refinements to create new refinement levels
#         ref_levels_y = [ ]
#         for j in range( y_new.shape[0] - 1 ):
#             # Apply new sub-region boundaries
#             fit_props       = copy.deepcopy( ref_level.fit_props )
#             fit_props.y_max = y_new[j+1]
#             fit_props.y_min = y_new[j]

#             # Create new refinement level
#             ref_level_i     = RefLevel( fit_props, parent=ref_level )

#             # Fit function over the interval
#             coeffs, fit_stats           = fit_integral_2d( fit_function, fit_props, show_figs=show_figs )
#             ref_level_i.cheby_coeffs    = coeffs
#             ref_level_i.fit_stats       = fit_stats

#             # # Check fit tolerances and refine if any
#             # if not ref_level_i.check_tolerances( ):
#             #     fit_residual_2D_adaptive_interface( fit_function, ref_level_i, show_figs=show_figs )
            
#             # Add ith refinement level region to the children list
#             ref_levels_y.append( ref_level_i )

#         # Check which refinement level is more expensive
#         num_coeffs_x = 0
#         for rl in ref_levels_x:
#             num_coeffs_x += np.array( rl.get_num_cheby_coeffs( ) ).sum( )
        
#         num_coeffs_y = 0
#         for rl in ref_levels_y:
#             num_coeffs_y += np.array( rl.get_num_cheby_coeffs( ) ).sum( )
        
#         if num_coeffs_x > num_coeffs_y:
#             for refli in ref_levels_y:
#                 if not refli.check_tolerances( ):
#                     fit_residual_2D_adaptive_interface( fit_function, refli, is_square_ref=is_square_ref, show_figs=show_figs )
            
#             ref_level.child.extend( ref_levels_y )


#         else:
#             for refli in ref_levels_x:
#                 if not refli.check_tolerances( ):
#                     fit_residual_2D_adaptive_interface( fit_function, refli, is_square_ref=is_square_ref, show_figs=show_figs )

#             ref_level.child.extend( ref_levels_x )


def fit_residual_region_11_dx(show_figs=False, show_summary_fig=False):
    # Define boundary values
    # x = array([1e-8, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1])
    # y = array([1e-8, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1])
    # x = 10**array( [ -5.0, 1.0, 2.0, 3.0, 4.0, 5.0 ] )
    # y = 10**array( [ -5.0, -1.0, log10( 1.0 ), log10( np.pi/2.0 - 1e-5 ) ] )
    x = 10**np.linspace( -5, 5, 60 )
    y = 10**np.linspace( -5, log10( np.pi/2.0 - 1e-5 ), 60 )
    results = []
    for i in range(x.shape[0]-1):
        for j in range(y.shape[0]-1):
            region_name = f"fit_residual_region_11_dx -> X: {x[i]:0.3E} - {x[i+1]:0.3E} | Y: {y[j]:0.3E} - {y[j+1]:0.3E}"

            # Set up domain extension and fit properties
            fit_props = FitProperties()
            fit_props.cheby_order_x = 5
            fit_props.cheby_order_y = 5
            fit_props.x_log_scale = True
            fit_props.x_max = x[i+1]
            fit_props.x_min = x[i]
            fit_props.y_log_scale = True
            fit_props.y_max = y[j+1]
            fit_props.y_min = y[j]
            fit_props.cheby_tol = FIT_TOL
            fit_props.x_map_fcn = lambda x: x
            fit_props.y_map_fcn = lambda y: y

            # Launch fit
            results.append(fit_integral_2d(residual_region_11_dx, fit_props, region_name, show_figs=show_figs, stop_to_show=False))

    # Plot results summay
    plot_results_summary(x, y, results, "region_11_dx", log_scale=True, show_figs=show_summary_fig)

    return x, y, results


def fit_residual_region_11A_dx(show_figs=False, show_summary_fig=False):
    # Define boundary values
    x = array([1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 2e-1, 3e-1, 5e-1])
    y = array([1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 2e-1, 3e-1, 5e-1])
    # x = array([1e-8, 1e-7])
    # y = array([1e-3, 3e-1])
    results = []
    for i in range(x.shape[0]-1):
        for j in range(y.shape[0]-1):
            region_name = f"fit_residual_region_11A_dx -> X: {x[i]:0.3E} - {x[i+1]:0.3E} | Y: {y[j]:0.3E} - {y[j+1]:0.3E}"

            # Set up domain extension and fit properties
            fit_props = FitProperties()
            fit_props.cheby_order_x = 20
            fit_props.cheby_order_y = 20
            fit_props.x_log_scale = True
            fit_props.x_max = x[i+1]
            fit_props.x_min = x[i]
            fit_props.y_log_scale = True
            fit_props.y_max = y[j+1]
            fit_props.y_min = y[j]
            fit_props.cheby_tol = FIT_TOL
            fit_props.x_map_fcn = lambda x: x
            fit_props.y_map_fcn = lambda y: y

            # Launch fit
            results.append(fit_integral_2d(residual_region_11A_dx, fit_props, region_name, show_figs=show_figs))

    # Plot results summay
    plot_results_summary(x, y, results, "region_11_dx", log_scale=True, show_figs=show_summary_fig)

    return x, y, results


def fit_residual_region_11B_dx(show_figs=False, show_summary_fig=False):
    # Define boundary values
    x = array([1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.3, 0.5, 1.0, 3.0])
    y = array([1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.3, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0])
    results = []
    for i in range(x.shape[0]-1):
        for j in range(y.shape[0]-1):
            region_name = f"fit_residual_region_11B_dx -> X: {x[i]:0.3E} - {x[i+1]:0.3E} | Y: {y[j]:0.3E} - {y[j+1]:0.3E}"

            # Set up domain extension and fit properties
            fit_props = FitProperties()
            fit_props.cheby_order_x = 20
            fit_props.cheby_order_y = 20
            fit_props.x_log_scale = False
            fit_props.x_max = x[i+1]
            fit_props.x_min = x[i]
            fit_props.y_log_scale = True
            fit_props.y_max = y[j+1]
            fit_props.y_min = y[j]
            fit_props.cheby_tol = FIT_TOL
            fit_props.x_map_fcn = lambda x: x
            fit_props.y_map_fcn = lambda y: y

            # Launch fit
            results_i = fit_integral_2d( residual_region_11B_dx, fit_props, region_name, show_figs=show_figs )

            # Overwrite coefficients data if any
            # if x[i] < 0.3 and y[j] < 0.3:
            #     results_i = [
            #                     array([0.0]),
            #                     array([0], dtype=int),
            #                     array([0], dtype=int),
            #                     results_i[3],
            #                     results_i[4]
            #                 ]
            
            results.append( results_i )

    # Plot results summay
    plot_results_summary(x, y, results, "region_11_dx", log_scale=True, show_figs=show_summary_fig)

    return x, y, results


def fit_residual_region_12(show_figs=False, show_summary_fig=False):
    # Define boundary values
    x = array([1e-12, 1e-3, 3.0])
    y = array([4.0, 12.0, 50.0, 200.0, 1000.0, 2000.0, 10000, 20000.0, 100000.0])
    results = []
    for i in range(x.shape[0]-1):
        for j in range(y.shape[0]-1):
            region_name = f"fit_residual_region_12 -> X: {x[i]:0.3E} - {x[i+1]:0.3E} | Y: {y[j]:0.3E} - {y[j+1]:0.3E}"

            # Set up domain extension and fit properties
            fit_props = FitProperties()
            fit_props.cheby_order_x = 20
            fit_props.cheby_order_y = 20
            fit_props.x_log_scale = False
            fit_props.x_max = x[i+1]
            fit_props.x_min = x[i]
            fit_props.y_log_scale = False
            fit_props.y_max = y[j+1]
            fit_props.y_min = y[j]
            fit_props.cheby_tol = FIT_TOL
            fit_props.x_map_fcn = lambda x: x
            fit_props.y_map_fcn = lambda y: y

            # Launch fit
            print("Working on: ", region_name)
            results.append(fit_integral_2d(residual_region_12, fit_props, region_name, show_figs=show_figs))

    # Plot results summay
    plot_results_summary(x, y, results, "region_12", log_scale=True, show_figs=show_summary_fig)

    return x, y, results


def fit_residual_region_12_dx(show_figs=False, show_summary_fig=False):
    # Define boundary values
    x = array([1e-8, 1e-7, 1e-6, 1e-7, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 3.0])
    y = array([4.0, 12.0, 37.0, 120.0, 200.0, 500.0, 1000.0])
    results = []
    for i in range(x.shape[0]-1):
        for j in range(y.shape[0]-1):
            region_name = f"fit_residual_region_12_dx -> X: {x[i]:0.3E} - {x[i+1]:0.3E} | Y: {y[j]:0.3E} - {y[j+1]:0.3E}"

            # Set up domain extension and fit properties
            fit_props = FitProperties()
            fit_props.cheby_order_x = 20
            fit_props.cheby_order_y = 20
            fit_props.x_log_scale = False
            fit_props.x_max = x[i+1]
            fit_props.x_min = x[i]
            fit_props.y_log_scale = False
            fit_props.y_max = y[j+1]
            fit_props.y_min = y[j]
            fit_props.cheby_tol = FIT_TOL
            fit_props.x_map_fcn = lambda x: x
            fit_props.y_map_fcn = lambda y: y

            # Launch fit
            print("Working on: ", region_name)
            results_i = fit_integral_2d(residual_region_12_dx, fit_props, region_name, show_figs=show_figs)

            # Check if the value of the function is above a prescribed limit or not
            if results_i[3].max_value < -8:
                print(" WARNING: \n - Maximum value is below the threshold to take into account the fit coeffs.")
                results_i = [
                            array([0.0]),
                            array([0], dtype=int),
                            array([0], dtype=int),
                            results_i[3],
                            results_i[4]
                            ]

            results.append(results_i)

    # Plot results summay
    plot_results_summary(x, y, results, "region_12_dx", log_scale=True, show_figs=show_summary_fig)

    return x, y, results


def fit_residual_region_21(show_figs=False, show_summary_fig=False):
    # Define boundary values
    x = array([3.0, 12.0, 50.0, 200.0, 1000.0, 2000.0, 10000.0, 20000.0, 100000.0])
    y = array([1e-12, 1e-3, 4.0])
    results = []
    for i in range(x.shape[0]-1):
        for j in range(y.shape[0]-1):
            region_name = f"fit_residual_region_21 -> X: {x[i]:0.3E} - {x[i+1]:0.3E} | Y: {y[j]:0.3E} - {y[j+1]:0.3E}"

            # Set up domain extension and fit properties
            fit_props = FitProperties()
            fit_props.cheby_order_x = 20
            fit_props.cheby_order_y = 20
            fit_props.x_log_scale = False
            fit_props.x_max = x[i+1]
            fit_props.x_min = x[i]
            fit_props.y_log_scale = False
            fit_props.y_max = y[j+1]
            fit_props.y_min = y[j]
            fit_props.cheby_tol = FIT_TOL
            fit_props.x_map_fcn = lambda x: x
            fit_props.y_map_fcn = lambda y: y

            # Launch fit
            print("Working on: ", region_name)
            results.append(fit_integral_2d(residual_region_21, fit_props, region_name, show_figs=show_figs))

    # Plot results summay
    plot_results_summary(x, y, results, "region_21", log_scale=True, show_figs=show_summary_fig)

    return x, y, results


def fit_residual_region_21_dx(show_figs=False, show_summary_fig=False):
    # Define boundary values
    x = array([3.0, 12.0, 50.0, 200.0, 500.0, 1000.0])
    y = array([1e-12, 1e-3, 4.0])
    results = []
    for i in range(x.shape[0]-1):
        for j in range(y.shape[0]-1):
            region_name = f"fit_residual_region_21_dx -> X: {x[i]:0.3E} - {x[i+1]:0.3E} | Y: {y[j]:0.3E} - {y[j+1]:0.3E}"

            # Set up domain extension and fit properties
            fit_props = FitProperties()
            fit_props.cheby_order_x = 20
            fit_props.cheby_order_y = 20
            fit_props.x_log_scale = False
            fit_props.x_max = x[i+1]
            fit_props.x_min = x[i]
            fit_props.y_log_scale = False
            fit_props.y_max = y[j+1]
            fit_props.y_min = y[j]
            fit_props.cheby_tol = FIT_TOL
            fit_props.x_map_fcn = lambda x: x
            fit_props.y_map_fcn = lambda y: y

            # Launch fit
            print("Working on: ", region_name)
            results.append(fit_integral_2d(residual_region_21_dx, fit_props, region_name, show_figs=show_figs))

    # Plot results summay
    plot_results_summary(x, y, results, "region_21_dx", log_scale=True, show_figs=show_summary_fig)

    return x, y, results


def fit_residual_region_22(show_figs=False, show_summary_fig=False):
    # Define boundary values
    x = array([3.0, 6.0, 12.0, 30.0, 50.0, 120.0, 200.0, 500.0, 1000.0, 2000.0, 10000.0, 20000.0, 100000.0])
    y = array([3.0, 6.0, 12.0, 30.0, 50.0, 120.0, 200.0, 500.0, 1000.0, 2000.0, 10000.0, 20000.0, 100000.0])
    results = []
    for i in range(x.shape[0]-1):
        for j in range(y.shape[0]-1):
            region_name = f"fit_residual_region_22 -> X: {x[i]:0.3E} - {x[i+1]:0.3E} | Y: {y[j]:0.3E} - {y[j+1]:0.3E}"

            # Set up domain extension and fit properties
            fit_props = FitProperties()
            fit_props.cheby_order_x = 20
            fit_props.cheby_order_y = 20
            fit_props.x_log_scale = False
            fit_props.x_max = x[i+1]
            fit_props.x_min = x[i]
            fit_props.y_log_scale = False
            fit_props.y_max = y[j+1]
            fit_props.y_min = y[j]
            fit_props.cheby_tol = FIT_TOL
            fit_props.x_map_fcn = lambda x: x
            fit_props.y_map_fcn = lambda y: y

            # Launch fit
            results.append(fit_integral_2d(residual_region_22, fit_props, region_name, show_figs=show_figs))

    # Plot results summay
    plot_results_summary(x, y, results, "region_22", log_scale=False, show_figs=show_summary_fig)

    return x, y, results


def fit_residual_region_22_dx(show_figs=False, show_summary_fig=False):
    # Define boundary values
    # x = array([3.0, 6.0, 12.0, 30.0, 50.0, 200.0, 500.0, 1000.0, 2000.0, 10000.0, 20000.0, 100000.0])
    # y = array([3.0, 6.0, 12.0, 30.0, 50.0, 200.0, 500.0, 1000.0, 2000.0, 10000.0, 20000.0, 100000.0])
    x = array([3.0, 12.0, 50.0, 200.0, 1000.0, 10000.0, 100000.0])
    y = array([3.0, 12.0, 50.0, 200.0, 1000.0, 10000.0, 100000.0])
    results = []
    for i in range(x.shape[0]-1):
        for j in range(y.shape[0]-1):
            region_name = f"fit_residual_region_22_dx -> X: {x[i]:0.3E} - {x[i+1]:0.3E} | Y: {y[j]:0.3E} - {y[j+1]:0.3E}"

            # Set up domain extension and fit properties
            fit_props = FitProperties()
            fit_props.cheby_order_x = 20
            fit_props.cheby_order_y = 20
            fit_props.x_log_scale = True
            fit_props.x_max = x[i+1]
            fit_props.x_min = x[i]
            fit_props.y_log_scale = True
            fit_props.y_max = y[j+1]
            fit_props.y_min = y[j]
            fit_props.cheby_tol = FIT_TOL
            fit_props.x_map_fcn = lambda x: x
            fit_props.y_map_fcn = lambda y: y

            # Launch fit
            print("Working on Region: ", region_name)
            results.append(fit_integral_2d(residual_region_22_dx, fit_props, region_name, show_figs=show_figs))

    # Plot results summay
    plot_results_summary(x, y, results, "region_22_dx", log_scale=True, show_figs=show_summary_fig)

    return x, y, results


def residual_region_0(X: float, Y: float)->float:
    return X**2.0+Y**2.0


def residual_region_11(X: float, Y: float)->float:
    R = sqrt(X**2.0+Y**2.0)
    c0 = (
        -2*jv(0, X)*log(R+Y)
        -(pi*yv(0, X)-2*jv(0, X)*log(X))
        )*exp(-Y)
    
    # return (fxy(X, Y)-c0)*exp(Y)/R
    return fxy(X, Y, only_int=True)
    # return np.log10( fxy(X, Y, only_int=True) )
    # return fxy_polar(X, Y, only_int=True)


def residual_region_11_dx(X: float, Y: float)->float:
    R = sqrt(X**2.0+Y**2.0)
    c0 = (
        +2*jv(1, X)*log(R+Y)
        -2*jv(0, X)*X/R/(R+Y)
        +pi*yv(1, X)
        -2*jv(1, X)*log(X)
        +2*jv(0, X)/X
        )*exp(-Y)
    
    # return (fxy_dx(X, Y)-c0)*exp(Y)/R
    return fxy_dx(X, Y, only_int=True)
    # return log10( fxy_dx_polar(X, Y, only_int=True) )
    # return log10(-fxy_dx(X, Y, only_int=False))
    # return fxy_dx(X, Y, only_int=False)


def residual_region_AA_dx(X: float, Y: float)->float:
    return X**2.0 + Y**2.0


def residual_region_11A_dx(X: float, Y: float)->float:
    return log10(fxy_dx(X, Y, only_int=True))


def residual_region_11B_dx(X: float, Y: float)->float:
    R = sqrt(X**2.0+Y**2.0)
    c0 = (
        +2*jv(1, X)*log(R+Y)
        -2*jv(0, X)*X/R/(R+Y)
        +pi*yv(1, X)
        -2*jv(1, X)*log(X)
        +2*jv(0, X)/X
        )*exp(-Y)
    
    return (fxy_dx(X, Y)-c0)*exp(Y)


def residual_region_12(X: float, Y: float)->float:
    return fxy(X, Y)


def residual_region_12_dx(X: float, Y: float)->float:
    return fxy_dx(X, Y, only_int=False)


def residual_region_21(X: float, Y: float)->float:
    return fxy(X, Y)+2*pi*exp(-Y)*yv(0, X)


def residual_region_21_dx(X: float, Y: float)->float:
    return fxy_dx(X, Y)-2*pi*exp(-Y)*yv(1, X)


def residual_region_22(X: float, Y: float)->float:
    R = sqrt(X**2.0+Y**2.0)
    costh = Y/R
    f = fxy(X, Y) + 2*pi*exp(-Y)*yv(0, X)
    for i in range(3):
        f -= factorial(i)*eval_legendre(i, costh)/R**(i+1)
    return f


def residual_region_22_dx(X: float, Y: float)->float:
    R = sqrt(X**2.0+Y**2.0)
    costh = Y/R
    costhd = -X*Y/R**3.0
    f = fxy_dx(X, Y) - 2*pi*exp(-Y)*yv(1, X)
    for i in range(3):
        poly = legendre(i)
        polyd = poly.deriv()
        f -= factorial(i)*(
                            polyval(polyd, costh)*costhd/R**(i+1)
                            -
                            eval_legendre(i, costh)*2*(i+1)*X/R**(i+2)
                            )
    return f


def plot_results_summary(x: ndarray, y: ndarray, results: list, region_name: str,
                        log_scale=False, show_figs=False)->None:
    def _plot_heat_map(axi, data, log_scale_i)->None:
        imis = axi.imshow(data, origin="lower", cmap="jet", extent=plot_extent)
        axi.set_xticks(x)
        axi.set_yticks(y)
        plt.colorbar(imis, ax=axi)

        if log_scale_i:
            axi.set_xscale("log")
            axi.set_yscale("log")

    # Arrange results in matrix order to plot as 
    # heat map
    num_x = x.shape[0]-1
    num_y = y.shape[0]-1
    err_max = zeros((num_x, num_y))
    err_mean = zeros((num_x, num_y))
    err_min = zeros((num_x, num_y))
    err_over_thr = zeros((num_x, num_y))
    num_coeffs = zeros((num_x, num_y))
    for i in range(num_x):
        for j in range(num_y):
            err_max[i, j] = results[i*num_y+j][3].max_err
            err_mean[i, j] = results[i*num_y+j][3].mean_err
            err_min[i, j] = results[i*num_y+j][3].min_err
            err_over_thr[i, j] = results[i*num_y+j][3].err_over_thr
            num_coeffs[i, j] = results[i*num_y+j][3].num_coeffs

    # Plot results
    plot_extent = [
                    x.min(),
                    x.max(),
                    y.min(),
                    y.max()
                    ]

    fig = plt.figure()
    fig.suptitle(region_name)
    ax0 = fig.add_subplot(231)
    ax0.set_title("Max. Error")
    ax1 = fig.add_subplot(234)
    ax1.set_title("Min. Error")
    ax2 = fig.add_subplot(232)
    ax2.set_title("Mean. Error")
    ax3 = fig.add_subplot(235)
    ax3.set_title("Values Over Threshold")
    ax4 = fig.add_subplot(233)
    ax4.set_title("Num. Coefficients")

    _plot_heat_map(ax0, err_max.T, log_scale)
    _plot_heat_map(ax1, err_min.T, log_scale)
    _plot_heat_map(ax2, err_mean.T, log_scale)
    _plot_heat_map(ax3, err_over_thr.T, log_scale)
    _plot_heat_map(ax4, num_coeffs.T, log_scale)

    if show_figs:
        plt.show()


if __name__ == "__main__":
    # fit_residual_region_test()

    # fit_residual_region_11(show_figs=True)
    # fit_residual_region_11_subregion_11(show_figs=True)
    # fit_residual_region_11_subregion_12(show_figs=True)
    # fit_residual_region_11_subregion_21(show_figs=True)
    # fit_residual_region_11_subregion_22(show_figs=True)

    # fit_residual_region_11_dx(show_figs=True)

    # fit_residual_region_11_subregion_11_dx(show_figs=True)
    # fit_residual_region_11_subregion_12_dx(show_figs=True)
    # fit_residual_region_11_subregion_21_dx(show_figs=True)
    # fit_residual_region_11_subregion_22_dx(show_figs=True)

    # fit_residual_region_12(show_figs=True)

    # fit_residual_region_12_subregion_11(show_figs=True)
    # fit_residual_region_12_subregion_21(show_figs=True)
    # fit_residual_region_12_subregion_12(show_figs=True)
    # fit_residual_region_12_subregion_22(show_figs=True)
    # fit_residual_region_12_subregion_13(show_figs=True)
    # fit_residual_region_12_subregion_23(show_figs=True)
    # fit_residual_region_12_subregion_14(show_figs=True)
    # fit_residual_region_12_subregion_15(show_figs=True)
    # fit_residual_region_12_subregion_24(show_figs=True)

    # fit_residual_region_12_dx(show_figs=True)

    # fit_residual_region_12_subregion_11_dx(show_figs=True)
    # fit_residual_region_12_subregion_12_dx(show_figs=True)
    # fit_residual_region_12_subregion_13_dx(show_figs=True)
    # fit_residual_region_12_subregion_14_dx(show_figs=True)
    # fit_residual_region_12_subregion_21_dx(show_figs=True)
    # fit_residual_region_12_subregion_22_dx(show_figs=True)
    # fit_residual_region_12_subregion_23_dx(show_figs=True)
    # fit_residual_region_12_subregion_24_dx(show_figs=True)

    # fit_residual_region_21(show_figs=True)

    # fit_residual_region_21_subregion_11(show_figs=True)
    # fit_residual_region_21_subregion_12(show_figs=True)
    # fit_residual_region_21_subregion_21(show_figs=True)
    # fit_residual_region_21_subregion_22(show_figs=True)
    # fit_residual_region_21_subregion_31(show_figs=True)
    # fit_residual_region_21_subregion_32(show_figs=True)
    # fit_residual_region_21_subregion_41(show_figs=True)
    # fit_residual_region_21_subregion_42(show_figs=True)

    # fit_residual_region_21_dx(show_figs=True)

    # fit_residual_region_21_subregion_11_dx(show_figs=True)
    # fit_residual_region_21_subregion_12_dx(show_figs=True)
    # fit_residual_region_21_subregion_21_dx(show_figs=True)
    # fit_residual_region_21_subregion_22_dx(show_figs=True)
    # fit_residual_region_21_subregion_31_dx(show_figs=True)
    # fit_residual_region_21_subregion_32_dx(show_figs=True)
    # fit_residual_region_21_subregion_41_dx(show_figs=True)
    # fit_residual_region_21_subregion_42_dx(show_figs=True)

    # fit_residual_region_11(show_figs=True)
    # fit_residual_region_11A_dx(show_figs=False)
    # fit_residual_region_11B_dx(show_figs=True)
    # fit_residual_region_22(show_figs=True)
    # fit_residual_region_12_dx(show_figs=True)
    # fit_residual_region_22_dx(show_figs=True)

    # fit_residual_region_22_subregion_11(show_figs=True)

    # fit_residual_region_12()
    # fit_residual_region_11_subregion_21()
    # fit_residual_region_11_subregion_22(show_figs=True)
    # fit_residual_region_21_subregion_11(show_figs=True)
    # fit_residual_region_21_subregion_12(show_figs=True)
    # fit_residual_region_21_subregion_21(show_figs=True)

    # X = 0.001
    # Y = 0.116384
    # ri = residual_region_11(X, Y)
    # ri = residual_region_11A_dx(X, Y)
    # ri = residual_region_22_dx(4.0, 5.0)
    # print("residual_region_11A_dx: ", ri, 10**ri)
    # fit_residual_region_11A_dx(show_figs=True, show_summary_fig=True)
    # residual_region_11B_dx( X, Y )
    
    # x       = array([1e-8, 2.0e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 3e-1])
    # y       = array([1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 3e-1])
    # y       = array([1e-7, 2e-7, 3e-7, 4e-7, 5e-7, 6e-7, 7e-7, 8e-7, 9e-7, 1e-6])
    # from pysr import PySRRegressor

    fit_residual_region_1d_test( )
    fit_residual_region_2d_test( )

    raise ValueError( "Stop by user" )
    


    # x       = array([1e-8, 1.25e-8, 1.5e-8, 1.75e-8, 2.0e-8])
    N       = 100
    x       = 10**linspace( -3, 5, N )
    y       = 10**linspace( -3, 5, N )
    # x       = 10**linspace( -5, 5, N )
    # y       = 10**linspace( -5, 5, N )
    data    = zeros( ( x.shape[0], y.shape[0] ) )

    X,Y     = np.meshgrid( x, y, indexing="ij" )

    fig     = plt.figure( )
    ax      = fig.add_subplot( 111 )
    for i in range( x.shape[0] ):
        for j in range( y.shape[0] ):
            # data[i, j] = fxy( x[i], y[j], only_int=True )
            data[i, j] = fxy_dx( x[i], y[j], only_int=True )
            # data[i, j] = residual_region_11_dx( x[i], y[j] )

        ax.plot( log10( y ), data[i, : ], label=f"X: {x[i]:0.2E}" )

    X_flat  = log10( X.ravel( ) )
    Y_flat  = log10( Y.ravel( ) )
    F_flat  = data.ravel( )
    XY      = np.column_stack( [ X_flat, Y_flat ] )
    # X = log10( np.column_stack( [ x, y ] ) )
    # ax.plot( log10( y ), -1/y**0.98, label="LOG" )
    # def inv(x):
    #     return 1.0 / x

    # model   = PySRRegressor(
    #                         niterations=2000,
    #                         binary_operators=["+", "-", "*", "/", "^"],
    #                         unary_operators=["exp", "log", "inv", "square", "cube", "sqrt", "cbrt"],
    #                         model_selection="best"
    #                         )

    # model.fit( log10( y.reshape(-1, 1) ), data[[0], :].T )
    # model.fit( XY, F_flat )

    # square = lambda x: x*x
    # cube = lambda x: x*x*x

    print( "Best Model:" )
    # print( model.get_best( ) )
    x0 = X_flat
    x1 = Y_flat
    # f = (x0 * (exp((sin(x0) * (x0 + 2.2302966)) + 7.5992866) + ((x0 * 13.14788) + x0))) * 6.1577554
    # # f = ((exp((x0 * -4.203177) + -17.442085) + 82.86912) + x0) * ((x0 + x0) + 1.552544)
    # f = x0 * (exp((x0 * -4.4348288) + -18.490313) + -41.41545)
    # f = x0 * ((inv(x1 + x0) * -3.4373405) * inv(x1 + x0))
    # f = exp(x1 * -1.8649621) * x0
    # f = (inv(exp(x0 + x0) + exp(((x1 + 3.9225373) * 0.87918144) * (x0 + 12.257254))) * x0) * 2.126287
    # f = (inv(exp(exp(x0) + (x1 * 5.726273))) * -450.02847) * (exp((((x0 + -3.1291173) * exp(((x1 * -0.83531713) + x0) + -1.1155022)) + x0) * 3.7475014) + 4.055301e-16)
    # f = (cube(x0) + 57.676624) * exp((sqrt(erfc(x1 + ((x0 * 0.108337566) - x0))) * x0) * -1.1298819)
    # f = (x0 + 1.2383021) * ((0.744047 - square(erf(0.03381948 - erf(x0 - x1)))) / (cube(erf(exp(x0))) + (125.58079 ** (x1 / 2.3112185))))
    # f = f.reshape( N, N )
    # model_id = model.get_best( )[ "complexity" ]
    # f = model.predict( XY ).reshape( N, N )
    # for i in range( x.shape[0] ):
    #     for j in range( y.shape[0] ):
    #         ax.plot( log10( y ), f[i, : ], "--", color="grey" )

    # fig = plt.figure( )
    # ax  = fig.add_subplot( 111 )
    
    # for i in range( x.shape[0] ):
    #     for j in range( y.shape[0] ):
    #         ax.plot( log10( y ), f[i, : ]-data[i, :], "--", color="grey" )
    # ax.plot( x0, f, label="PySR" )

    # X = 2 * np.random.randn(100, 5)
    # y = 1 / X[:, [0, 1, 2]]
    # model = PySRRegressor(
    #     binary_operators=["+", "*"],
    #     unary_operators=["inv(x) = 1/x"],
    #     extra_sympy_mappings={"inv": lambda x: 1/x},
    # )
    # model.fit(X, y)

    # print( model )
    # print( "Best Model:" )
    # print( model.get_best( )[ "complexity" ] )
    
    ax.legend( )

    fig = plt.figure( )
    ax  = fig.add_subplot( 111 )
    cnf = ax.contourf( log10( X ), log10( Y ), data )
    plt.colorbar( cnf, ax=ax ) 

    fig = plt.figure( )
    fig.suptitle( "LOG Scale" )
    ax  = fig.add_subplot( 111 )
    cnf = ax.contourf( log10( X ), log10( Y ), log10( data ) )
    plt.colorbar( cnf, ax=ax ) 

    fig = plt.figure( )
    ax  = fig.add_subplot( 111, projection="3d" )
    
    psf = ax.plot_surface( log10( X ), log10( Y ), data, cmap="jet")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    plt.colorbar(psf, ax=ax)
    
    plt.show( )
    