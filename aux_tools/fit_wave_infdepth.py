
#
# Copyright (c) 2025 Sergio Fernández Ruano / IHCantabria
#
# This file is part of SeaMotions Software.
#
# SeaMotions is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SeaMotions is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#

# Import general usage libraries
import copy
from typing import Callable

# Import general usage scientific libraries
import matplotlib.pyplot as plt
import numpy as np

# Import local modules
from base_integrals import fxy, fxy_dx
from fit_cheby import FitProperties, fit_integral_2d, RefLevel, fit_residual_2D_adaptive_interface


def fit_residual_region_00_dx_adaptive( folder_path: str, show_figs=False, is_square_ref=True, show_summary_fig=False ) -> None:
    # Define boundary values
    fit_props               = FitProperties( )
    fit_props.dims          = 2
    fit_props.region_name   = "region_00_dx"
    fit_props.cheby_order_x = 10
    fit_props.cheby_order_y = 10
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = True
    fit_props.x_max         = 1.0
    fit_props.x_min         = -5.0
    fit_props.y_log_scale   = True
    fit_props.y_max         = 1.0
    fit_props.y_min         = -5.0
    fit_props.cheby_abs_tol = 5E-7
    fit_props.cheby_rel_tol = 1E-14
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.max_ref_level = 10

    fit_props.num_x         = 10
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = 10
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1

    fit_props.generate_fitting_matrix( )

    fit_function            = residual_region_00_dx

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


def fit_residual_region_11_adaptive( folder_path: str, show_figs=False, is_square_ref=True, show_summary_fig=False ) -> None:
    # Define boundary values
    fit_props               = FitProperties( )
    fit_props.dims          = 2
    fit_props.region_name   = "region_11"
    fit_props.cheby_order_x = 10
    fit_props.cheby_order_y = 10
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = True
    fit_props.x_max         = 5.0
    fit_props.x_min         = -8.0
    fit_props.y_log_scale   = True
    fit_props.y_max         = 5.0
    fit_props.y_min         = -8.0
    fit_props.cheby_abs_tol = 5E-7
    fit_props.cheby_rel_tol = 1E-14
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.max_ref_level = 10

    fit_props.num_x         = 10
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = 10
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1
    
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
    fit_props.cheby_order_x = 10
    fit_props.cheby_order_y = 10
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = True
    fit_props.x_max         = 5.0
    fit_props.x_min         = -5.0
    fit_props.y_log_scale   = True
    fit_props.y_max         = 5.0
    fit_props.y_min         = -5.0
    fit_props.cheby_abs_tol = 5E-7
    fit_props.cheby_rel_tol = 1E-14
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.max_ref_level = 10

    fit_props.num_x         = 10
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = 10
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1

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


def residual_region_00_dx( X: float, Y: float ) -> float:
    f = fxy_dx( X, Y )
    f = 1e-7*np.sign( f ) if np.abs( f ) < 1e-7 else f
    return f


def residual_region_11( X: float, Y: float ) -> float:
    return fxy( X, Y, only_int=True )


def residual_region_11_dx( X: float, Y: float ) -> float:
    return fxy_dx( X, Y, only_int=True )
