
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
import h5py
import os

# Import general usage scientific libraries
from numpy import argsort, array, concatenate, log, ndarray
from numpy import abs as np_abs
import numpy as np

# Import local modules
from base_integrals import L0, L0_dA, L0_dB, L1, L1_dA, L1_dB, L2, L3, L3_dA, L3_dB, Linf, Linf_dA, Linf_dB, M1, M1_dA, M1_dB, M2, M3, M3_dA, M3_dB
from fit_cheby import fit_integral_1d, fit_integral_2d, fit_integral_3d, FitProperties, FitStats, fit_residual_1D_adaptive_interface, fit_residual_2D_adaptive_interface, fit_residual_3D_adaptive_interface
from fit_tools import RefLevel, write_coeffs_module_adaptive_1d_only_header, write_coeffs_module_adaptive_2d_only_header, write_coeffs_module_adaptive_3d_only_header


FIT_TOL = 1E-5


def data_to_dict(
                fit_props: FitProperties,
                cheby_coeffs: ndarray,
                ncx: ndarray,
                ncy: ndarray,
                ncz: ndarray
                ):
    # Generate dict with the results data
    data = {
            "fit_props": fit_props,
            "cheby_coeffs": cheby_coeffs,
            "ncx": ncx,
            "ncy": ncy,
            "ncz": ncz
            }
    
    return data


def fit_L0( fopath: str ) -> None:
    ref_level = fit_residual_region_L0_adaptive( )

    write_coeffs_module_adaptive_2d_only_header( ref_level, fopath, "L0" )


def fit_L0_dA( fopath: str ) -> None:
    ref_level = fit_residual_region_L0_dA_adaptive( )

    write_coeffs_module_adaptive_2d_only_header( ref_level, fopath, "L0_dA" )


def fit_L0_dB( fopath: str ) -> None:
    ref_level = fit_residual_region_L0_dB_adaptive( )

    write_coeffs_module_adaptive_2d_only_header( ref_level, fopath, "L0_dB" )


def fit_L1( fopath: str ) -> None:
    ref_level = fit_residual_region_L1_adaptive( )

    write_coeffs_module_adaptive_3d_only_header( ref_level, fopath, "L1" )


def fit_L1_dA( fopath: str ) -> None:
    ref_level = fit_residual_region_L1_dA_adaptive( )

    write_coeffs_module_adaptive_3d_only_header( ref_level, fopath, "L1_dA" )


def fit_L1_dB( fopath: str ) -> None:
    ref_level = fit_residual_region_L1_dB_adaptive( )

    write_coeffs_module_adaptive_3d_only_header( ref_level, fopath, "L1_dB" )


def fit_L2( fopath: str ) -> None:
    ref_level = fit_residual_region_L2_adaptive( )

    write_coeffs_module_adaptive_1d_only_header( ref_level, fopath, "L2" )


def fit_L3( fopath: str ) -> None:
    ref_level = fit_residual_region_L3_adaptive( )

    write_coeffs_module_adaptive_3d_only_header( ref_level, fopath, "L3" )


def fit_L3_dA( fopath: str ) -> None:
    ref_level = fit_residual_region_L3_dA_adaptive( )

    write_coeffs_module_adaptive_3d_only_header( ref_level, fopath, "L3_dA" )


def fit_L3_dB( fopath: str ) -> None:
    ref_level = fit_residual_region_L3_dB_adaptive( )

    write_coeffs_module_adaptive_3d_only_header( ref_level, fopath, "L3_dB" )


def fit_Linf( fopath: str ) -> None:
    ref_level = fit_residual_region_Linf_adaptive( )

    write_coeffs_module_adaptive_2d_only_header( ref_level, fopath, "Linf" )


def fit_Linf_dA( fopath: str ) -> None:
    ref_level = fit_residual_region_Linf_dA_adaptive( )
    write_coeffs_module_adaptive_2d_only_header( ref_level, fopath, "Linf_dA" )


def fit_Linf_dB( fopath: str ) -> None:
    ref_level = fit_residual_region_Linf_dB_adaptive( )
    write_coeffs_module_adaptive_2d_only_header( ref_level, fopath, "Linf_dB" )


def fit_M1( fopath: str ) -> None:
    ref_level = fit_residual_region_M1_adaptive( )

    write_coeffs_module_adaptive_3d_only_header( ref_level, fopath, "M1" )


def fit_M1_dA( fopath: str ) -> None:
    ref_level = fit_residual_region_M1_dA_adaptive( )

    write_coeffs_module_adaptive_3d_only_header( ref_level, fopath, "M1_dA" )


def fit_M1_dB( fopath: str ) -> None:
    ref_level = fit_residual_region_M1_dB_adaptive( )

    write_coeffs_module_adaptive_3d_only_header( ref_level, fopath, "M1_dB" )


def fit_M2( fopath: str ) -> None:
    ref_level = fit_residual_region_M2_adaptive( )

    write_coeffs_module_adaptive_1d_only_header( ref_level, fopath, "M2" )


def fit_M3( fopath: str ) -> None:
    ref_level = fit_residual_region_M3_adaptive( )

    write_coeffs_module_adaptive_3d_only_header( ref_level, fopath, "M3" )


def fit_M3_dA( fopath: str ) -> None:
    ref_level = fit_residual_region_M3_dA_adaptive( )

    write_coeffs_module_adaptive_3d_only_header( ref_level, fopath, "M3_dA" )


def fit_M3_dB( fopath: str ) -> None:
    ref_level = fit_residual_region_M3_dB_adaptive( )

    write_coeffs_module_adaptive_3d_only_header( ref_level, fopath, "M3_dB" )


def fit_residual_region_L0_adaptive( ) -> None:
    fit_props               = FitProperties( )
    fit_props.dims          = 2
    fit_props.region_name   = "L0"
    fit_props.cheby_order_x = 10
    fit_props.cheby_order_y = 10
    fit_props.cheby_order_z = 50
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = False
    fit_props.x_max         = 1.0
    fit_props.x_min         = 0.0
    fit_props.y_log_scale   = False
    fit_props.y_max         = 1.0
    fit_props.y_min         = 0.0
    fit_props.z_log_scale   = True
    fit_props.z_max         = 0.0
    fit_props.z_min         = -5.0
    fit_props.cheby_abs_tol = 5E-7
    fit_props.cheby_rel_tol = 1E-14
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.z_map_fcn     = lambda z: z
    fit_props.max_ref_level = 10

    fit_props.num_x         = fit_props.cheby_order_x + 1
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = fit_props.cheby_order_y + 1
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1
    fit_props.num_z         = fit_props.cheby_order_z + 1
    fit_props.num_z_fit     = fit_props.cheby_order_z + 1

    fit_props.generate_fitting_matrix( )

    fit_function            = lambda x, y: L0(x, y)[0].real

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_2d( fit_function, fit_props ) )

    # Start adaptive refinement loop
    fit_residual_2D_adaptive_interface( fit_function, ref_level )

    # Set starting position
    ref_level.set_start_index( 0 )
    ref_level.set_start_index_folded( 0 )


    return ref_level


def fit_residual_region_L0_dA_adaptive( ) -> None:
    fit_props               = FitProperties( )
    fit_props.dims          = 2
    fit_props.region_name   = "L0_dA"
    fit_props.cheby_order_x = 10
    fit_props.cheby_order_y = 10
    fit_props.cheby_order_z = 50
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = False
    fit_props.x_max         = 1.0
    fit_props.x_min         = 0.0
    fit_props.y_log_scale   = False
    fit_props.y_max         = 1.0
    fit_props.y_min         = 0.0
    fit_props.z_log_scale   = True
    fit_props.z_max         = 0.0
    fit_props.z_min         = -5.0
    fit_props.cheby_abs_tol = 5E-7
    fit_props.cheby_rel_tol = 1E-14
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.z_map_fcn     = lambda z: z
    fit_props.max_ref_level = 10

    fit_props.num_x         = fit_props.cheby_order_x + 1
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = fit_props.cheby_order_y + 1
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1
    fit_props.num_z         = fit_props.cheby_order_z + 1
    fit_props.num_z_fit     = fit_props.cheby_order_z + 1

    fit_props.generate_fitting_matrix( )

    fit_function            = lambda x, y: L0_dA(x, y)[0].real

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_2d( fit_function, fit_props ) )

    # Start adaptive refinement loop
    fit_residual_2D_adaptive_interface( fit_function, ref_level )

    # Set starting position
    ref_level.set_start_index( 0 )
    ref_level.set_start_index_folded( 0 )


    return ref_level


def fit_residual_region_L0_dB_adaptive( ) -> None:
    fit_props               = FitProperties( )
    fit_props.dims          = 2
    fit_props.region_name   = "L0_dB"
    fit_props.cheby_order_x = 10
    fit_props.cheby_order_y = 10
    fit_props.cheby_order_z = 50
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = False
    fit_props.x_max         = 1.0
    fit_props.x_min         = 0.0
    fit_props.y_log_scale   = False
    fit_props.y_max         = 1.0
    fit_props.y_min         = 0.0
    fit_props.z_log_scale   = True
    fit_props.z_max         = 0.0
    fit_props.z_min         = -5.0
    fit_props.cheby_abs_tol = 5E-7
    fit_props.cheby_rel_tol = 1E-14
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.z_map_fcn     = lambda z: z
    fit_props.max_ref_level = 10

    fit_props.num_x         = fit_props.cheby_order_x + 1
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = fit_props.cheby_order_y + 1
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1
    fit_props.num_z         = fit_props.cheby_order_z + 1
    fit_props.num_z_fit     = fit_props.cheby_order_z + 1

    fit_props.generate_fitting_matrix( )

    fit_function            = lambda x, y: L0_dB(x, y)[0].real

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_2d( fit_function, fit_props ) )

    # Start adaptive refinement loop
    fit_residual_2D_adaptive_interface( fit_function, ref_level )

    # Set starting position
    ref_level.set_start_index( 0 )
    ref_level.set_start_index_folded( 0 )


    return ref_level


def fit_residual_region_L1_adaptive( ) -> None:
    fit_props               = FitProperties( )
    fit_props.dims          = 3
    fit_props.region_name   = "L1"
    fit_props.cheby_order_x = 10
    fit_props.cheby_order_y = 10
    fit_props.cheby_order_z = 50
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = False
    fit_props.x_max         = 1.0
    fit_props.x_min         = 0.0
    fit_props.y_log_scale   = False
    fit_props.y_max         = 1.0
    fit_props.y_min         = 0.0
    fit_props.z_log_scale   = True
    fit_props.z_max         = 0.0
    fit_props.z_min         = -5.0
    fit_props.cheby_abs_tol = 5E-7
    fit_props.cheby_rel_tol = 1E-14
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.z_map_fcn     = lambda z: z
    fit_props.max_ref_level = 10

    fit_props.num_x         = fit_props.cheby_order_x + 1
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = fit_props.cheby_order_y + 1
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1
    fit_props.num_z         = fit_props.cheby_order_z + 1
    fit_props.num_z_fit     = fit_props.cheby_order_z + 1

    fit_props.generate_fitting_matrix( )

    fit_function            = lambda x, y, z: L1(x, y, z)[0].real

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_3d( fit_function, fit_props ) )

    # Start adaptive refinement loop
    fit_residual_3D_adaptive_interface( fit_function, ref_level )

    # Set starting position
    ref_level.set_start_index( 0 )
    ref_level.set_start_index_folded( 0 )


    return ref_level


def fit_residual_region_LXX_adaptive( ) -> None:
    fit_props               = FitProperties( )
    fit_props.dims          = 3
    fit_props.region_name   = "LXX"
    fit_props.cheby_order_x = 10
    fit_props.cheby_order_y = 10
    fit_props.cheby_order_z = 50
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = False
    fit_props.x_max         = 3.906E-01
    fit_props.x_min         = 3.867E-01
    fit_props.y_log_scale   = False
    fit_props.y_max         = 8.477E-01
    fit_props.y_min         = 8.438E-01
    fit_props.z_log_scale   = True
    fit_props.z_max         = np.log10( 4.068E+00 )
    fit_props.z_min         = np.log10( 3.854E+00 )
    fit_props.cheby_abs_tol = 5E-7
    fit_props.cheby_rel_tol = 1E-14
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.z_map_fcn     = lambda z: z
    fit_props.max_ref_level = 10

    fit_props.num_x         = fit_props.cheby_order_x + 1
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = fit_props.cheby_order_y + 1
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1
    fit_props.num_z         = fit_props.cheby_order_z + 1
    fit_props.num_z_fit     = fit_props.cheby_order_z + 1

    fit_props.generate_fitting_matrix( )

    fit_function            = lambda x, y, z: L1(x, y, z)[0].real

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_3d( fit_function, fit_props ) )

    # Start adaptive refinement loop
    fit_residual_3D_adaptive_interface( fit_function, ref_level )

    # Set starting position
    ref_level.set_start_index( 0 )
    ref_level.set_start_index_folded( 0 )


    return ref_level


def fit_residual_region_L1_dA_adaptive( ) -> None:
    fit_props               = FitProperties( )
    fit_props.dims          = 3
    fit_props.region_name   = "L1_dA"
    fit_props.cheby_order_x = 10
    fit_props.cheby_order_y = 10
    fit_props.cheby_order_z = 50
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = False
    fit_props.x_max         = 1.0
    fit_props.x_min         = 0.0
    fit_props.y_log_scale   = False
    fit_props.y_max         = 1.0
    fit_props.y_min         = 0.0
    fit_props.z_log_scale   = True
    fit_props.z_max         = 0.0
    fit_props.z_min         = -5.0
    fit_props.cheby_abs_tol = 5E-7
    fit_props.cheby_rel_tol = 1E-14
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.z_map_fcn     = lambda z: z
    fit_props.max_ref_level = 10

    fit_props.num_x         = fit_props.cheby_order_x + 1
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = fit_props.cheby_order_y + 1
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1
    fit_props.num_z         = fit_props.cheby_order_z + 1
    fit_props.num_z_fit     = fit_props.cheby_order_z + 1

    fit_props.generate_fitting_matrix( )

    fit_function            = lambda x, y, z: L1_dA(x, y, z)[0].real

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_3d( fit_function, fit_props ) )

    # Start adaptive refinement loop
    fit_residual_3D_adaptive_interface( fit_function, ref_level )

    # Set starting position
    ref_level.set_start_index( 0 )
    ref_level.set_start_index_folded( 0 )


    return ref_level


def fit_residual_region_L1_dB_adaptive( ) -> None:
    fit_props               = FitProperties( )
    fit_props.dims          = 3
    fit_props.region_name   = "L1_dB"
    fit_props.cheby_order_x = 10
    fit_props.cheby_order_y = 10
    fit_props.cheby_order_z = 50
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = False
    fit_props.x_max         = 1.0
    fit_props.x_min         = 0.0
    fit_props.y_log_scale   = False
    fit_props.y_max         = 1.0
    fit_props.y_min         = 0.0
    fit_props.z_log_scale   = True
    fit_props.z_max         = 0.0
    fit_props.z_min         = -5.0
    fit_props.cheby_abs_tol = 5E-7
    fit_props.cheby_rel_tol = 1E-14
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.z_map_fcn     = lambda z: z
    fit_props.max_ref_level = 10

    fit_props.num_x         = fit_props.cheby_order_x + 1
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = fit_props.cheby_order_y + 1
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1
    fit_props.num_z         = fit_props.cheby_order_z + 1
    fit_props.num_z_fit     = fit_props.cheby_order_z + 1

    fit_props.generate_fitting_matrix( )

    fit_function            = lambda x, y, z: L1_dB(x, y, z)[0].real

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_3d( fit_function, fit_props ) )

    # Start adaptive refinement loop
    fit_residual_3D_adaptive_interface( fit_function, ref_level )

    # Set starting position
    ref_level.set_start_index( 0 )
    ref_level.set_start_index_folded( 0 )


    return ref_level


def fit_residual_region_L2_adaptive( ) -> None:
    fit_props               = FitProperties( )
    fit_props.dims          = 1
    fit_props.region_name   = "L2"
    fit_props.cheby_order_x = 20
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = True
    fit_props.x_max         = 0.0
    fit_props.x_min         = -5.0
    fit_props.cheby_abs_tol = 5E-7
    fit_props.cheby_rel_tol = 1E-14
    fit_props.x_map_fcn     = lambda x: x
    fit_props.max_ref_level = 10

    fit_props.num_x         = fit_props.cheby_order_x + 1
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1

    fit_props.generate_fitting_matrix( )

    fit_function            = lambda x: L2( x )[0].real

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_1d( fit_function, fit_props ) )

    # Start adaptive refinement loop
    fit_residual_1D_adaptive_interface( fit_function, ref_level )

    # Set starting position
    ref_level.set_start_index( 0 )


    return ref_level


def fit_residual_region_L3_adaptive( ) -> None:
    fit_props               = FitProperties( )
    fit_props.dims          = 3
    fit_props.region_name   = "L3"
    fit_props.cheby_order_x = 10
    fit_props.cheby_order_y = 10
    fit_props.cheby_order_z = 50
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = False
    fit_props.x_max         = 1.0
    fit_props.x_min         = 0.0
    fit_props.y_log_scale   = False
    fit_props.y_max         = 1.0
    fit_props.y_min         = 0.0
    fit_props.z_log_scale   = True
    fit_props.z_max         = 3.0
    fit_props.z_min         = 0.0
    fit_props.cheby_abs_tol = 5E-7
    fit_props.cheby_rel_tol = 1E-14
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.z_map_fcn     = lambda z: z
    fit_props.max_ref_level = 10

    fit_props.num_x         = fit_props.cheby_order_x + 1
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = fit_props.cheby_order_y + 1
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1
    fit_props.num_z         = fit_props.cheby_order_z + 1
    fit_props.num_z_fit     = fit_props.cheby_order_z + 1

    fit_props.generate_fitting_matrix( )

    fit_function            = lambda x, y, z: L3(x, y, z)[0].real

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_3d( fit_function, fit_props ) )

    # Start adaptive refinement loop
    fit_residual_3D_adaptive_interface( fit_function, ref_level )

    # Set starting position
    ref_level.set_start_index( 0 )
    ref_level.set_start_index_folded( 0 )


    return ref_level


def fit_residual_region_L3_dA_adaptive( ) -> None:
    fit_props               = FitProperties( )
    fit_props.dims          = 3
    fit_props.region_name   = "L3_dA"
    fit_props.cheby_order_x = 10
    fit_props.cheby_order_y = 10
    fit_props.cheby_order_z = 50
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = False
    fit_props.x_max         = 1.0
    fit_props.x_min         = 0.0
    fit_props.y_log_scale   = False
    fit_props.y_max         = 1.0
    fit_props.y_min         = 0.0
    fit_props.z_log_scale   = True
    fit_props.z_max         = 3.0
    fit_props.z_min         = 0.0
    fit_props.cheby_abs_tol = 5E-6
    fit_props.cheby_rel_tol = 1E-14
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.z_map_fcn     = lambda z: z
    fit_props.max_ref_level = 10

    fit_props.num_x         = fit_props.cheby_order_x + 1
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = fit_props.cheby_order_y + 1
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1
    fit_props.num_z         = fit_props.cheby_order_z + 1
    fit_props.num_z_fit     = fit_props.cheby_order_z + 1

    fit_props.generate_fitting_matrix( )

    fit_function            = lambda x, y, z: L3_dA(x, y, z)[0].real

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_3d( fit_function, fit_props ) )

    # Start adaptive refinement loop
    fit_residual_3D_adaptive_interface( fit_function, ref_level )

    # Set starting position
    ref_level.set_start_index( 0 )
    ref_level.set_start_index_folded( 0 )


    return ref_level


def fit_residual_region_L3_dB_adaptive( ) -> None:
    fit_props               = FitProperties( )
    fit_props.dims          = 3
    fit_props.region_name   = "L3_dB"
    fit_props.cheby_order_x = 10
    fit_props.cheby_order_y = 10
    fit_props.cheby_order_z = 50
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = False
    fit_props.x_max         = 1.0
    fit_props.x_min         = 0.0
    fit_props.y_log_scale   = False
    fit_props.y_max         = 1.0
    fit_props.y_min         = 0.0
    fit_props.z_log_scale   = True
    fit_props.z_max         = 3.0
    fit_props.z_min         = 0.0
    fit_props.cheby_abs_tol = 5E-6
    fit_props.cheby_rel_tol = 1E-14
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.z_map_fcn     = lambda z: z
    fit_props.max_ref_level = 10

    fit_props.num_x         = fit_props.cheby_order_x + 1
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = fit_props.cheby_order_y + 1
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1
    fit_props.num_z         = fit_props.cheby_order_z + 1
    fit_props.num_z_fit     = fit_props.cheby_order_z + 1

    fit_props.generate_fitting_matrix( )

    fit_function            = lambda x, y, z: L3_dB(x, y, z)[0].real

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_3d( fit_function, fit_props ) )

    # Start adaptive refinement loop
    fit_residual_3D_adaptive_interface( fit_function, ref_level )

    # Set starting position
    ref_level.set_start_index( 0 )
    ref_level.set_start_index_folded( 0 )


    return ref_level


def fit_residual_region_Linf_adaptive( ) -> None:
    fit_props               = FitProperties( )
    fit_props.dims          = 2
    fit_props.region_name   = "Linf"
    fit_props.cheby_order_x = 10
    fit_props.cheby_order_y = 10
    fit_props.cheby_order_z = 50
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = False
    fit_props.x_max         = 1.0
    fit_props.x_min         = 0.0
    fit_props.y_log_scale   = False
    fit_props.y_max         = 1.0
    fit_props.y_min         = 0.0
    fit_props.z_log_scale   = True
    fit_props.z_max         = 0.0
    fit_props.z_min         = -5.0
    fit_props.cheby_abs_tol = 5E-7
    fit_props.cheby_rel_tol = 1E-14
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.z_map_fcn     = lambda z: z
    fit_props.max_ref_level = 10

    fit_props.num_x         = fit_props.cheby_order_x + 1
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = fit_props.cheby_order_y + 1
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1
    fit_props.num_z         = fit_props.cheby_order_z + 1
    fit_props.num_z_fit     = fit_props.cheby_order_z + 1

    fit_props.generate_fitting_matrix( )

    fit_function            = lambda x, y: Linf(x, y)[0].real

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_2d( fit_function, fit_props ) )

    # Start adaptive refinement loop
    fit_residual_2D_adaptive_interface( fit_function, ref_level )

    # Set starting position
    ref_level.set_start_index( 0 )
    ref_level.set_start_index_folded( 0 )


    return ref_level


def fit_residual_region_Linf_dA_adaptive( ) -> None:
    fit_props               = FitProperties( )
    fit_props.dims          = 2
    fit_props.region_name   = "Linf_dA"
    fit_props.cheby_order_x = 10
    fit_props.cheby_order_y = 10
    fit_props.cheby_order_z = 50
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = False
    fit_props.x_max         = 1.0
    fit_props.x_min         = 0.0
    fit_props.y_log_scale   = False
    fit_props.y_max         = 1.0
    fit_props.y_min         = 0.0
    fit_props.z_log_scale   = True
    fit_props.z_max         = 0.0
    fit_props.z_min         = -5.0
    fit_props.cheby_abs_tol = 5E-7
    fit_props.cheby_rel_tol = 1E-14
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.z_map_fcn     = lambda z: z
    fit_props.max_ref_level = 10

    fit_props.num_x         = fit_props.cheby_order_x + 1
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = fit_props.cheby_order_y + 1
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1
    fit_props.num_z         = fit_props.cheby_order_z + 1
    fit_props.num_z_fit     = fit_props.cheby_order_z + 1

    fit_props.generate_fitting_matrix( )

    fit_function            = lambda x, y: Linf_dA(x, y)[0].real

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_2d( fit_function, fit_props ) )

    # Start adaptive refinement loop
    fit_residual_2D_adaptive_interface( fit_function, ref_level )

    # Set starting position
    ref_level.set_start_index( 0 )
    ref_level.set_start_index_folded( 0 )


    return ref_level


def fit_residual_region_Linf_dB_adaptive( ) -> None:
    fit_props               = FitProperties( )
    fit_props.dims          = 2
    fit_props.region_name   = "Linf_dB"
    fit_props.cheby_order_x = 10
    fit_props.cheby_order_y = 10
    fit_props.cheby_order_z = 50
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = False
    fit_props.x_max         = 1.0
    fit_props.x_min         = 0.0
    fit_props.y_log_scale   = False
    fit_props.y_max         = 1.0
    fit_props.y_min         = 0.0
    fit_props.z_log_scale   = True
    fit_props.z_max         = 0.0
    fit_props.z_min         = -5.0
    fit_props.cheby_abs_tol = 5E-7
    fit_props.cheby_rel_tol = 1E-14
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.z_map_fcn     = lambda z: z
    fit_props.max_ref_level = 10

    fit_props.num_x         = fit_props.cheby_order_x + 1
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = fit_props.cheby_order_y + 1
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1
    fit_props.num_z         = fit_props.cheby_order_z + 1
    fit_props.num_z_fit     = fit_props.cheby_order_z + 1

    fit_props.generate_fitting_matrix( )

    fit_function            = lambda x, y: Linf_dB(x, y)[0].real

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_2d( fit_function, fit_props ) )

    # Start adaptive refinement loop
    fit_residual_2D_adaptive_interface( fit_function, ref_level )

    # Set starting position
    ref_level.set_start_index( 0 )
    ref_level.set_start_index_folded( 0 )


    return ref_level


def fit_residual_region_M1_adaptive( ) -> None:
    fit_props               = FitProperties( )
    fit_props.dims          = 3
    fit_props.region_name   = "M1"
    fit_props.cheby_order_x = 10
    fit_props.cheby_order_y = 10
    fit_props.cheby_order_z = 50
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = False
    fit_props.x_max         = 1.0
    fit_props.x_min         = 0.0
    fit_props.y_log_scale   = False
    fit_props.y_max         = 2.0
    fit_props.y_min         = 1.0
    fit_props.z_log_scale   = True
    fit_props.z_max         = 0.0
    fit_props.z_min         = -5.0
    fit_props.cheby_abs_tol = 5E-7
    fit_props.cheby_rel_tol = 1E-14
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.z_map_fcn     = lambda z: z
    fit_props.max_ref_level = 10

    fit_props.num_x         = fit_props.cheby_order_x + 1
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = fit_props.cheby_order_y + 1
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1
    fit_props.num_z         = fit_props.cheby_order_z + 1
    fit_props.num_z_fit     = fit_props.cheby_order_z + 1

    fit_props.generate_fitting_matrix( )

    fit_function            = lambda x, y, z: M1(x, y, z)[0].real

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_3d( fit_function, fit_props ) )

    # Start adaptive refinement loop
    fit_residual_3D_adaptive_interface( fit_function, ref_level )

    # Set starting position
    ref_level.set_start_index( 0 )
    ref_level.set_start_index_folded( 0 )


    return ref_level


def fit_residual_region_M1_dA_adaptive( ) -> None:
    fit_props               = FitProperties( )
    fit_props.dims          = 3
    fit_props.region_name   = "M1_dA"
    fit_props.cheby_order_x = 10
    fit_props.cheby_order_y = 10
    fit_props.cheby_order_z = 50
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = False
    fit_props.x_max         = 1.0
    fit_props.x_min         = 0.0
    fit_props.y_log_scale   = False
    fit_props.y_max         = 2.0
    fit_props.y_min         = 1.0
    fit_props.z_log_scale   = True
    fit_props.z_max         = 0.0
    fit_props.z_min         = -5.0
    fit_props.cheby_abs_tol = 5E-7
    fit_props.cheby_rel_tol = 1E-14
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.z_map_fcn     = lambda z: z
    fit_props.max_ref_level = 10

    fit_props.num_x         = fit_props.cheby_order_x + 1
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = fit_props.cheby_order_y + 1
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1
    fit_props.num_z         = fit_props.cheby_order_z + 1
    fit_props.num_z_fit     = fit_props.cheby_order_z + 1

    fit_props.generate_fitting_matrix( )

    fit_function            = lambda x, y, z: M1_dA(x, y, z)[0].real

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_3d( fit_function, fit_props ) )

    # Start adaptive refinement loop
    fit_residual_3D_adaptive_interface( fit_function, ref_level )

    # Set starting position
    ref_level.set_start_index( 0 )
    ref_level.set_start_index_folded( 0 )


    return ref_level


def fit_residual_region_M1_dB_adaptive( ) -> None:
    fit_props               = FitProperties( )
    fit_props.dims          = 3
    fit_props.region_name   = "M1_dB"
    fit_props.cheby_order_x = 10
    fit_props.cheby_order_y = 10
    fit_props.cheby_order_z = 50
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = False
    fit_props.x_max         = 1.0
    fit_props.x_min         = 0.0
    fit_props.y_log_scale   = False
    fit_props.y_max         = 2.0
    fit_props.y_min         = 1.0
    fit_props.z_log_scale   = True
    fit_props.z_max         = 0.0
    fit_props.z_min         = -5.0
    fit_props.cheby_abs_tol = 5E-7
    fit_props.cheby_rel_tol = 1E-14
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.z_map_fcn     = lambda z: z
    fit_props.max_ref_level = 10

    fit_props.num_x         = fit_props.cheby_order_x + 1
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = fit_props.cheby_order_y + 1
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1
    fit_props.num_z         = fit_props.cheby_order_z + 1
    fit_props.num_z_fit     = fit_props.cheby_order_z + 1

    fit_props.generate_fitting_matrix( )

    fit_function            = lambda x, y, z: M1_dB(x, y, z)[0].real

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_3d( fit_function, fit_props ) )

    # Start adaptive refinement loop
    fit_residual_3D_adaptive_interface( fit_function, ref_level )

    # Set starting position
    ref_level.set_start_index( 0 )
    ref_level.set_start_index_folded( 0 )


    return ref_level


def fit_residual_region_M2_adaptive( ) -> None:
    fit_props               = FitProperties( )
    fit_props.dims          = 1
    fit_props.region_name   = "M2"
    fit_props.cheby_order_x = 20
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = True
    fit_props.x_max         = 0.0
    fit_props.x_min         = -5.0
    fit_props.cheby_abs_tol = 5E-7
    fit_props.cheby_rel_tol = 1E-14
    fit_props.x_map_fcn     = lambda x: x
    fit_props.max_ref_level = 10

    fit_props.num_x         = fit_props.cheby_order_x + 1
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1

    fit_props.generate_fitting_matrix( )

    fit_function            = lambda x: M2( x )[0].real

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_1d( fit_function, fit_props ) )

    # Start adaptive refinement loop
    fit_residual_1D_adaptive_interface( fit_function, ref_level )

    # Set starting position
    ref_level.set_start_index( 0 )


    return ref_level


def fit_residual_region_M3_adaptive( ) -> None:
    fit_props               = FitProperties( )
    fit_props.dims          = 3
    fit_props.region_name   = "M3"
    fit_props.cheby_order_x = 10
    fit_props.cheby_order_y = 10
    fit_props.cheby_order_z = 50
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = False
    fit_props.x_max         = 1.0
    fit_props.x_min         = 0.0
    fit_props.y_log_scale   = False
    fit_props.y_max         = 2.0
    fit_props.y_min         = 1.0
    fit_props.z_log_scale   = True
    fit_props.z_max         = 3.0
    fit_props.z_min         = 0.0
    fit_props.cheby_abs_tol = 5E-7
    fit_props.cheby_rel_tol = 1E-14
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.z_map_fcn     = lambda z: z
    fit_props.max_ref_level = 10

    fit_props.num_x         = fit_props.cheby_order_x + 1
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = fit_props.cheby_order_y + 1
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1
    fit_props.num_z         = fit_props.cheby_order_z + 1
    fit_props.num_z_fit     = fit_props.cheby_order_z + 1

    fit_props.generate_fitting_matrix( )

    fit_function            = lambda x, y, z: M3(x, y, z)[0].real

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_3d( fit_function, fit_props ) )

    # Start adaptive refinement loop
    fit_residual_3D_adaptive_interface( fit_function, ref_level )

    # Set starting position
    ref_level.set_start_index( 0 )
    ref_level.set_start_index_folded( 0 )


    return ref_level


def fit_residual_region_M3_dA_adaptive( ) -> None:
    fit_props               = FitProperties( )
    fit_props.dims          = 3
    fit_props.region_name   = "M3_dA"
    fit_props.cheby_order_x = 10
    fit_props.cheby_order_y = 10
    fit_props.cheby_order_z = 50
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = False
    fit_props.x_max         = 1.0
    fit_props.x_min         = 0.0
    fit_props.y_log_scale   = False
    fit_props.y_max         = 2.0
    fit_props.y_min         = 1.0
    fit_props.z_log_scale   = True
    fit_props.z_max         = 3.0
    fit_props.z_min         = 0.0
    fit_props.cheby_abs_tol = 5E-6
    fit_props.cheby_rel_tol = 1E-14
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.z_map_fcn     = lambda z: z
    fit_props.max_ref_level = 10

    fit_props.num_x         = fit_props.cheby_order_x + 1
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = fit_props.cheby_order_y + 1
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1
    fit_props.num_z         = fit_props.cheby_order_z + 1
    fit_props.num_z_fit     = fit_props.cheby_order_z + 1

    fit_props.generate_fitting_matrix( )

    fit_function            = lambda x, y, z: M3_dA(x, y, z)[0].real

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_3d( fit_function, fit_props ) )

    # Start adaptive refinement loop
    fit_residual_3D_adaptive_interface( fit_function, ref_level )

    # Set starting position
    ref_level.set_start_index( 0 )
    ref_level.set_start_index_folded( 0 )


    return ref_level


def fit_residual_region_M3_dB_adaptive( ) -> None:
    fit_props               = FitProperties( )
    fit_props.dims          = 3
    fit_props.region_name   = "M3_dB"
    fit_props.cheby_order_x = 10
    fit_props.cheby_order_y = 10
    fit_props.cheby_order_z = 50
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = False
    fit_props.x_max         = 1.0
    fit_props.x_min         = 0.0
    fit_props.y_log_scale   = False
    fit_props.y_max         = 2.0
    fit_props.y_min         = 1.0
    fit_props.z_log_scale   = True
    fit_props.z_max         = 3.0
    fit_props.z_min         = 0.0
    fit_props.cheby_abs_tol = 5E-6
    fit_props.cheby_rel_tol = 1E-14
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.z_map_fcn     = lambda z: z
    fit_props.max_ref_level = 10

    fit_props.num_x         = fit_props.cheby_order_x + 1
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = fit_props.cheby_order_y + 1
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1
    fit_props.num_z         = fit_props.cheby_order_z + 1
    fit_props.num_z_fit     = fit_props.cheby_order_z + 1

    fit_props.generate_fitting_matrix( )

    fit_function            = lambda x, y, z: M3_dB(x, y, z)[0].real

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_3d( fit_function, fit_props ) )

    # Start adaptive refinement loop
    fit_residual_3D_adaptive_interface( fit_function, ref_level )

    # Set starting position
    ref_level.set_start_index( 0 )
    ref_level.set_start_index_folded( 0 )


    return ref_level


if __name__ == "__main__":
    import sys

    # folder_path = sys.argv[1]
    # option      = sys.argv[2]
    folder_path = r"D:\sergio\developments\SeaMotions\aux_data"
    option      = "L0"
    # fit_L1_Hfix()
    # fit_M1_Hfix()
    # fit_L1_Afix_Bfix()
    # fit_L1()
    # fit_L3()
    # fit_L3_Hfix()
    # fit_L3_Afix_Bfix()
    # fit_L1_P0()
    # fit_L1_P3()
    # fit_L2()
    # fit_L2_P1()
    # file_path = r"E:\sergio\developments\SeaMotions\aux_tools\test_coeffs.py"
    # fid = h5py.File(file_path, "w")
    # fid.create_dataset("C_filter", data=C_filter)
    # fid.create_dataset("NCX_filter", data=NCX_filter)
    # fid.create_dataset("NCY_filter", data=NCY_filter)
    # fid.create_dataset("NCZ_filter", data=NCZ_filter)
    # fid.close()
    # fit_M2()
    # fit_L3_dA()
    # fit_L3_dB()
    # fit_L3_P2()
    # folder_path = r"D:\sergio\developments\SeaMotions\aux_tools\0_databases\1_finite_water_depth\Prec_1Em6"
    # fit_L1( folder_path )
    # fit_L2( folder_path )
    # fit_L3( folder_path )
    # fit_M1( folder_path )
    # fit_M2( folder_path )
    # fit_M3( folder_path )
    # fit_L1_dA( folder_path )
    # fit_L1_dB( folder_path )
    # fit_L3_dA( folder_path )
    # fit_L3_dB( folder_path )
    # fit_M1_dA( folder_path )
    # fit_M1_dB( folder_path )
    # fit_M3_dA( folder_path )
    # fit_M3_dB( folder_path )

    if option == "L0":
        fit_L0( folder_path )

    elif option == "L0_dA":
        fit_L0_dA( folder_path )

    elif option == "L0_dB":
        fit_L0_dB( folder_path )

    elif option == "L1":
        fit_L1( folder_path )

    elif option == "L1_dA":
        fit_L1_dA( folder_path )

    elif option == "L1_dB":
        fit_L1_dB( folder_path )

    elif option == "L2":
        fit_L2( folder_path )

    elif option == "L3":
        fit_L3( folder_path )

    elif option == "L3_dA":
        fit_L3_dA( folder_path )

    elif option == "L3_dB":
        fit_L3_dB( folder_path )

    elif option == "Linf":
        fit_Linf( folder_path )

    elif option == "Linf_dA":
        fit_Linf_dA( folder_path )

    elif option == "Linf_dB":
        fit_Linf_dB( folder_path )

    elif option == "M1":
        fit_M1( folder_path )

    elif option == "M1_dA":
        fit_M1_dA( folder_path )

    elif option == "M1_dB":
        fit_M1_dB( folder_path )

    elif option == "M2":
        fit_M2( folder_path )

    elif option == "M3":
        fit_M3( folder_path )

    elif option == "M3_dA":
        fit_M3_dA( folder_path )

    elif option == "M3_dB":
        fit_M3_dB( folder_path )