
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
import os
from typing import Callable

# Import general usage scientific libraries
from numpy import linspace, ndarray

# Import local modules
from base_integrals import (fxy, fxy_dx, fxy_dy,
                            L1, L1_dA, L1_dB, L2, L3, L3_dA, L3_dB,
                            M1, M1_dA, M1_dB, M2, M3, M3_dA, M3_dB)


# Define global module variables
COL_WIDTH = 25


class Bounds:
    def __init__(self) -> None:
        self.a_max = 0.0
        self.a_min = 0.0
        self.b_max = 0.0
        self.b_min = 0.0
        self.h_max = 0.0
        self.h_min = 0.0
        self.num_a = 10
        self.num_b = 10
        self.num_h = 10


def create_reference_database( folder_path: str, mode: int )->None:
    # Create reference database for infinite water depth integrals
    if mode == 0 or mode == 1:
        # Create R11 reference data
        bnd_r11         = Bounds()
        bnd_r11.a_max   = 5.0
        bnd_r11.a_min   = -2.0
        bnd_r11.b_max   = 5.0
        bnd_r11.b_min   = -5.0
        bnd_r11.num_a   = 100
        bnd_r11.num_b   = 100

        f_def           = lambda x, y: fxy( x, y, only_int=True )
        f_dx_def        = lambda x, y: fxy_dx( x, y, only_int=True )
        write_reference_file_2d( folder_path, "R11", f_def, f_dx_def, bnd_r11 )

    # Create reference database for finite water depth integrals
    if mode == 0 or mode == 2:
        # Create L1 reference data
        bnd_l1          = Bounds()
        bnd_l1.a_max    = 1.0
        bnd_l1.a_min    = 0.0
        bnd_l1.b_max    = 1.0
        bnd_l1.b_min    = 0.0
        bnd_l1.h_max    = 0.0
        bnd_l1.h_min    = -5.0
        write_reference_file_3d( folder_path, "L1", L1, L1_dA, L1_dB, bnd_l1 )

        # Create L2 reference data
        bnd_l2          = Bounds()
        bnd_l2.h_max    = 0.0
        bnd_l2.h_min    = -5.0
        write_reference_file_1d( folder_path, "L2", L2, bnd_l2 )

        # Create L3 reference data
        bnd_l3          = Bounds()
        bnd_l3.a_max    = 1.0
        bnd_l3.a_min    = 0.0
        bnd_l3.b_max    = 1.0
        bnd_l3.b_min    = 0.0
        bnd_l3.h_max    = 3.0
        bnd_l3.h_min    = 0.0
        write_reference_file_3d( folder_path, "L3", L3, L3_dA, L3_dB, bnd_l3 )

        # Create M1 reference data
        bnd_m1          = Bounds()
        bnd_m1.a_max    = 1.0
        bnd_m1.a_min    = 0.0
        bnd_m1.b_max    = 2.0
        bnd_m1.b_min    = 1.0
        bnd_m1.h_max    = 0.0
        bnd_m1.h_min    = -5.0
        write_reference_file_3d( folder_path, "M1", M1, M1_dA, M1_dB, bnd_m1 )

        # Create M2 reference data
        bnd_m2          = Bounds()
        bnd_m2.h_max    = 0.0
        bnd_m2.h_min    = -5.0
        write_reference_file_1d(folder_path, "M2", M2, bnd_m2)

        # Create M3 reference data
        bnd_m3          = Bounds()
        bnd_m3.a_max    = 1.0
        bnd_m3.a_min    = 0.0
        bnd_m3.b_max    = 2.0
        bnd_m3.b_min    = 1.0
        bnd_m3.h_max    = 3.0
        bnd_m3.h_min    = 0.0
        write_reference_file_3d( folder_path, "M3", M3, M3_dA, M3_dB, bnd_m3 )


def write_reference_file_1d(
                                folder_path: str, 
                                file_name: str, 
                                f_def: Callable, 
                                bounds: Bounds
                            ) -> None:
    # Define database points
    h_values    = 10**linspace( bounds.h_min, bounds.h_max, bounds.num_h )

    # Compose total file path
    file_path   = os.path.join( folder_path, file_name + ".dat" )

    # Open file unit to write data
    with open(file_path, "w") as fid:
        # Write database definition points
        fid.writelines( f"Num.H: {bounds.num_h:d}\n" )
        fid.writelines( " ".join( f"{v:0.6f}" for v in h_values ) + "\n" )
        
        # Write header for the columns
        header_names    = [ "G" ]
        header_label    = "".join( hn.center( COL_WIDTH, " " ) for hn in header_names ) + "\n"
        fid.writelines( header_label )

        # Loop over bounds to write the reference data
        for k in range(bounds.num_h):
            vtw     =  [
                            f_def( h_values[k] )[0].real
                        ]
            vtws    =  [ f"{fi:^0.16f}" for fi in vtw ]
            fid.writelines( "".join( fis.center( COL_WIDTH, " " ) for fis in vtws ) + "\n" ) 


def write_reference_file_2d( 
                                folder_path: str, 
                                file_name: str, 
                                f_def: Callable,
                                f_dx_def: Callable, 
                                bounds: Bounds 
                            ) -> None:
    # Define database points
    x_values = 10**linspace( bounds.a_min, bounds.a_max, bounds.num_a )
    y_values = 10**linspace( bounds.b_min, bounds.b_max, bounds.num_b )

    # Compose total file path
    file_path = os.path.join( folder_path, file_name + ".dat" )
    
    # Open file unit to write data
    with open(file_path, "w") as fid:
        # Write database definition points
        fid.writelines( f"Num.X: {bounds.num_a:d}\n" )
        fid.writelines( " ".join( f"{v:0.6f}" for v in x_values ) + "\n" )
        fid.writelines( f"Num.Y: {bounds.num_b:d}\n" )
        fid.writelines( " ".join( f"{v:0.6f}" for v in y_values ) + "\n" )

        # Write header for the columns
        header_names    = [ "G", "G_dX" ]
        header_label    = "".join( hn.center( COL_WIDTH, " " ) for hn in header_names ) + "\n"
        fid.writelines( header_label )

        # Loop over bounds to write the reference data
        for i in range(bounds.num_a):
            for j in range(bounds.num_b):
                vtw = [
                            f_def( x_values[i], y_values[j] ),
                            f_dx_def(x_values[i], y_values[j] )
                        ]
                vtws = [ f"{fi:^0.16f}" for fi in vtw ]
                fid.writelines( "".join( fis.center( COL_WIDTH, " " ) for fis in vtws ) + "\n" )


def write_reference_file_3d( 
                                folder_path: str, 
                                file_name: str, 
                                f_def: Callable,
                                f_da_def: Callable, 
                                f_db_def: Callable, 
                                bounds: Bounds 
                            ) -> None:
    # Define database points
    a_values = linspace( bounds.a_min, bounds.a_max, bounds.num_a )
    b_values = linspace( bounds.b_min, bounds.b_max, bounds.num_b )
    h_values = 10**linspace( bounds.h_min, bounds.h_max, bounds.num_h )

    # Compose total file path
    file_path = os.path.join( folder_path, file_name + ".dat" )
    
    # Open file unit to write data
    with open(file_path, "w") as fid:
        # Write database definition points
        fid.writelines( f"Num.A: {bounds.num_a:d}\n" )
        fid.writelines( " ".join( f"{v:0.6f}" for v in a_values ) + "\n" )
        fid.writelines( f"Num.B: {bounds.num_b:d}\n" )
        fid.writelines( " ".join( f"{v:0.6f}" for v in b_values ) + "\n" )
        fid.writelines( f"Num.H: {bounds.num_h:d}\n" )
        fid.writelines( " ".join( f"{v:0.6f}" for v in h_values ) + "\n" )

        # Write header for the columns
        header_names    = [ "G", "G_dA", "G_dB" ]
        header_label    = "".join( hn.center( COL_WIDTH, " " ) for hn in header_names ) + "\n"
        fid.writelines( header_label )

        # Loop over bounds to write the reference data
        for i in range(bounds.num_a):
            for j in range(bounds.num_b):
                for k in range(bounds.num_h):
                    vtw = [
                                f_def(a_values[i], b_values[j], h_values[k])[0].real,
                                f_da_def(a_values[i], b_values[j], h_values[k])[0].real,
                                f_db_def(a_values[i], b_values[j], h_values[k])[0].real
                            ]
                    vtws = [ f"{fi:^0.16f}" for fi in vtw ]
                    fid.writelines( "".join( fis.center( COL_WIDTH, " " ) for fis in vtws ) + "\n" )


if __name__ == "__main__":
    this_path = os.path.dirname( os.path.abspath(__file__) )
    target_folder = os.path.join( os.path.dirname(this_path), "tests", "tests_data", "green_depth_tables" )
    create_reference_database(target_folder, 1)