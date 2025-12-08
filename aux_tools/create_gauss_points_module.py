
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

# Import general usage scientific libraries
import numpy as np
import scipy as sp


def main( fipath: str ) -> None:
    with open( fipath, "w" ) as fid:
        fid.writelines( "\n" )
        fid.writelines( "#pragma once" )
        fid.writelines( "\n" )
        fid.writelines( "\n" )

        fid.writelines( "template<std::size_t NDim, std::size_t Order>\n" )
        fid.writelines( "struct GaussPointsT\n" )
        fid.writelines( "{\n" )
        fid.writelines( "    static_assert( Order >= 1 && Order < 20 );\n" )
        fid.writelines( "    static_assert( NDim >= 1 && NDim <= 3 );\n" )
        fid.writelines( "};\n" )

        fid.writelines( "\n" )
        fid.writelines( "\n" )

        for ndim in range( 1, 3 ):
            for i in range( 1, 20 ):
                write_class( fid, i, ndim )


def write_class( fid, i, dim ):

    roots, weights = sp.special.roots_legendre(i)

    if dim == 2:
        R2DX, R2DY = np.meshgrid( roots, roots, indexing="ij" )
        W2DX, W2DY = np.meshgrid( weights, weights, indexing="ij" )
    elif dim == 3:
        R3DX, R3DY, R3DZ = np.meshgrid( roots, roots, roots, indexing="ij" )
        W3DX, W3DY, W3DZ = np.meshgrid( weights, weights, weights, indexing="ij" )

    fid.writelines( "template<>\n" )
    fid.writelines( f"struct GaussPointsT<{dim:d},{i:d}>\n" )
    fid.writelines( "{\n" )
    fid.writelines( "public:\n" )
    fid.writelines( f"    MEMALINGR static constexpr std::size_t    N                 = {i:d};\n" )
    fid.writelines( "\n" )

    if dim == 1:
        write_vector( fid, roots, "roots_x" )
        fid.writelines( "\n" )
        write_vector( fid, weights, "weights_x" )
        fid.writelines( "\n" )

    elif dim == 2:
        write_vector( fid, R2DX.flatten( ), "roots_x" )
        fid.writelines( "\n" )
        write_vector( fid, R2DY.flatten( ), "roots_y" )
        fid.writelines( "\n" )

        write_vector( fid, W2DX.flatten( ), "weights_x" )
        fid.writelines( "\n" )
        write_vector( fid, W2DY.flatten( ), "weights_y" )
        fid.writelines( "\n" )

    elif dim == 3:
        write_vector( fid, R3DX.flatten( ), "roots_x" )
        fid.writelines( "\n" )
        write_vector( fid, R3DY.flatten( ), "roots_y" )
        fid.writelines( "\n" )
        write_vector( fid, R3DZ.flatten( ), "roots_z" )
        fid.writelines( "\n" )

        write_vector( fid, W3DX.flatten( ), "weights_x" )
        fid.writelines( "\n" )
        write_vector( fid, W3DY.flatten( ), "weights_y" )
        fid.writelines( "\n" )
        write_vector( fid, W3DZ.flatten( ), "weights_z" )
        fid.writelines( "\n" )

    fid.writelines( "};" )
    fid.writelines( "\n" )
    fid.writelines( "\n" )


def write_vector( fid, data, field_name ):
    header_str = f"MEMALINGR static constexpr cusfloat       {field_name:s}[{data.shape[0]:d}]"
    fid.writelines( "    " + f"{header_str:<60}"  + "= {\n" )
    for i in range( data.shape[0] ):
        number_str = f"{data[i]:+0.12E}"
        if i < data.shape[0]-1:
            number_str += ",\n"
        else:
            number_str += "\n"
        fid.writelines( f"{"":>72}" + number_str )
    close_bracket = "};\n"
    fid.writelines( f"{"":>68}" + close_bracket )


if __name__ == "__main__":
    this_path = os.path.dirname( os.path.abspath( __file__ ) )
    file_path = os.path.join( this_path, "gauss_t.hpp" )

    main( file_path )