
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
from matplotlib import cm
import matplotlib.pyplot as plt
from numpy import linspace, meshgrid, ndarray, ones

# Import local modules
from hydlib.math.integration.gauss import gauss_points


def plot_gauss_points(  ) -> None:
    N = 7

    aux_vec     = ones( ( N, ) )
    gpr, gpw    = gauss_points( N )

    gprx,gpry   = meshgrid( gpr, gpr )
    xi1         = ( gprx + 1.0 ) * ( 1.0 - gpry ) / 2.0 - 1.0
    xi2         = gpry

    fig = plt.figure( figsize=( 15, 15 ) )
    ax0 = fig.add_subplot( 211 )
    ax1 = fig.add_subplot( 212 )

    ax0.scatter( gprx.ravel( ), gpry.ravel( ) )
    ax1.scatter( xi1.ravel( ), xi2.ravel( ) )

    plt.show( )


def plot_shape_fcns(  ) -> None:
    # Generate mesh for calculation
    x   = linspace( -0.99, 0.99, 100 )
    y   = linspace( -0.99, 0.99, 100 )
    X,Y = meshgrid( x, y )

    plot_gauss_points( )

    # Plot three node shape functions
    plot_triangle_shape_fcn( X, Y )

    # Plot triangle shapes using karniadakis expansion
    plot_triangle_karniadakis( X, Y )

    # Plot four node shape functions
    plot_rectangle_shape_fcn( X, Y )


def plot_rectangle_shape_fcn( X: ndarray, Y: ndarray ) -> None:
    N0 = ( 1 - X ) / 2.0 * ( 1 - Y ) / 2.0
    N1 = ( 1 - X ) / 2.0 * ( 1 + Y ) / 2.0
    N2 = ( 1 + X ) / 2.0 * ( 1 - Y ) / 2.0
    N3 = ( 1 + X ) / 2.0 * ( 1 + Y ) / 2.0

    fig = plt.figure( figsize=( 15, 15 ) )
    ax0 = fig.add_subplot( 221, projection="3d" )
    ax1 = fig.add_subplot( 222, projection="3d" )
    ax2 = fig.add_subplot( 223, projection="3d" )
    ax3 = fig.add_subplot( 224, projection="3d" )

    _plot_surf( ax0, X, Y, N0, "N0" )
    _plot_surf( ax1, X, Y, N1, "N1" )
    _plot_surf( ax2, X, Y, N2, "N2" )
    _plot_surf( ax3, X, Y, N3, "N3" )

    plt.show( )


def plot_triangle_karniadakis( Xi1: ndarray, Xi2: ndarray ) -> None:
    # Calculate shape functions value
    N0 = ( 1 - Xi1 ) * ( 1 - Xi2 ) / 4.0
    N1 = ( 1 + Xi1 ) * ( 1 - Xi2 ) / 4.0
    N2 = ( 1 - Xi1 ) * ( 1 + Xi2 ) / 4.0
    N2 += ( 1 + Xi1 ) * ( 1 + Xi2 ) / 4.0

    X = -N0 + N1 - N2
    Y = -N0 - N1 + N2

    fig = plt.figure( figsize=( 15, 15 ) )
    ax0 = fig.add_subplot( 221, projection="3d" )
    ax1 = fig.add_subplot( 222, projection="3d" )
    ax2 = fig.add_subplot( 223, projection="3d" )
    ax3 = fig.add_subplot( 224, projection="3d" )

    _plot_surf( ax0, X, Y, N0, "N0" )
    _plot_surf( ax1, X, Y, N1, "N1" )
    _plot_surf( ax2, X, Y, N2, "N2" )
    plt.show( )


def plot_triangle_shape_fcn( X: ndarray, Y: ndarray ) -> None:
    Xm = ( 1 + X ) * ( 1 - Y ) / 2.0 - 1.0
    Ym = Y
    N0 = ( 1 - Xm ) / 2.0 * ( 1 - Ym ) / 2.0
    N1 = ( 1 - Xm ) / 2.0 * ( 1 + Ym ) / 2.0
    N2 = ( 1 + Xm ) / 2.0 * ( 1 - Ym ) / 2.0

    fig = plt.figure( figsize=( 15, 15 ) )
    ax0 = fig.add_subplot( 221, projection="3d" )
    ax1 = fig.add_subplot( 222, projection="3d" )
    ax2 = fig.add_subplot( 223, projection="3d" )
    ax3 = fig.add_subplot( 224, projection="3d" )

    _plot_surf( ax0, Xm, Ym, N0, "N0" )
    _plot_surf( ax1, Xm, Ym, N1, "N1" )
    _plot_surf( ax2, Xm, Ym, N2, "N2" )

    plt.show( )


def _plot_surf( ax, Xi: ndarray, Yi: ndarray, Zi: ndarray, title: str ):
    surf_i = ax.plot_surface( Xi, Yi, Zi, cmap=cm.jet )
    ax.set_title( title )
    plt.colorbar( surf_i )


if __name__ == "__main__":
    plot_shape_fcns( )