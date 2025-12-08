
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
import os
from typing import Callable

# Import general usage scientific libraries
from bokeh.plotting import figure, show, save, output_file
from bokeh.models import ColumnDataSource, LinearColorMapper, ColorBar, HoverTool
from bokeh.transform import transform
from bokeh.palettes import Viridis256
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
from numpy import array, linspace, log10, meshgrid, ndarray, sqrt, zeros
from numpy import abs as np_abs
from numpy import max as np_max
from numpy.linalg import solve as np_solve
from scipy.special import eval_chebyt, roots_chebyt

from fit_tools import RefLevel


class FitProperties:

    def __init__( self ) -> None:
        self.alpha_shift        = np.ndarray
        self.cheby_order_x      = 0
        self.cheby_order_y      = 0
        self.cheby_order_z      = 0
        self.cheby_abs_tol      = 0.0
        self.cheby_abs_tol_f    = 2
        self.cheby_rel_tol      = 0.0
        self.cheby_rel_tol_f    = 2
        self.dims               = 0
        self.fcn_log_scale      = False
        self.G0_cheby_intv      = [ ]
        self.G0_cheby_np        = 0
        self.G0_coeffs          = np.ndarray
        self.G0_mui             = np.ndarray
        self.G0_mui_log         = np.ndarray
        self.G0_norm_f          = 0.0
        self.max_ref_level      = 0
        self.num_x              = 7
        self.num_x_fit          = 7
        self.num_y              = 7
        self.num_y_fit          = 7
        self.num_z              = 50
        self.num_z_fit          = 50
        self.region_name        = ""
        self.sysmat_S_inv       = np.ndarray
        self.sysmat_U           = np.ndarray
        self.sysmat_VT          = np.ndarray
        self.vander_mat         = np.ndarray
        self.x_hpatch_np        = 1
        self.x_log_scale        = False
        self.x_max              = 0
        self.x_map_fcn          = lambda x: x
        self.x_min              = 0.0
        self.y_hpatch_np        = 1
        self.y_log_scale        = False
        self.y_max              = 0.0
        self.y_map_fcn          = lambda y: y
        self.y_min              = 0.0
        self.z_hpatch_np        = 1
        self.z_log_scale        = False
        self.z_max              = 0.0
        self.z_map_fcn          = lambda z: z
        self.z_min              = 0.0

    def fit_points_to_order(self)->None:
        self.num_x      = self.cheby_order_x
        self.num_x_fit  = self.cheby_order_x
        self.num_y      = self.cheby_order_y
        self.num_y_fit  = self.cheby_order_y
        self.num_z      = self.cheby_order_z
        self.num_z_fit  = self.cheby_order_z

    def generate_fitting_matrix( self ) -> None:
        if self.dims == 0:
            raise ValueError( "Problem dimensions not provided" )
        
        if self.dims == 1:
            # Get polynomial dimensions
            cheby_order_x       = self.cheby_order_x

            # Define parametric space for function fit
            num_x_fit           = self.num_x_fit
            x_fit_poly          = roots_chebyt( num_x_fit )[0]

            # Calculate polynomial fit coefficients
            V                   = np.polynomial.chebyshev.chebvander( x_fit_poly, cheby_order_x-1 )
            U, S, VT            = np.linalg.svd( V, full_matrices=False )
            S_inv               = np.diag( 1 / S )

            self.sysmat_U       = U
            self.sysmat_VT      = VT
            self.sysmat_S_inv   = S_inv

        elif self.dims == 2:
            # Get polynomial dimensions
            cheby_order_x       = self.cheby_order_x
            cheby_order_y       = self.cheby_order_y

            # Define parametric space for function fit
            num_x_fit           = self.num_x_fit
            num_y_fit           = self.num_y_fit
            x_fit_poly          = roots_chebyt( num_x_fit )[0]
            y_fit_poly          = roots_chebyt( num_y_fit )[0]
            Xfp,Yfp             = meshgrid( x_fit_poly, y_fit_poly, indexing="ij" )

            # Calculate polynomial fit coefficients
            X2                  = Xfp.ravel( )
            Y2                  = Yfp.ravel( )
            V                   = np.polynomial.chebyshev.chebvander2d( X2, Y2, [ cheby_order_x-1, cheby_order_y-1 ] )
            U, S, VT            = np.linalg.svd( V, full_matrices=False )
            S_inv               = np.diag( 1 / S )

            self.sysmat_U       = U
            self.sysmat_VT      = VT
            self.sysmat_S_inv   = S_inv

        elif self.dims == 3:
            # Get polynomial dimensions
            cheby_order_x       = self.cheby_order_x
            cheby_order_y       = self.cheby_order_y
            cheby_order_z       = self.cheby_order_z

            # Define parametric space for function fit
            num_x_fit           = self.num_x_fit
            num_y_fit           = self.num_y_fit
            num_z_fit           = self.num_z_fit
            x_fit_poly          = roots_chebyt( num_x_fit )[0]
            y_fit_poly          = roots_chebyt( num_y_fit )[0]
            z_fit_poly          = roots_chebyt( num_z_fit )[0]
            Xfp,Yfp,Zfp         = meshgrid( x_fit_poly, y_fit_poly, z_fit_poly, indexing="ij" )

            # Calculate polynomial fit coefficients
            X2                  = Xfp.ravel( )
            Y2                  = Yfp.ravel( )
            Z2                  = Zfp.ravel( )
            V                   = np.polynomial.chebyshev.chebvander3d( X2, Y2, Z2, [ cheby_order_x-1, cheby_order_y-1, cheby_order_z-1 ] )
            U, S, VT            = np.linalg.svd( V, full_matrices=False )
            S_inv               = np.diag( 1 / S )

            self.sysmat_U       = U
            self.sysmat_VT      = VT
            self.sysmat_S_inv   = S_inv

        else:
            raise ValueError( f"Problem dimensions: {self.dims} not available" )

    def x_map(self, x: ndarray)->ndarray:
        return self.x_map_fcn(self.x_map_lin(x))
        # return self.x_map_fcn(x)

    def x_map_lin(self, x: ndarray)->ndarray:
        xmap = (x+1)*(self.x_max-self.x_min)/2 + self.x_min
        if self.x_log_scale:
            xmap = 10**xmap
        
        return xmap

    def y_map(self, y: ndarray)->ndarray:
        # return self.y_map_fcn(self.y_map_lin(y))
        return self.y_map_fcn(y)
    
    def y_map_lin(self, y: ndarray)->ndarray:
        ymap = (y+1)*(self.y_max-self.y_min)/2 + self.y_min
        if self.y_log_scale:
            ymap = 10**ymap
        
        return ymap

    def z_map(self, z: ndarray)->ndarray:
        # return self.y_map_fcn(self.y_map_lin(y))
        return self.z_map_fcn(z)
    
    def z_map_lin(self, z: ndarray)->ndarray:
        zmap = (z+1)*(self.z_max-self.z_min)/2 + self.z_min
        if self.z_log_scale:
            zmap = 10**zmap

        return zmap

class FitStats:

    def __init__(self) -> None:
        self.abs_err_over_thr   = 0.0
        self.abs_err_tol        = 0.0
        self.ref_level          = 0
        self.max_abs_err        = 0.0
        self.max_cheby_order    = 0
        self.max_rel_err        = 0.0
        self.max_value          = 0.0
        self.mean_abs_err       = 0.0
        self.mean_rel_err       = 0.0
        self.min_abs_err        = 0.0
        self.min_rel_err        = 0.0
        self.num_coeffs         = 0
        self.rel_err_over_thr   = 0.0
        self.rel_err_tol        = 0.0


def eval_chebyshev_1d(x: ndarray, n: int, c: ndarray)->ndarray:
    sol = 0.0
    for i in range(n):
        sol += c[i]*eval_chebyt(i, x)

    return sol


def eval_chebyshev_1d_filter(x: ndarray, c: ndarray, ncx: ndarray)->ndarray:
    sol = zeros(x.shape)
    for i in range(c.shape[0]):
        sol += c[i]*eval_chebyt(ncx[i], x)

    return sol


def eval_chebyshev_2d(x: ndarray, y: ndarray, n: int, c: ndarray)->ndarray:
    sol = 0.0
    for i in range(n):
        for j in range(n):
            sol += c[i*n+j]*eval_chebyt(i, x)*eval_chebyt(j, y)

    return sol


def eval_chebyshev_2d_filter(x: ndarray, y: ndarray, c: ndarray, ncx: ndarray, ncy: ndarray)->ndarray:
    sol = 0.0
    for i in range(c.shape[0]):
        sol += c[i]*eval_chebyt(ncx[i], x)*eval_chebyt(ncy[i], y)

    return sol


def eval_chebyshev_3d_filter(x: ndarray, y: ndarray, z: ndarray, c: ndarray, ncx: ndarray, ncy: ndarray,
                            ncz: ndarray)->ndarray:
    sol = 0.0
    for i in range(c.shape[0]):
        sol += c[i]*eval_chebyt(ncx[i], x)*eval_chebyt(ncy[i], y)*eval_chebyt(ncz[i], z)

    return sol


def fit_chebyshev_1d(x: ndarray, f: ndarray, n: int)->ndarray:
    num_points = x.shape[0]
    A = zeros((num_points, n))
    for i in range(n):
            A[:, i] = eval_chebyt(i, x)

    At = A.T.dot(A)
    F = A.T.dot(f)
    C = np_solve(At, F)

    return C


def fit_chebyshev_1d_robust( x: ndarray, f: ndarray, nx: int )->ndarray:
    num_points  = x.shape[0]
    V           = zeros( ( num_points, nx ) )
    for i in range( nx ):
        V[:, i] = eval_chebyt(i, x)
    
    coeffs, _, _, _ = np.linalg.lstsq(V, f, rcond=None)

    return coeffs


def fit_chebyshev_2d(x: ndarray, y: ndarray, f: ndarray, nx: int, ny: int)->ndarray:
    num_points = x.shape[0]
    A = zeros((num_points, nx*ny))
    for i in range(nx):
        for j in range(ny):
            A[:, i*ny+j] = eval_chebyt(i, x)*eval_chebyt(j, y)

    At = A.T.dot(A)
    F = A.T.dot(f)
    C = np_solve(At, F)

    return C


def fit_chebyshev_2d_b(x: ndarray, y: ndarray, f: ndarray, nx: int, ny: int)->ndarray:
    num_points = x.shape[0]
    A = zeros((num_points, nx*ny))
    for i in range(nx):
        for j in range(ny):
            A[:, i*ny+j] = eval_chebyt(i, x)*eval_chebyt(j, y)

    At = A.T.dot(A)
    F = A.T.dot(f)
    C = np_solve(At, F)

    return C, A, At, F


def fit_chebyshev_2d_robust(x: ndarray, y: ndarray, f: ndarray, nx: int, ny: int)->ndarray:
    V               = np.polynomial.chebyshev.chebvander2d( x, y, [ nx, ny ] )
    coeffs, _, _, _ = np.linalg.lstsq(V, f, rcond=None)

    return coeffs


def fit_chebyshev_3d(x: ndarray, y: ndarray, z: ndarray, f: ndarray, nx: int, ny: int, nz: int)->ndarray:
    num_points = x.shape[0]
    A = zeros((num_points, nx*ny*nz))
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                A[:, i*(ny*nz)+j*nz+k] = eval_chebyt(i, x)*eval_chebyt(j, y)*eval_chebyt(k, z)

    At = A.T.dot(A)
    F = A.T.dot(f)
    C = np_solve(At, F)

    return C


def fit_chebyshev_3d_robust(x: ndarray, y: ndarray, z: ndarray, f: ndarray, nx: int, ny: int, nz: int)->ndarray:
    print( "nx: ", nx, " - ny: ", ny, " - nz: ", nz )
    print( "Generating vandermonde matrix..." )
    V               = np.polynomial.chebyshev.chebvander3d( x, y, z, [ nx, ny, nz ] )
    print( "Solving system matrix...", end="" )
    coeffs, _, _, _ = np.linalg.lstsq(V, f, rcond=None)
    # U, S, VT = np.linalg.svd(V, full_matrices=False)
    # S_inv = np.diag(1 / S)
    # coeffs = VT.T @ S_inv @ U.T @ f  # Chebyshev coefficients
    print( " -> Done" )

    return coeffs


def fit_integral_1d(
                        f_residual: Callable,
                        fit_props: FitProperties,
                        show_stats= False,
                        create_figs = False,
                        show_figs = False
                    )->None:
    # Define parametric space limits
    x_max               = 10**fit_props.x_max if fit_props.x_log_scale else fit_props.x_max
    x_min               = 10**fit_props.x_min if fit_props.x_log_scale else fit_props.x_min
    cheby_order_x       = fit_props.cheby_order_x
    cheby_abs_tol       = fit_props.cheby_abs_tol
    cheby_abs_tol_f     = fit_props.cheby_abs_tol_f
    cheby_rel_tol       = fit_props.cheby_rel_tol
    cheby_rel_tol_f     = fit_props.cheby_rel_tol_f

    # Define parametric space for function fit
    num_x_fit           = fit_props.num_x_fit
    x_fit_poly          = roots_chebyt( num_x_fit )[0]
    x_fit               = fit_props.x_map_lin( x_fit_poly )

    # Calculate residual function over the fit points
    Ff = zeros( (num_x_fit, ) )
    for i in range( num_x_fit ):
        Ff[i] = f_residual( x_fit[i] )

    # Calculate polynomial fit coefficients
    C = fit_props.sysmat_VT.T @ fit_props.sysmat_S_inv @ fit_props.sysmat_U.T @ Ff
    
    # Filter coefficients to the precision required
    n_cheby         = linspace( 0, cheby_order_x-1, cheby_order_x, dtype=int )
    pos             = ( np_abs(C) > np_max( np_abs(C)) * cheby_rel_tol / cheby_rel_tol_f ) & ( np_abs(C) > cheby_abs_tol / cheby_abs_tol_f )
    C_filter        = C[pos]
    NCX_filter      = n_cheby[pos]


    if C_filter.shape[0] == 0:
        C_filter    = array( [ 0.0 ] )
        NCX_filter  = array( [ 0 ], dtype=int )

    # Define parametric space
    num_x   = fit_props.num_x
    xe      = linspace( -1.0, 1.0, num_x )
    x       = fit_props.x_map_lin( xe )

    # Calculate residual function over the evaluation points
    Fr = zeros( (num_x, ) )
    for i in range( num_x ):
        Fr[i] = f_residual( x[i] )
    
    Ffe = eval_chebyshev_1d_filter( xe, C_filter, NCX_filter )
    
    # Calculate and storage fit statistics
    err_abs                     = np_abs( Ffe - Fr )
    err_rel                     = err_abs  / np_abs( Fr ).max( )
    pos_abs_above               = ( np.abs( Fr ) > np.log10( cheby_abs_tol ) ) if fit_props.fcn_log_scale else ( np.abs( Fr ) > cheby_abs_tol )
    pos_abs                     = ( err_abs > cheby_abs_tol ) & pos_abs_above
    pos_rel                     = ( err_rel > cheby_rel_tol ) & pos_abs_above
    err_abs_over_thr            = pos_abs.sum( ) / err_abs.shape[0] * 100
    err_rel_over_thr            = pos_rel.sum( ) / err_rel.shape[0] * 100
    fit_stats                   = FitStats( )
    fit_stats.abs_err_tol       = cheby_abs_tol
    fit_stats.rel_err_tol       = cheby_rel_tol
    fit_stats.max_abs_err       = err_abs.max( )
    fit_stats.max_cheby_order   = NCX_filter.max( )
    fit_stats.max_rel_err       = err_rel.max( )
    fit_stats.max_value         = Fr.max( )
    fit_stats.mean_abs_err      = err_abs.mean( )
    fit_stats.mean_rel_err      = err_rel.mean( )
    fit_stats.min_abs_err       = err_abs.min( )
    fit_stats.min_rel_err       = err_rel.min( )
    fit_stats.abs_err_over_thr  = err_abs_over_thr
    fit_stats.rel_err_over_thr  = err_rel_over_thr
    fit_stats.num_coeffs        = C_filter.shape[0]

    # Print error statistics
    region_name = f"{fit_props.region_name}: -> X: {x_min:0.3E} - {x_max:0.3E}"
    print( f"Error Statistics - {region_name}:", flush=True )
    print( f" -> Num.Coeffs: {fit_stats.num_coeffs:d}", flush=True )
    print( f" -> Absolute Error Statistics:", flush=True )
    print( f"     -> Maximum:   {fit_stats.max_abs_err}", flush=True )
    print( f"     -> Minimum:   {fit_stats.min_abs_err}", flush=True )
    print( f"     -> Mean:      {fit_stats.mean_abs_err}", flush=True )
    print( f"     -> V.O.T:     {fit_stats.abs_err_over_thr} % [Threshold: {cheby_abs_tol:0.1E}]", flush=True )
    print( f" -> Relative Error Statistics:", flush=True )
    print( f"     -> Maximum:   {fit_stats.max_rel_err}", flush=True )
    print( f"     -> Minimum:   {fit_stats.min_rel_err}", flush=True )
    print( f"     -> Mean:      {fit_stats.mean_rel_err}", flush=True )
    print( f"     -> V.O.T:     {fit_stats.rel_err_over_thr} % [Threshold: {cheby_rel_tol:0.1E}]", flush=True )

    # Plot residual function
    if show_figs:
        create_figs = True
    
    if create_figs:
        fig = plt.figure()
        fig.suptitle(region_name + f" - Total Coeffs: {C_filter.shape[0]:d}")
        ax0 = fig.add_subplot(221)
        ax1 = fig.add_subplot(223)
        ax2 = fig.add_subplot(222)
        ax3 = fig.add_subplot(224)

        ax0.plot(x, Fr)
        ax0.set_xlabel("X")
        if fit_props.x_log_scale:
            ax0.set_xscale("log")

        ax1.plot(x, Fr, label="Target")
        ax1.plot(x, Ffe, label="Fit")
        ax1.set_xlabel("X")
        if fit_props.x_log_scale:
            ax1.set_xscale("log")
        ax1.legend()

        ax2.set_title("log10(|Ffe-Fr|)")
        ax2.plot(x, log10(err_abs))
        ax2.set_xlabel("X")
        if fit_props.x_log_scale:
            ax2.set_xscale("log")

        ax3.set_title("log10(|Ffe-Fr|/|Fr|)")
        ax3.plot(x, log10(np_abs(Ffe-Fr)/np_abs(Fr)))
        ax3.set_xlabel("X")
        if fit_props.x_log_scale:
            ax3.set_xscale("log")

        cheby_order_ticks = linspace(0, cheby_order_x-1, cheby_order_x)
        fig = plt.figure()
        fig.suptitle("Fit coefficient matrix")
        ax = fig.add_subplot(111)
        ax.plot(n_cheby, C, "-x", label="Raw")
        ax.plot(NCX_filter, C_filter, "o", label="Filter")
        ax.set_xlabel("X")
        ax.set_xticks(cheby_order_ticks)

    if show_figs:
        plt.show()

    return np.vstack( ( C_filter, NCX_filter ) ).T, fit_stats


def fit_integral_2d(
                            f_residual: Callable,
                            fit_props: FitProperties,
                            show_figs = False,
                            stop_to_show = True,
                            ref_level=None
                        ) -> None:
    # Define parametric space limits
    num_cross_sections  = 5
    x_max               = 10**fit_props.x_max if fit_props.x_log_scale else fit_props.x_max
    x_min               = 10**fit_props.x_min if fit_props.x_log_scale else fit_props.x_min
    y_max               = 10**fit_props.y_max if fit_props.y_log_scale else fit_props.y_max
    y_min               = 10**fit_props.y_min if fit_props.y_log_scale else fit_props.y_min
    cheby_order_x       = fit_props.cheby_order_x
    cheby_order_y       = fit_props.cheby_order_y
    cheby_abs_tol       = fit_props.cheby_abs_tol
    cheby_abs_tol_f     = fit_props.cheby_abs_tol_f
    cheby_rel_tol       = fit_props.cheby_rel_tol
    cheby_rel_tol_f     = fit_props.cheby_rel_tol_f

    # Define parametric space for function fit
    num_x_fit           = fit_props.num_x_fit
    num_y_fit           = fit_props.num_y_fit
    x_fit_poly          = roots_chebyt( num_x_fit )[0]
    y_fit_poly          = roots_chebyt( num_y_fit )[0]
    x_fit               = fit_props.x_map_lin( x_fit_poly )
    y_fit               = fit_props.y_map_lin( y_fit_poly )

    # Calculate residual function over the fit points
    Ff = zeros((num_x_fit, num_y_fit))
    for i in range(num_x_fit):
        for j in range(num_y_fit):
            Ff[i, j] = f_residual(x_fit[i], y_fit[j])
    
    F2  = Ff.flatten( )
    C   = fit_props.sysmat_VT.T @ fit_props.sysmat_S_inv @ fit_props.sysmat_U.T @ F2

    # Filter coefficients to the precision required
    n_cheby_x       = linspace(0, cheby_order_x-1, cheby_order_x, dtype=int)
    n_cheby_y       = linspace(0, cheby_order_y-1, cheby_order_y, dtype=int)
    NCX,NCY         = meshgrid(n_cheby_x, n_cheby_y, indexing="ij")
    pos             = ( np_abs(C) > np_max( np_abs(C)) * cheby_rel_tol / cheby_rel_tol_f ) & ( np_abs(C) > cheby_abs_tol / cheby_abs_tol_f )
    C_filter        = C[pos]
    NCX_filter      = ( NCX.ravel( ) )[pos]
    NCY_filter      = ( NCY.ravel( ) )[pos]

    if C_filter.shape[0] == 0:
        C_filter    = array( [ 0.0 ] )
        NCX_filter  = array( [ 0 ], dtype=int )
        NCY_filter  = array( [ 0 ], dtype=int )

    # Define parametric space
    num_x   = fit_props.num_x
    num_y   = fit_props.num_y
    xe      = linspace( -1.0, 1.0, num_x )
    ye      = linspace( -1.0, 1.0, num_y )
    x       = fit_props.x_map_lin( xe )
    y       = fit_props.y_map_lin( ye )
    Xe,Ye   = meshgrid( xe, ye, indexing="ij" )
    X,Y     = meshgrid( fit_props.x_map_fcn( x ), fit_props.y_map_fcn( y ), indexing="ij" )

    # Calculate residual function over the evaluation points
    Fr = zeros( ( num_x, num_y ) )
    for i in range( num_x ):
        for j in range( num_y ):
            Fr[i, j] = f_residual( X[i, j], Y[i, j] )

    Ffe = eval_chebyshev_2d_filter( Xe, Ye, C_filter, NCX_filter, NCY_filter )
    Ffe = Ffe.reshape(num_x, num_y)

    y0 = linspace(fit_props.y_min, fit_props.y_max, num_cross_sections)
    if fit_props.y_log_scale:
        y0 = 10**y0

    FY0 = zeros((num_cross_sections, num_x))
    for i in range(num_cross_sections):
        for j in range(num_x):
            FY0[i, j] = f_residual(x[j], y0[i])

    x0 = linspace(fit_props.x_min, fit_props.x_max, num_cross_sections)
    if fit_props.x_log_scale:
        x0 = 10**x0

    FX0 = zeros((num_cross_sections, num_y))
    for i in range(num_cross_sections):
        for j in range(num_y):
            FX0[i, j] = f_residual(x0[i], y[j])
    


    # Calculate and storage fit statistics
    err_abs                     = np_abs( Ffe - Fr )
    err_rel                     = err_abs  / np_abs( Fr ).max( )
    pos_abs_above               = ( np.abs( Fr ) > np.log10( cheby_abs_tol ) ) if fit_props.fcn_log_scale else ( np.abs( Fr ) > cheby_abs_tol )
    pos_abs                     = ( err_abs > cheby_abs_tol ) & pos_abs_above
    pos_rel                     = ( err_rel > cheby_rel_tol ) & pos_abs_above
    err_abs_over_thr            = pos_abs.sum( ) / err_abs.shape[0] / err_abs.shape[1] * 100
    err_rel_over_thr            = pos_rel.sum( ) / err_rel.shape[0] / err_rel.shape[1] * 100
    fit_stats                   = FitStats( )
    fit_stats.abs_err_tol       = cheby_abs_tol
    fit_stats.rel_err_tol       = cheby_rel_tol
    fit_stats.max_abs_err       = err_abs.max( )
    fit_stats.max_cheby_order   = max( NCX_filter.max( ), NCY_filter.max( ) )
    fit_stats.max_rel_err       = err_rel.max( )
    fit_stats.max_value         = Fr.max( )
    fit_stats.mean_abs_err      = err_abs.mean( )
    fit_stats.mean_rel_err      = err_rel.mean( )
    fit_stats.min_abs_err       = err_abs.min( )
    fit_stats.min_rel_err       = err_rel.min( )
    fit_stats.abs_err_over_thr  = err_abs_over_thr
    fit_stats.rel_err_over_thr  = err_rel_over_thr
    fit_stats.num_coeffs        = C_filter.shape[0]

    # Print error statistics
    region_name = f"{fit_props.region_name}: -> X: {x_min:0.3E} - {x_max:0.3E} | Y: {y_min:0.3E} - {y_max:0.3E}"
    print( f"Error Statistics - {region_name}:", flush=True )
    print( f" -> Num.Coeffs: {fit_stats.num_coeffs:d}", flush=True )
    print( f" -> Absolute Error Statistics:", flush=True )
    print( f"     -> Maximum:   {fit_stats.max_abs_err}", flush=True )
    print( f"     -> Minimum:   {fit_stats.min_abs_err}", flush=True )
    print( f"     -> Mean:      {fit_stats.mean_abs_err}", flush=True )
    print( f"     -> V.O.T:     {fit_stats.abs_err_over_thr} % [Threshold: {cheby_abs_tol:0.1E}]", flush=True )
    print( f" -> Relative Error Statistics:", flush=True )
    print( f"     -> Maximum:   {fit_stats.max_rel_err}", flush=True )
    print( f"     -> Minimum:   {fit_stats.min_rel_err}", flush=True )
    print( f"     -> Mean:      {fit_stats.mean_rel_err}", flush=True )
    print( f"     -> V.O.T:     {fit_stats.rel_err_over_thr} % [Threshold: {cheby_rel_tol:0.1E}]", flush=True )

    # Plot residual function
    if ref_level is not None:
        if ref_level > 10:
            if show_figs or stop_to_show:
                fig = plt.figure()
                fig.suptitle(region_name + f" - Total Coeffs: {C_filter.shape[0]:d}")
                ax0 = fig.add_subplot(231)
                ax1 = fig.add_subplot(232, projection='3d')
                ax4 = fig.add_subplot(233)
                ax2 = fig.add_subplot(234)
                ax3 = fig.add_subplot(235)
                ax5 = fig.add_subplot(236)

                cfx = ax0.contourf(X, Y, Fr)
                ax0.set_xlabel("X")
                ax0.set_ylabel("Y")
                if fit_props.x_log_scale:
                    ax0.set_xscale("log")
                if fit_props.y_log_scale:
                    ax0.set_yscale("log")
                plt.colorbar(cfx, ax=ax0)

                psf = ax1.plot_surface(X, Y, Ffe, cmap="jet")
                ax1.set_xlabel("X")
                ax1.set_ylabel("Y")
                plt.colorbar(psf, ax=ax1)

                ax2.set_title(f"Y = {y_min:0.1E}")
                for i in range(num_cross_sections):
                    ax2.plot(fit_props.x_map_fcn(x), FY0[i, :], label=f"Y = {y0[i]:0.2E}")
                ax2.legend()
                if fit_props.x_log_scale:
                    ax2.set_xscale("log")
                ax2.set_xlabel("X")

                ax3.set_title(f"X = {x_min:0.1E}")
                for i in range(num_cross_sections):
                    ax3.plot(fit_props.y_map_fcn(y), FX0[i, :], label=f"X = {x0[i]:0.2E}")
                ax3.legend()
                if fit_props.y_log_scale:
                    ax3.set_xscale("log")
                ax3.set_xlabel("Y")

                ax4.set_title("log10(|Ffe-Fr|)")
                cntf = ax4.contourf(X, Y, log10(err_abs), cmap="jet")
                ax4.set_xlabel("X")
                ax4.set_ylabel("Y")
                plt.colorbar(cntf, ax=ax4)

                ax5.set_title("log10(|Ffe-Fr|/|Fr|)")
                cntf = ax5.contourf(X, Y, log10(np_abs(Ffe-Fr)/np_abs(Fr)))
                ax5.set_xlabel("X")
                ax5.set_ylabel("Y")
                plt.colorbar(cntf, ax=ax5)

                cheby_order_x_ticks = linspace(0, cheby_order_x-1, cheby_order_x)
                cheby_order_y_ticks = linspace(0, cheby_order_y-1, cheby_order_y)
                fig = plt.figure()
                fig.suptitle("Fit coefficient matrix")
                ax = fig.add_subplot(111)
                imc = ax.imshow(log10(C.reshape(cheby_order_x, cheby_order_y).T), origin="lower")
                ax.set_xlabel("X")
                ax.set_xticks(cheby_order_x_ticks)
                ax.set_ylabel("Y")
                ax.set_yticks(cheby_order_y_ticks)
                plt.colorbar(imc, ax=ax)
    
    if stop_to_show:
        plt.show()

    return np.vstack( ( C_filter, NCX_filter, NCY_filter ) ).T, fit_stats


def fit_integral_3d(
                        f_residual: Callable,
                        fit_props: FitProperties,
                    )->None:
    # Define parametric space limits
    num_cross_sections  = 5
    x_max               = 10**fit_props.x_max if fit_props.x_log_scale else fit_props.x_max
    x_min               = 10**fit_props.x_min if fit_props.x_log_scale else fit_props.x_min
    y_max               = 10**fit_props.y_max if fit_props.y_log_scale else fit_props.y_max
    y_min               = 10**fit_props.y_min if fit_props.y_log_scale else fit_props.y_min
    z_max               = 10**fit_props.z_max if fit_props.z_log_scale else fit_props.z_max
    z_min               = 10**fit_props.z_min if fit_props.z_log_scale else fit_props.z_min
    cheby_order_x       = fit_props.cheby_order_x
    cheby_order_y       = fit_props.cheby_order_y
    cheby_order_z       = fit_props.cheby_order_z
    cheby_abs_tol       = fit_props.cheby_abs_tol
    cheby_abs_tol_f     = fit_props.cheby_abs_tol_f
    cheby_rel_tol       = fit_props.cheby_rel_tol
    cheby_rel_tol_f     = fit_props.cheby_rel_tol_f

    # Define parametric space for function fit
    num_x_fit           = fit_props.num_x_fit
    num_y_fit           = fit_props.num_y_fit
    num_z_fit           = fit_props.num_z_fit
    x_fit_poly          = roots_chebyt( num_x_fit )[0]
    y_fit_poly          = roots_chebyt( num_y_fit )[0]
    z_fit_poly          = roots_chebyt( num_z_fit )[0]
    x_fit               = fit_props.x_map_lin( x_fit_poly )
    y_fit               = fit_props.y_map_lin( y_fit_poly )
    z_fit               = fit_props.z_map_lin( z_fit_poly )
    Xfp,Yfp,Zfp         = meshgrid( x_fit_poly, y_fit_poly, z_fit_poly, indexing="ij" )
    Xf,Yf,Zf            = meshgrid( x_fit, y_fit, z_fit, indexing="ij" )

    # Calculate residual function over the fit points
    Ff = zeros( ( num_x_fit, num_y_fit, num_z_fit ) )
    for i in range( num_x_fit ):
        for j in range( num_y_fit ):
            for k in range( num_z_fit ):
                Ff[i, j, k] = f_residual( x_fit[i], y_fit[j], z_fit[k] )

    # Calculate polynomial fit coefficients
    F2  = Ff.ravel( )
    C   = fit_props.sysmat_VT.T @ fit_props.sysmat_S_inv @ fit_props.sysmat_U.T @ F2
    
    # Filter coefficients to the precision required
    n_cheby_x       = linspace( 0, cheby_order_x-1, cheby_order_x, dtype=int )
    n_cheby_y       = linspace( 0, cheby_order_y-1, cheby_order_y, dtype=int )
    n_cheby_z       = linspace( 0, cheby_order_z-1, cheby_order_z, dtype=int )
    NCX,NCY,NCZ     = meshgrid( n_cheby_x, n_cheby_y, n_cheby_z, indexing="ij" )
    pos             = ( np_abs(C) > np_max( np_abs(C)) * cheby_rel_tol / cheby_rel_tol_f ) & ( np_abs(C) > cheby_abs_tol / cheby_abs_tol_f )
    C_filter        = C[pos]
    NCX_filter      = ( NCX.ravel( ) )[pos]
    NCY_filter      = ( NCY.ravel( ) )[pos]
    NCZ_filter      = ( NCZ.ravel( ) )[pos]

    if C_filter.shape[0] == 0:
        C_filter    = array( [ 0.0 ] )
        NCX_filter  = array( [ 0 ], dtype=int )
        NCY_filter  = array( [ 0 ], dtype=int )
        NCZ_filter  = array( [ 0 ], dtype=int )

    # Define parametric space
    num_x           = fit_props.num_x
    num_y           = fit_props.num_y
    num_z           = fit_props.num_z
    xe              = linspace( -1.0, 1.0, num_x )
    ye              = linspace( -1.0, 1.0, num_y )
    ze              = linspace( -1.0, 1.0, num_z )
    x               = fit_props.x_map_lin( xe )
    y               = fit_props.y_map_lin( ye )
    z               = fit_props.z_map_lin( ze )
    Xe,Ye,Ze        = meshgrid( xe, ye, ze, indexing="ij" )
    X,Y,Z           = meshgrid( fit_props.x_map_fcn( x ), fit_props.y_map_fcn( y ), fit_props.z_map_fcn( z ), indexing="ij" )
    X0              = x_min
    Y0              = y_min
    Z0              = z_min

    # Calculate residual function over the evaluation points
    Fr = zeros( ( num_x, num_y, num_z ) )
    for i in range( num_x ):
        for j in range( num_y ):
            for k in range( num_z ):
                Fr[i, j, k] = f_residual( x[i], y[j], z[k] )
    
    Ffe = eval_chebyshev_3d_filter( Xe.ravel(), Ye.ravel(), Ze.ravel(), C_filter, NCX_filter, NCY_filter, NCZ_filter )
    Ffe = Ffe.reshape( num_x, num_y, num_z )

    xe2 = np.array( [ 0.0 ] )
    ye2 = np.array( [ 0.0 ] )
    ze2 = np.array( [ ze[0] ] )
    FF1 = eval_chebyshev_3d_filter( xe2, ye2, ze2, C_filter, NCX_filter, NCY_filter, NCZ_filter )
    print( xe2, ye2, ze2, FF1 )

    # Calculate and storage fit statistics
    err_abs                     = np_abs( Ffe - Fr )
    err_rel                     = err_abs  / np_abs( Fr ).max( )
    pos_abs_above               = ( np.abs( Fr ) > np.log10( cheby_abs_tol ) ) if fit_props.fcn_log_scale else ( np.abs( Fr ) > cheby_abs_tol )
    pos_abs                     = ( err_abs > cheby_abs_tol ) & pos_abs_above
    pos_rel                     = ( err_rel > cheby_rel_tol ) & pos_abs_above
    err_abs_over_thr            = pos_abs.sum( ) / err_abs.shape[0] / err_abs.shape[1] / err_abs.shape[2] * 100
    err_rel_over_thr            = pos_rel.sum( ) / err_rel.shape[0] / err_rel.shape[1] / err_rel.shape[2] * 100
    fit_stats                   = FitStats( )
    fit_stats.abs_err_tol       = cheby_abs_tol
    fit_stats.rel_err_tol       = cheby_rel_tol
    fit_stats.max_abs_err       = err_abs.max( )
    fit_stats.max_cheby_order   = max( NCX_filter.max( ), NCY_filter.max( ), NCZ_filter.max( ) )
    fit_stats.max_rel_err       = err_rel.max( )
    fit_stats.max_value         = Fr.max( )
    fit_stats.mean_abs_err      = err_abs.mean( )
    fit_stats.mean_rel_err      = err_rel.mean( )
    fit_stats.min_abs_err       = err_abs.min( )
    fit_stats.min_rel_err       = err_rel.min( )
    fit_stats.abs_err_over_thr  = err_abs_over_thr
    fit_stats.rel_err_over_thr  = err_rel_over_thr
    fit_stats.num_coeffs        = C_filter.shape[0]

    # Print error statistics
    region_name = f"{fit_props.region_name}: -> X: {x_min:0.3E} - {x_max:0.3E} | Y: {y_min:0.3E} - {y_max:0.3E} | Z: {z_min:0.3E} - {z_max:0.3E}"
    print( f"Error Statistics - {region_name}:", flush=True )
    print( f" -> Num.Coeffs: {fit_stats.num_coeffs:d}", flush=True )
    print( f" -> Absolute Error Statistics:", flush=True )
    print( f"     -> Maximum:   {fit_stats.max_abs_err}", flush=True )
    print( f"     -> Minimum:   {fit_stats.min_abs_err}", flush=True )
    print( f"     -> Mean:      {fit_stats.mean_abs_err}", flush=True )
    print( f"     -> V.O.T:     {fit_stats.abs_err_over_thr} % [Threshold: {cheby_abs_tol:0.1E}]", flush=True )
    print( f" -> Relative Error Statistics:", flush=True )
    print( f"     -> Maximum:   {fit_stats.max_rel_err}", flush=True )
    print( f"     -> Minimum:   {fit_stats.min_rel_err}", flush=True )
    print( f"     -> Mean:      {fit_stats.mean_rel_err}", flush=True )
    print( f"     -> V.O.T:     {fit_stats.rel_err_over_thr} % [Threshold: {cheby_rel_tol:0.1E}]", flush=True )

    return np.vstack( ( C_filter, NCX_filter, NCY_filter, NCZ_filter ) ).T, fit_stats


def fit_residual_1D_adaptive_interface( fit_function, ref_level: RefLevel ) -> None:
    # Define boundary values
    x_new = np.array( [ 
                            ref_level.fit_props.x_min, 
                            ( ref_level.fit_props.x_min + ref_level.fit_props.x_max ) / 2.0,
                            ref_level.fit_props.x_max
                        ] )
    
    # Loop over refinements to create new refinement levels
    for i in range( x_new.shape[0] - 1 ):
        # Apply new sub-region boundaries
        fit_props       = copy.copy( ref_level.fit_props )
        fit_props.x_max = x_new[i+1]
        fit_props.x_min = x_new[i]

        # Create new refinement level
        ref_level_i                 = RefLevel( fit_props, parent=ref_level )

        # Check for refinement level
        if ref_level.level == fit_props.max_ref_level:
            print( "--> MAXIMUM REFINEMENT LEVEL REACHED" )
            print( f" - X: {fit_props.x_min:0.3E}, {fit_props.x_max:0.3E}" )
        
        # Fit function over the interval
        ref_level_i.add_data( *fit_integral_1d( fit_function, fit_props ) )

        print( f" -> Ref.Level:     ", ref_level_i.level, flush=True )
        print( f" -> Passed?:       ", ref_level_i.check_tolerances( ), flush=True )
        print( "\n" )

        # Check fit tolerances and refine if any
        if not ref_level_i.check_tolerances( ) and ref_level_i.level < ref_level_i.fit_props.max_ref_level:
            fit_residual_1D_adaptive_interface(fit_function, ref_level_i )
        
        # Add ith refinement level region to the children list
        ref_level.child.append( ref_level_i )


def fit_residual_2D_adaptive_interface(fit_function, ref_level: RefLevel, is_square_ref=True, show_figs=False) -> None:
    # Define boundary values
    x_new = np.array( [ 
                            ref_level.fit_props.x_min, 
                            ( ref_level.fit_props.x_min + ref_level.fit_props.x_max ) / 2.0,
                            ref_level.fit_props.x_max
                        ] )
    y_new = np.array( [ 
                            ref_level.fit_props.y_min, 
                            ( ref_level.fit_props.y_min + ref_level.fit_props.y_max ) / 2.0,
                            ref_level.fit_props.y_max
                        ] )
    
    if is_square_ref:
        # Loop over refinements to create new refinement levels
        for i in range( x_new.shape[0] - 1 ):
            for j in range( y_new.shape[0] - 1 ):
                # Apply new sub-region boundaries
                fit_props       = copy.copy( ref_level.fit_props )
                fit_props.x_max = x_new[i+1]
                fit_props.x_min = x_new[i]
                fit_props.y_max = y_new[j+1]
                fit_props.y_min = y_new[j]

                # Create new refinement level
                ref_level_i     = RefLevel( fit_props, parent=ref_level )

                # Check for refinement level
                if ref_level.level == fit_props.max_ref_level:
                    print( "--> MAXIMUM REFINEMENT LEVEL REACHED" )
                    print( f" - X: {fit_props.x_min:0.3E}, {fit_props.x_max:0.3E}" )
                    print( f" - Y: {fit_props.y_min:0.3E}, {fit_props.y_max:0.3E}" )

                # Fit function over the interval
                ref_level_i.add_data( *fit_integral_2d( fit_function, fit_props, show_figs=show_figs, ref_level=ref_level_i.level ) )

                print( f" -> Ref.Level:     ", ref_level_i.level, flush=True )
                print( f" -> Passed?:       ", ref_level_i.check_tolerances( ), flush=True )
                print( "\n" )

                # Check fit tolerances and refine if any
                if not ref_level_i.check_tolerances( ) and ref_level_i.level < ref_level_i.fit_props.max_ref_level:
                    fit_residual_2D_adaptive_interface(fit_function, ref_level_i, show_figs=show_figs )
                
                # Add ith refinement level region to the children list
                ref_level.child.append( ref_level_i )

    else:
        # Loop over refinements to create new refinement levels
        ref_levels_x = [ ]
        for i in range( x_new.shape[0] - 1 ):
            # Apply new sub-region boundaries
            fit_props       = copy.copy( ref_level.fit_props )
            fit_props.x_max = x_new[i+1]
            fit_props.x_min = x_new[i]

            # Create new refinement level
            ref_level_i     = RefLevel( fit_props, parent=ref_level )

            # Fit function over the interval
            coeffs, fit_stats           = fit_integral_2d( fit_function, fit_props, show_figs=show_figs )
            ref_level_i.cheby_coeffs    = coeffs
            ref_level_i.fit_stats       = fit_stats

            # # Check fit tolerances and refine if any
            # if not ref_level_i.check_tolerances( ):
            #     fit_residual_2D_adaptive_interface( fit_function, ref_level_i, show_figs=show_figs )
            
            # Add ith refinement level region to the children list
            ref_levels_x.append( ref_level_i )

        # Loop over refinements to create new refinement levels
        ref_levels_y = [ ]
        for j in range( y_new.shape[0] - 1 ):
            # Apply new sub-region boundaries
            fit_props       = copy.copy( ref_level.fit_props )
            fit_props.y_max = y_new[j+1]
            fit_props.y_min = y_new[j]

            # Create new refinement level
            ref_level_i     = RefLevel( fit_props, parent=ref_level )

            # Fit function over the interval
            coeffs, fit_stats           = fit_integral_2d( fit_function, fit_props, show_figs=show_figs )
            ref_level_i.cheby_coeffs    = coeffs
            ref_level_i.fit_stats       = fit_stats

            # # Check fit tolerances and refine if any
            # if not ref_level_i.check_tolerances( ):
            #     fit_residual_2D_adaptive_interface( fit_function, ref_level_i, show_figs=show_figs )
            
            # Add ith refinement level region to the children list
            ref_levels_y.append( ref_level_i )

        # Check which refinement level is more expensive
        num_coeffs_x = 0
        for rl in ref_levels_x:
            num_coeffs_x += np.array( rl.get_num_cheby_coeffs( ) ).sum( )
        
        num_coeffs_y = 0
        for rl in ref_levels_y:
            num_coeffs_y += np.array( rl.get_num_cheby_coeffs( ) ).sum( )
        
        if num_coeffs_x > num_coeffs_y:
            for refli in ref_levels_y:
                if not refli.check_tolerances( ):
                    fit_residual_2D_adaptive_interface( fit_function, refli, is_square_ref=is_square_ref, show_figs=show_figs )
            
            ref_level.child.extend( ref_levels_y )


        else:
            for refli in ref_levels_x:
                if not refli.check_tolerances( ):
                    fit_residual_2D_adaptive_interface( fit_function, refli, is_square_ref=is_square_ref, show_figs=show_figs )

            ref_level.child.extend( ref_levels_x )

def fit_residual_3D_adaptive_interface( fit_function, ref_level: RefLevel ) -> None:
    # Define boundary values
    x_new = np.array( [ 
                            ref_level.fit_props.x_min, 
                            ( ref_level.fit_props.x_min + ref_level.fit_props.x_max ) / 2.0,
                            ref_level.fit_props.x_max
                        ] )
    y_new = np.array( [ 
                            ref_level.fit_props.y_min, 
                            ( ref_level.fit_props.y_min + ref_level.fit_props.y_max ) / 2.0,
                            ref_level.fit_props.y_max
                        ] )
    z_new = np.array( [ 
                            ref_level.fit_props.z_min, 
                            ( ref_level.fit_props.z_min + ref_level.fit_props.z_max ) / 2.0,
                            ref_level.fit_props.z_max
                        ] )
    
    # Loop over refinements to create new refinement levels
    if ref_level.level < 3:
        for i in range( x_new.shape[0] - 1 ):
            for j in range( y_new.shape[0] - 1 ):
                for k in range( z_new.shape[0] - 1 ):
                    # Apply new sub-region boundaries
                    fit_props       = copy.copy( ref_level.fit_props )
                    fit_props.x_max = x_new[i+1]
                    fit_props.x_min = x_new[i]
                    fit_props.y_max = y_new[j+1]
                    fit_props.y_min = y_new[j]
                    fit_props.z_max = z_new[k+1]
                    fit_props.z_min = z_new[k]

                    # Create new refinement level
                    ref_level_i                 = RefLevel( fit_props, parent=ref_level )

                    # Check for refinement level
                    if ref_level.level == fit_props.max_ref_level:
                        print( "--> MAXIMUM REFINEMENT LEVEL REACHED" )
                        print( f" - X: {fit_props.x_min:0.3E}, {fit_props.x_max:0.3E}" )
                        print( f" - Y: {fit_props.y_min:0.3E}, {fit_props.y_max:0.3E}" )
                        print( f" - Z: {fit_props.z_min:0.3E}, {fit_props.z_max:0.3E}" )

                    # Fit function over the interval
                    ref_level_i.add_data( *fit_integral_3d( fit_function, fit_props ) )

                    print( f" -> Ref.Level:     ", ref_level_i.level, flush=True )
                    print( f" -> Passed?:       ", ref_level_i.check_tolerances( ), flush=True )
                    print( "\n" )

                    # Check fit tolerances and refine if any
                    if not ref_level_i.check_tolerances( ) and ref_level_i.level < ref_level_i.fit_props.max_ref_level:
                        fit_residual_3D_adaptive_interface(fit_function, ref_level_i )
                    
                    # Add ith refinement level region to the children list
                    ref_level.child.append( ref_level_i )

    else:
        for k in range( z_new.shape[0] - 1 ):
            # Apply new sub-region boundaries
            fit_props       = copy.copy( ref_level.fit_props )
            fit_props.x_max = ref_level.fit_props.x_max
            fit_props.x_min = ref_level.fit_props.x_min
            fit_props.y_max = ref_level.fit_props.y_max
            fit_props.y_min = ref_level.fit_props.y_min
            fit_props.z_max = z_new[k+1]
            fit_props.z_min = z_new[k]

            # Create new refinement level
            ref_level_i                 = RefLevel( fit_props, parent=ref_level )

            # Check for refinement level
            if ref_level.level == fit_props.max_ref_level:
                print( "--> MAXIMUM REFINEMENT LEVEL REACHED" )
                print( f" - X: {fit_props.x_min:0.3E}, {fit_props.x_max:0.3E}" )
                print( f" - Y: {fit_props.y_min:0.3E}, {fit_props.y_max:0.3E}" )
                print( f" - Z: {fit_props.z_min:0.3E}, {fit_props.z_max:0.3E}" )

            # Fit function over the interval
            ref_level_i.add_data( *fit_integral_3d( fit_function, fit_props ) )

            print( f" -> Ref.Level:     ", ref_level_i.level, flush=True )
            print( f" -> Passed?:       ", ref_level_i.check_tolerances( ), flush=True )
            print( "\n" )

            # Check fit tolerances and refine if any
            if not ref_level_i.check_tolerances( ) and ref_level_i.level < ref_level_i.fit_props.max_ref_level:
                fit_residual_3D_adaptive_interface(fit_function, ref_level_i )
            
            # Add ith refinement level region to the children list
            ref_level.child.append( ref_level_i )


def fit_residual_3D_adaptive_interface_random( fit_function, ref_level: RefLevel ) -> None:
    # Define boundary values
    x_new = np.array( [ 
                            ref_level.fit_props.x_min, 
                            ( ref_level.fit_props.x_min + ref_level.fit_props.x_max ) / 2.0,
                            ref_level.fit_props.x_max
                        ] )
    y_new = np.array( [ 
                            ref_level.fit_props.y_min, 
                            ( ref_level.fit_props.y_min + ref_level.fit_props.y_max ) / 2.0,
                            ref_level.fit_props.y_max
                        ] )
    z_new = np.array( [ 
                            ref_level.fit_props.z_min, 
                            ( ref_level.fit_props.z_min + ref_level.fit_props.z_max ) / 2.0,
                            ref_level.fit_props.z_max
                        ] )
    
    # Loop over refinements to create new refinement levels
    for i in range( x_new.shape[0] - 1 ):
        for j in range( y_new.shape[0] - 1 ):
            for k in range( z_new.shape[0] - 1 ):
                # Apply new sub-region boundaries
                fit_props       = copy.copy( ref_level.fit_props )
                fit_props.x_max = x_new[i+1]
                fit_props.x_min = x_new[i]
                fit_props.y_max = y_new[j+1]
                fit_props.y_min = y_new[j]
                fit_props.z_max = z_new[k+1]
                fit_props.z_min = z_new[k]

                # Create new refinement level
                ref_level_i                 = RefLevel( fit_props, parent=ref_level )

                # Check for refinement level
                if ref_level.level == fit_props.max_ref_level:
                    print( "--> MAXIMUM REFINEMENT LEVEL REACHED" )
                    print( f" - X: {fit_props.x_min:0.3E}, {fit_props.x_max:0.3E}" )
                    print( f" - Y: {fit_props.y_min:0.3E}, {fit_props.y_max:0.3E}" )
                    print( f" - Z: {fit_props.z_min:0.3E}, {fit_props.z_max:0.3E}" )

                # Fit function over the interval
                ref_level_i.add_data( *fit_integral_3d( fit_function, fit_props ) )

                print( f" -> Ref.Level:     ", ref_level_i.level, flush=True )
                print( f" -> Passed?:       ", ref_level_i.check_tolerances( ), flush=True )
                print( "\n" )

                # Check fit tolerances and refine if any
                if np.random.rand( ) > 0.5 and ref_level_i.level < 5:
                # if not ref_level_i.check_tolerances( ) and ref_level_i.level < 1:
                    fit_residual_3D_adaptive_interface(fit_function, ref_level_i )
                
                # Add ith refinement level region to the children list
                ref_level.child.append( ref_level_i )
