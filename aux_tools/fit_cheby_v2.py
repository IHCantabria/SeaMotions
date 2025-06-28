
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
        self.cheby_order_x  = 0
        self.cheby_order_y  = 0
        self.cheby_order_z  = 0
        self.cheby_abs_tol  = 0.0
        self.dims           = 0
        self.fcn_log_scale  = False
        self.max_ref_level  = 0
        self.cheby_rel_tol  = 0.0
        self.num_x          = 20
        self.num_x_fit      = 20
        self.num_y          = 20
        self.num_y_fit      = 20
        self.num_z          = 50
        self.num_z_fit      = 50
        self.region_name    = ""
        self.sysmat_S_inv   = np.ndarray
        self.sysmat_U       = np.ndarray
        self.sysmat_VT      = np.ndarray
        self.vander_mat     = np.ndarray
        self.x_log_scale    = False
        self.x_max          = 0
        self.x_map_fcn      = lambda x: x
        self.x_min          = 0.0
        self.y_log_scale    = False
        self.y_max          = 0.0
        self.y_map_fcn      = lambda y: y
        self.y_min          = 0.0
        self.z_log_scale    = False
        self.z_max          = 0.0
        self.z_map_fcn      = lambda z: z
        self.z_min          = 0.0

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


# class StatPatch:

#     def __init__( self, name: str ) -> None:
#         self.dx     = [ ]
#         self.dy     = [ ]
#         self.dz     = [ ]
#         self.name   = name
#         self.value  = [ ]
#         self.x      = [ ]
#         self.y      = [ ]


# class RefLevel:

#     def __init__( self, fit_props: FitProperties, parent=None ) -> None:
#         self.cheby_coeffs           = np.ndarray
#         self.child                  = [ ]
#         self.fit_props              = fit_props
#         self.fit_stats              = FitProperties( )
#         self.level                  = 0
#         self.parent                 = parent
#         self.start_index_global     = 0

#         if parent is not None:
#             self.level = parent.level + 1

#     def add_data( self, coeffs: np.ndarray, fit_stats: FitStats ) -> None:
#         # Storage data
#         self.cheby_coeffs           = coeffs
#         self.fit_stats              = fit_stats

#         # Set level value to fit_stats
#         self.fit_stats.ref_level    = self.level

#     def check_position_interval( self, pos: np.ndarray ):
#         if pos.shape[0] == 1:
#             is_pass =   ( 
#                             ( self.fit_props.x_max > pos[0] ) & ( self.fit_props.x_min < pos[0] ) 
#                         )
#         elif pos.shape[0] == 2:
#             is_pass =   ( 
#                             ( self.fit_props.x_max > pos[0] ) & ( self.fit_props.x_min < pos[0] )
#                             &
#                             ( self.fit_props.y_max > pos[1] ) & ( self.fit_props.y_min < pos[1] )
#                         )
            
#         elif pos.shape[0] == 3:
#             is_pass =   ( 
#                             ( self.fit_props.x_max > pos[0] ) & ( self.fit_props.x_min < pos[0] )
#                             &
#                             ( self.fit_props.y_max > pos[1] ) & ( self.fit_props.y_min < pos[1] )
#                             &
#                             ( self.fit_props.z_max > pos[2] ) & ( self.fit_props.z_min < pos[2] )
#                         )
#         else:
#             raise ValueError( "Could not check position in more that 3 dimensions." )
        
#         return is_pass

#     def check_tolerances( self ) -> bool:
#         is_pass_abs         = self.fit_stats.max_abs_err < self.fit_stats.abs_err_tol
#         is_pass_rel         = self.fit_stats.max_rel_err < self.fit_stats.rel_err_tol
#         # is_pass_thr_abs     = self.fit_stats.abs_err_over_thr < 1.0
#         is_pass_mean_abs    = self.fit_stats.mean_abs_err < self.fit_stats.abs_err_tol
#         is_pass_loose_abs   = self.fit_stats.max_abs_err < 10 * self.fit_stats.abs_err_tol
#         is_pass             = is_pass_abs or is_pass_rel or ( is_pass_mean_abs and is_pass_loose_abs )

#         return is_pass
    
#     def get_cheby_coeffs( self ) -> np.ndarray:
#         if self.child:
#             coeffs_i = [ ]
#             for c in self.child:
#                 coeffs_i.append( c.get_cheby_coeffs( ) )

#             coeffs_i = np.vstack( coeffs_i )

#         else:
#             coeffs_i = self.cheby_coeffs.copy( )

#         return coeffs_i
    
#     def get_max_level( self ) -> int:
#         if self.child:
#             ref_level = -1
#             for c in self.child:
#                 li = c.get_max_level( )
#                 if li > ref_level:
#                     ref_level = li
#         else:
#             ref_level = self.level

#         return ref_level

#     def get_num_cheby_coeffs( self ) -> np.ndarray:
#         if self.child:
#             cum_coeffs  = [ ]
#             last        = 0
#             for c in self.child:
#                 for cc in c.get_num_cheby_coeffs( ):
#                     last += cc
#                     cum_coeffs.append( cc )

#         else:
#             cum_coeffs = [ self.cheby_coeffs.shape[0] ]

#         return cum_coeffs

#     def get_start_index( self, pos: np.ndarray ) -> int:
#         if self.child:
#             for c in self.child:
#                 index = c.get_start_index( pos )
#                 if index is not None:
#                     break

#         else:
#             is_pos_interval = self.check_position_interval( pos )
#             if is_pos_interval:
#                 index = (self.start_index_global, self.cheby_coeffs.shape[0], self.fit_props.x_min, self.fit_props.x_max, self.fit_props.y_min, self.fit_props.y_max)
#                 print( "Region -> ", index )
            
#             else:
#                 index = None

#         return index
    
#     def get_stat_patch( self, stat_patch: StatPatch ) -> None:
#         if self.child:
#             for c in self.child:
#                 c.get_stat_patch( stat_patch )

#         else:
#             stat_patch.x.append( ( self.fit_props.x_max + self.fit_props.x_min ) / 2.0 )
#             stat_patch.y.append( ( self.fit_props.y_max + self.fit_props.y_min ) / 2.0 )
#             stat_patch.dx.append( self.fit_props.x_max - self.fit_props.x_min )
#             stat_patch.dy.append( self.fit_props.y_max - self.fit_props.y_min )
#             stat_patch.dz.append( self.fit_props.z_max - self.fit_props.z_min )
#             stat_patch.value.append( getattr( self.fit_stats, stat_patch.name ) )

#     def set_start_index( self, start_index: int ) -> int:
#         if self.child:
#             for c in self.child:
#                 start_index = c.set_start_index( start_index )
        
#         else:
#             self.start_index_global     = start_index + 0
#             start_index                 += self.cheby_coeffs.shape[0]
#             print( "Ref.Level: ", self.level, " - x_min: ", self.fit_props.x_min, " - x_max: ", self.fit_props.x_max, " - y_min: ", self.fit_props.y_min, " - y_max: ", self.fit_props.y_max, " - start_index: ", self.start_index_global )

#         return start_index
    
#     def show_summary( self, folder_path: str ) -> None:
#         # Define output coefficients
#         num_coeffs_sp       = StatPatch( "num_coeffs" )
#         max_err_sp          = StatPatch( "max_abs_err" )
#         ref_level_sp        = StatPatch( "ref_level" )
#         max_cheby_order_sp  = StatPatch( "max_cheby_order" )
#         max_abs_err_sp      = StatPatch( "max_abs_err" )
#         mean_abs_err_sp     = StatPatch( "mean_abs_err" )
#         min_abs_err_sp      = StatPatch( "min_abs_err" )

#         self.get_stat_patch( num_coeffs_sp )
#         self.get_stat_patch( max_err_sp )
#         self.get_stat_patch( ref_level_sp )
#         self.get_stat_patch( max_cheby_order_sp )
#         self.get_stat_patch( max_abs_err_sp )
#         self.get_stat_patch( mean_abs_err_sp )
#         self.get_stat_patch( min_abs_err_sp )

#         # Ouput summary data
#         print( "Num.Coeffs: ", np.array( num_coeffs_sp.value ).sum( ) )
#         print( "Num.Patches: ", len( max_err_sp.value ) )
#         print( "Max.Cheyby Order: ", np.array( max_cheby_order_sp.value ).max( ) )
#         with open( os.path.join( folder_path, "summary_stats.dat" ), "w" ) as fid:
#             fid.writelines( f"Num.Coeffs:            {np.array( num_coeffs_sp.value ).sum( )}\n"  )
#             fid.writelines( f"Num.Patches:           {len( max_err_sp.value )}\n" )
#             fid.writelines( f"Max.Cheyby Order:      {np.array( max_cheby_order_sp.value ).max( )}" )

#         plot_patch( num_coeffs_sp, folder_path )
#         plot_patch( max_err_sp, folder_path )
#         plot_patch( ref_level_sp, folder_path )
#         plot_patch( max_cheby_order_sp, folder_path )
#         plot_patch( max_abs_err_sp, folder_path )
#         plot_patch( mean_abs_err_sp, folder_path )
#         plot_patch( min_abs_err_sp, folder_path )


# def plot_patch( stat_patch: StatPatch, folder_path: str ) -> None:
#     # Create output file
#     fipath = os.path.join( folder_path, stat_patch.name + ".html" )
#     output_file( fipath, stat_patch.name.upper( ) )

#     # Sample data
#     x       = stat_patch.x
#     y       = stat_patch.y
#     width   = stat_patch.dx
#     height  = stat_patch.dy
#     value   = stat_patch.value

#     # Prepare main rectangle source
#     rect_data   = dict(x=x, y=y, width=width, height=height, value=value)
#     rect_source = ColumnDataSource(rect_data)

#     # Compute corner coordinates for each rectangle
#     corner_xs = []
#     corner_ys = []

#     for xi, yi, wi, hi in zip(x, y, width, height):
#         half_w = wi / 2
#         half_h = hi / 2
#         corners = [
#             (xi - half_w, yi - half_h),
#             (xi - half_w, yi + half_h),
#             (xi + half_w, yi - half_h),
#             (xi + half_w, yi + half_h),
#         ]
#         for cx, cy in corners:
#             corner_xs.append(cx)
#             corner_ys.append(cy)

#     corner_source = ColumnDataSource(data=dict(x=corner_xs, y=corner_ys))

#     # Set up color mapper
#     color_mapper = LinearColorMapper(palette=Viridis256, low=min(value), high=max(value))

#     # Create plot
#     p = figure(
#                     title=stat_patch.name.upper( ),
#                     x_axis_label="X",
#                     y_axis_label="Y",
#                     tools="pan,box_zoom,wheel_zoom,reset",
#                     sizing_mode="scale_height"
#                 )

#     # Draw rectangles with hover
#     rects = p.rect(x='x', y='y', width='width', height='height', 
#                 fill_color=transform('value', color_mapper), line_color='black', source=rect_source)

#     # Add hover tool for rectangles
#     hover = HoverTool(renderers=[rects],
#                     tooltips=[
#                         ("Value", "@value"),
#                         ("Center (x, y)", "(@x, @y)")
#                     ])
#     p.add_tools(hover)

#     # Draw corner circles
#     p.circle(x='x', y='y', size=6, color='black', source=corner_source)

#     # Add color bar
#     color_bar = ColorBar(color_mapper=color_mapper, label_standoff=12, location=(0, 0))
#     p.add_layout(color_bar, 'right')

#     save( p )
#     show( p )


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
    cheby_rel_tol       = fit_props.cheby_rel_tol

    # Define parametric space for function fit
    num_x_fit           = fit_props.num_x_fit
    x_fit_poly          = roots_chebyt( num_x_fit )[0]
    x_fit               = fit_props.x_map_lin( x_fit_poly )

    # Calculate residual function over the fit points
    Ff = zeros( (num_x_fit, ) )
    for i in range( num_x_fit ):
        Ff[i] = f_residual( x_fit[i] )

    # Calculate polynomial fit coefficients
    # C = fit_chebyshev_1d_robust( x_fit_poly, Ff, cheby_order_x )
    C = fit_props.sysmat_VT.T @ fit_props.sysmat_S_inv @ fit_props.sysmat_U.T @ Ff
    
    # Filter coefficients to the precision required
    n_cheby         = linspace( 0, cheby_order_x-1, cheby_order_x, dtype=int )
    pos             = ( np_abs(C) > np_max( np_abs(C)) * cheby_rel_tol ) & ( np_abs(C) > cheby_abs_tol )
    C_filter        = C[pos]
    NCX_filter      = n_cheby[pos]


    if C_filter.shape[0] == 0:
        C_filter    = array( [ 0.0 ] )
        NCX_filter  = array( [ 0 ], dtype=int )

    # Define parametric space
    num_x   = fit_props.num_x
    xe      = linspace( -1.0, 1.0, num_x )
    x       = fit_props.x_map_lin( xe )
    X0      = x_min

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
                            stop_to_show = False
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
    cheby_rel_tol       = fit_props.cheby_rel_tol

    # Define parametric space for function fit
    num_x_fit           = fit_props.num_x_fit
    num_y_fit           = fit_props.num_y_fit
    x_fit_poly          = roots_chebyt( num_x_fit )[0]
    y_fit_poly          = roots_chebyt( num_y_fit )[0]
    x_fit               = fit_props.x_map_lin( x_fit_poly )
    y_fit               = fit_props.y_map_lin( y_fit_poly )
    Xfp,Yfp             = meshgrid( x_fit_poly, y_fit_poly, indexing="ij" )
    Xf,Yf               = meshgrid( x_fit, y_fit, indexing="ij" )

    # Calculate residual function over the fit points
    Ff = zeros((num_x_fit, num_y_fit))
    for i in range(num_x_fit):
        for j in range(num_y_fit):
            Ff[i, j] = f_residual(x_fit[i], y_fit[j])
    
    # X2  = Xfp.flatten( )
    # Y2  = Yfp.flatten( )
    F2  = Ff.flatten( )
    # C   = fit_chebyshev_2d_robust( X2, Y2, F2, cheby_order_x-1, cheby_order_y-1 )
    C = fit_props.sysmat_VT.T @ fit_props.sysmat_S_inv @ fit_props.sysmat_U.T @ F2

    # Filter coefficients to the precision required
    n_cheby_x       = linspace(0, cheby_order_x-1, cheby_order_x, dtype=int)
    n_cheby_y       = linspace(0, cheby_order_y-1, cheby_order_y, dtype=int)
    NCX,NCY         = meshgrid(n_cheby_x, n_cheby_y, indexing="ij")
    pos             = ( np_abs(C) > np_max( np_abs(C)) * cheby_rel_tol ) & ( np_abs(C) > cheby_abs_tol )
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
    X0      = x_min
    Y0      = y_min

    # Calculate residual function over the evaluation points
    Fr = zeros( ( num_x, num_y ) )
    for i in range( num_x ):
        for j in range( num_y ):
            Fr[i, j] = f_residual( x[i], y[j] )

    
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

    # print( "C: ", C )
    # print( "C_filter: ", C_filter )
    # print( "NCX: ", NCX_filter )
    # print( "NCY: ", NCY_filter )
    # fig = plt.figure( )
    # ax0 = fig.add_subplot( 251 )
    # ax1 = fig.add_subplot( 252, projection="3d" )
    # ax2 = fig.add_subplot( 253, projection="3d" )
    # ax3 = fig.add_subplot( 254 )
    # ax4 = fig.add_subplot( 255 )

    # cnf = ax0.contourf( Xe, Ye, Fr )
    # plt.colorbar( cnf, ax=ax0 )

    # cnf = ax1.plot_surface( Xe, Ye, Fr, cmap="viridis" )
    # plt.colorbar( cnf, ax=ax1 )

    # cnf = ax2.plot_surface( Xe, Ye, np.abs( Fr - Ffe ), cmap="viridis" )
    # plt.colorbar( cnf, ax=ax2 )

    # for i in range( Xe.shape[1] ):
    #     ax3.plot( xe, np.abs( Fr[:, i] ) )

    # for i in range( Ye.shape[0] ):
    #     ax4.plot( ye, np.abs( Fr[i, :] ) )

    # plt.show( )
    # raise ValueError( "Stop by user" )

    # print( "C_filter" )
    # print( C_filter )
    # print( "NCX_filter" )
    # print( NCX_filter )
    # print( "NCY_filter" )
    # print( NCY_filter )

    # Plot residual function
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
    cheby_rel_tol       = fit_props.cheby_rel_tol

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
    X2  = Xfp.ravel( )
    Y2  = Yfp.ravel( )
    Z2  = Zfp.ravel( )
    F2  = Ff.ravel( )
    # C   = fit_chebyshev_3d_robust( X2, Y2, Z2, F2, cheby_order_x-1, cheby_order_y-1, cheby_order_z-1)
    # C   = fit_chebyshev_3d( X2, Y2, Z2, F2, cheby_order_x, cheby_order_y, cheby_order_z)
    C = fit_props.sysmat_VT.T @ fit_props.sysmat_S_inv @ fit_props.sysmat_U.T @ F2
    
    # Filter coefficients to the precision required
    n_cheby_x       = linspace( 0, cheby_order_x-1, cheby_order_x, dtype=int )
    n_cheby_y       = linspace( 0, cheby_order_y-1, cheby_order_y, dtype=int )
    n_cheby_z       = linspace( 0, cheby_order_z-1, cheby_order_z, dtype=int )
    NCX,NCY,NCZ     = meshgrid( n_cheby_x, n_cheby_y, n_cheby_z, indexing="ij" )
    pos             = ( np_abs(C) > np_max( np_abs(C)) * cheby_rel_tol ) & ( np_abs(C) > cheby_abs_tol )
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
                ref_level_i                 = RefLevel( fit_props, parent=ref_level )

                # Check for refinement level
                if ref_level.level == fit_props.max_ref_level:
                    print( "--> MAXIMUM REFINEMENT LEVEL REACHED" )
                    print( f" - X: {fit_props.x_min:0.3E}, {fit_props.x_max:0.3E}" )
                    print( f" - Y: {fit_props.y_min:0.3E}, {fit_props.y_max:0.3E}" )

                # Fit function over the interval
                ref_level_i.add_data( *fit_integral_2d( fit_function, fit_props, show_figs=show_figs ) )

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


if __name__ == "__main__":
    # from fit_wave_infdepth_v2 import residual_region_11

    # C_filter = np.array( [
    #                         2.1632010219008144, # c[50]
    #                         0.6850473111366150, # c[51]
    #                         0.0101001298748349, # c[52]
    #                         -0.0018367082614362,#  c[53]
    #                         0.0001855458802984, # c[54]
    #                         -0.0000011136074591,#  c[55]
    #                         -0.0000028109106051,#  c[56]
    #                         -0.6851414361552296,#  c[57]
    #                         -0.0371467747792773,#  c[58]
    #                         0.0104120650434155, # c[59]
    #                         -0.0016397349226908,#  c[60]
    #                         0.0000920407770437, # c[61]
    #                         0.0000212334862666, # c[62]
    #                         -0.0000057128667158,#  c[63]
    #                         0.0101355646066731, # c[64]
    #                         -0.0104038972233652,#  c[65]
    #                         0.0024990294273947, # c[66]
    #                         -0.0002245780952191,#  c[67]
    #                         -0.0000431096175548,#  c[68]
    #                         0.0000174154078048, # c[69]
    #                         -0.0000021441787856,#  c[70]
    #                         0.0018307053603520, # c[71]
    #                         -0.0016408146885545,#  c[72]
    #                         0.0002245856053710, # c[73]
    #                         0.0000539995125187, # c[74]
    #                         -0.0000290438133964,#  c[75]
    #                         0.0000048409862906, # c[76]
    #                         0.0001861857882361, # c[77]
    #                         -0.0000919652227995,#  c[78]
    #                         -0.0000431086310019,#  c[79]
    #                         0.0000290439008627, # c[80]
    #                         -0.0000062405939360,#  c[81]
    #                         0.0000010598454522, # c[82]
    #                         0.0000212258828154, # c[83]
    #                         -0.0000174155261781,#  c[84]
    #                         0.0000048410534390, # c[85]
    #                         -0.0000028069843610,#  c[86]
    #                         0.0000057132649092, # c[87]
    #                         -0.0000021442416209,#  c[88]
    #                     ] )
    # NCX_filter  = np.array([
    #                             0,  # ncx[50]
    #                             0,  # ncx[51]
    #                             0,  # ncx[52]
    #                             0,  # ncx[53]
    #                             0,  # ncx[54]
    #                             0,  # ncx[55]
    #                             0,  # ncx[56]
    #                             1,  # ncx[57]
    #                             1,  # ncx[58]
    #                             1,  # ncx[59]
    #                             1,  # ncx[60]
    #                             1,  # ncx[61]
    #                             1,  # ncx[62]
    #                             1,  # ncx[63]
    #                             2,  # ncx[64]
    #                             2,  # ncx[65]
    #                             2,  # ncx[66]
    #                             2,  # ncx[67]
    #                             2,  # ncx[68]
    #                             2,  # ncx[69]
    #                             2,  # ncx[70]
    #                             3,  # ncx[71]
    #                             3,  # ncx[72]
    #                             3,  # ncx[73]
    #                             3,  # ncx[74]
    #                             3,  # ncx[75]
    #                             3,  # ncx[76]
    #                             4,  # ncx[77]
    #                             4,  # ncx[78]
    #                             4,  # ncx[79]
    #                             4,  # ncx[80]
    #                             4,  # ncx[81]
    #                             5,  # ncx[82]
    #                             5,  # ncx[83]
    #                             5,  # ncx[84]
    #                             5,  # ncx[85]
    #                             6,  # ncx[86]
    #                             6,  # ncx[87]
    #                             6,  # ncx[88]
    #                         ])
    
    # NCY_filter  = np.array( [
    #                             0,  # ncy[50]
    #                             1,  # ncy[51]
    #                             2,  # ncy[52]
    #                             3,  # ncy[53]
    #                             4,  # ncy[54]
    #                             5,  # ncy[55]
    #                             6,  # ncy[56]
    #                             0,  # ncy[57]
    #                             1,  # ncy[58]
    #                             2,  # ncy[59]
    #                             3,  # ncy[60]
    #                             4,  # ncy[61]
    #                             5,  # ncy[62]
    #                             6,  # ncy[63]
    #                             0,  # ncy[64]
    #                             1,  # ncy[65]
    #                             2,  # ncy[66]
    #                             3,  # ncy[67]
    #                             4,  # ncy[68]
    #                             5,  # ncy[69]
    #                             6,  # ncy[70]
    #                             0,  # ncy[71]
    #                             1,  # ncy[72]
    #                             2,  # ncy[73]
    #                             3,  # ncy[74]
    #                             4,  # ncy[75]
    #                             5,  # ncy[76]
    #                             0,  # ncy[77]
    #                             1,  # ncy[78]
    #                             2,  # ncy[79]
    #                             3,  # ncy[80]
    #                             4,  # ncy[81]
    #                             0,  # ncy[82]
    #                             1,  # ncy[83]
    #                             2,  # ncy[84]
    #                             3,  # ncy[85]
    #                             0,  # ncy[86]
    #                             1,  # ncy[87]
    #                             2,  # ncy[88]
    #                         ] )
    
    # fit_props               = FitProperties( )
    # fit_props.region_name   = "region_11"
    # fit_props.cheby_order_x = 20
    # fit_props.cheby_order_y = 20
    # fit_props.fcn_log_scale = False
    # fit_props.x_log_scale   = True
    # fit_props.x_max         = -4.375
    # fit_props.x_min         = -5.0
    # fit_props.y_log_scale   = True
    # fit_props.y_max         = -3.75
    # fit_props.y_min         = -4.375
    # fit_props.cheby_abs_tol = 1E-6
    # fit_props.cheby_rel_tol = 1E-14
    # fit_props.x_map_fcn     = lambda x: x
    # fit_props.y_map_fcn     = lambda y: y
    # fit_props.max_ref_level = 10

    # A, fit_stats = fit_integral_2d( residual_region_11, fit_props, show_figs=True, stop_to_show=True )
    # Xe  = -1.0
    # Ye  = -0.993417
    # Ye  = 0.9
    # X   = fit_props.x_map_lin( Xe )
    # Y   = fit_props.y_map_lin( Ye )
    # F   = eval_chebyshev_2d_filter(Xe, Ye, C_filter, NCX_filter, NCY_filter)
    # # F   = eval_chebyshev_2d_filter(Xe, Ye, A[:, 0], A[:, 1], A[:, 2])
    # F2  = residual_region_11( X, Y )
    # print( "X: ", X, " - Y: ", Y, " - Xe: ", Xe, " - Ye: ", Ye )
    # print( "Value: ", F )
    # print( "Value 2: ", F2 )
    # print( "Diff: ", F-F2 )

    class A:

        def __init__(self):
            self.x_min = 0.0
            self.x_max = 0.0
            self.V      = np.ndarray

    a = A( )
    a.x_min = 2.0
    a.x_max = 3.0
    a.V     = np.linspace( 0, 2, 3 )
    b = copy.copy( a )

    b.x_min = 4.0
    b.x_max = 5.0
    b.V[0]  = 999

    print( a.x_min, a.x_max, a.V )
    print( b.x_min, b.x_max, b.V )

