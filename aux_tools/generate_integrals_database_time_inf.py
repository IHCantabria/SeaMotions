
# Import general usage libraries
import copy
import h5py

# Import general usage scientific libraries
import numpy as np
import scipy as sp
from scipy.special import jv

# Import general usage scientific plotting libraries
import matplotlib.pyplot as plt
from matplotlib import cm
import mplcursors

# Import local modules
from fit_cheby_v2 import FitProperties, fit_integral_1d, fit_integral_2d, RefLevel, fit_residual_1D_adaptive_interface, fit_residual_2D_adaptive_interface
from fit_tools import RefLevel, write_coeffs_module_adaptive_1d_only_header, write_coeffs_module_adaptive_2d_only_header


class ChebyFit:

    def __init__( self, rhs: np.ndarray, show_fit=False ) -> None:
        # Define class properties
        self._coeffs    = np.ndarray
        self._vandermat = np.ndarray
        self._rhs       = rhs
        self._show_fit  = show_fit

        # Initialize class by calculating the interpolation coefficients
        self._initialize( )

    def __call__( self, x: np.ndarray ) -> np.ndarray:
        return self.evaluate( x )

    def _initialize( self ) -> None:
        t               = sp.special.roots_chebyt( self._rhs.shape[0] )[0]
        self._vandermat = np.polynomial.chebyshev.chebvander( t, t.shape[0]-1 )
        self._coeffs    = np.linalg.solve( self._vandermat, self._rhs )

        if self._show_fit:
            t_eval      = np.linspace( -1, 1, 100 )
            recon_fit   = np.polynomial.chebyshev.chebval( t, self._coeffs )
            recon       = np.polynomial.chebyshev.chebval( t_eval, self._coeffs )

            fig, ax = plt.subplots( 2, 1, figsize= ( 8, 8 ) )
            ax[0].plot( t, self._rhs, label='original' )
            ax[0].plot( t, recon_fit, '--', label=f'recon fit' )
            ax[0].plot( t_eval, recon, '--', label=f'recon random' )
            ax[0].plot( t, self._rhs - recon_fit, '--', label=f'Difference (recon fit)' )
            ax[0].legend()
            ax[0].set_xlabel( "t" )
            ax[0].set_ylabel( 'Fit' )
            ax[0].grid(True)

            ax[1].semilogy(np.abs(self._coeffs), 'o-')
            ax[1].set_xlabel('order n')
            ax[1].set_ylabel('|a_n|')
            ax[1].set_title('Coefficients by order')
            ax[1].grid(True)

            plt.show()

    def evaluate( self, x: np.ndarray ) -> np.ndarray:
        return np.polynomial.chebyshev.chebval( x, self._coeffs )
    

class ChebyFitHP:

    def __init__( self, sp_fit ) -> None:
        # Define class attribues
        self._fit_fcn   = RefLevel
        self._spf       = sp_fit

        # Start fit
        self._initialize( )

    def __call__( self, xe: float ) -> float:
        return self.evaluate( xe )

    def _fit( self ) -> RefLevel:
        # Define fit properties
        fit_props               = FitProperties( )
        fit_props.dims          = 1
        fit_props.region_name   = "ChebyFitHP"
        fit_props.cheby_order_x = 5
        fit_props.fcn_log_scale = False
        fit_props.x_log_scale   = False
        fit_props.x_max         = self._spf._x[-1]
        fit_props.x_min         = self._spf._x[0]
        fit_props.cheby_abs_tol = 1E-4
        fit_props.cheby_rel_tol = 1E-14
        fit_props.x_map_fcn     = lambda x: x
        fit_props.y_map_fcn     = lambda y: y
        fit_props.max_ref_level = 2

        fit_props.num_x         = fit_props.cheby_order_x + 1
        fit_props.num_x_fit     = fit_props.cheby_order_x + 1
        fit_props.num_y         = fit_props.cheby_order_y + 1
        fit_props.num_y_fit     = fit_props.cheby_order_y + 1
        
        fit_props.generate_fitting_matrix( )

        fit_function            = self._spf

        # Create root FitRegion
        ref_level               = RefLevel( copy.copy( fit_props ) )

        # Create first fit to feed root refinement level
        ref_level.add_data( *fit_integral_1d( fit_function, fit_props, show_figs=False ) )

        # Start adaptive refinement loop
        fit_residual_1D_adaptive_interface( fit_function, ref_level )

        # Set starting position
        ref_level.set_start_index( 0 )

        return ref_level

    def _initialize( self ) -> None:
        # Fit to target function by using adaptive interface
        self._fit_fcn = self._fit( )

        # Recompose the adaptive fit to have a working interpolation scheme
        self._recompose_fit( self._fit_fcn )

    def _recompose_fit( self, ref_level: RefLevel ) -> None:
        # Get cumulative data
        self.cum_coeffs      = ref_level.get_cheby_coeffs( )
        self.cheby_coeffs    = self.cum_coeffs[:, 0]
        self.ncx             = self.cum_coeffs[:, 1].astype( int )
        self.max_cheby_order = int( self.cum_coeffs[:, 1:].max( ) )

        # Get intervals size
        self.max_ref_level   = ref_level.get_max_level( )
        self.x_min           = ref_level.fit_props.x_min
        self.x_max           = ref_level.fit_props.x_max

        # Get log scales
        self.x_log_scale     = ref_level.fit_props.x_log_scale
        self.fcn_log_scale   = ref_level.fit_props.fcn_log_scale

        # Calculate hash table
        self.intervals_np    = 2**(self.max_ref_level)
        self.dx              = ( self.x_max - self.x_min ) / self.intervals_np
        self.x               = np.arange( self.x_min, self.x_max+self.dx, self.dx )
        self.xm              = ( self.x[1:] + self.x[:-1] ) / 2.0

        self.blocks_np                  = self.xm.shape[0]
        self.blocks_start               = np.zeros( ( self.blocks_np, ), dtype=int )
        self.blocks_coeffs_np           = np.zeros( ( self.blocks_np, ), dtype=int )
        self.blocks_start_f             = np.zeros( ( self.blocks_np, ), dtype=int )
        self.blocks_coeffs_np_f         = np.zeros( ( self.blocks_np, ), dtype=int )
        self.blocks_max_cheby_order     = np.zeros( ( self.blocks_np, ), dtype=int )
        self.blocks_max_cheby_order_f   = np.zeros( ( self.blocks_np, ), dtype=int )
        self.x_min_vec                  = np.zeros( ( self.blocks_np, ), dtype=float )
        self.x_max_vec                  = np.zeros( ( self.blocks_np, ), dtype=float )
        count                           = 0
        for i in range( self.xm.shape[0] ):
            xmi = self.xm[i]
            print( "I: ", i, xmi )
            bs, bc, bsf, bcf, mco, mcof, x_min_i, x_max_i, _, _, _, _   = ref_level.get_start_index( np.array( [ xmi ] ) )
            self.blocks_start[count]                                    = bs
            self.blocks_coeffs_np[count]                                = bc
            self.blocks_start_f[count]                                  = bsf
            self.blocks_coeffs_np_f[count]                              = bcf
            self.blocks_max_cheby_order[count]                          = mco
            self.blocks_max_cheby_order_f[count]                        = mcof
            self.x_min_vec[count]                                       = x_min_i
            self.x_max_vec[count]                                       = x_max_i
            count                                                       += 1

    def evaluate( self, xe: float ) -> float:
        # Find interval
        npos = int( ( xe - self.x_min ) / self.dx )
        if npos > ( self.x_min_vec.shape[0] - 1 ): 
            npos = ( self.x_min_vec.shape[0] - 1 )
         
        # Convert Xe to local coordinates [-1, 1]
        xeu = map_to_unit( self.x_min_vec[ npos ], self.x_max_vec[ npos ], xe )

        # Calculate interpolated value
        coeffs_i    = np.zeros( ( int( self.blocks_max_cheby_order[ npos ]+1 ), ) )
        for i in range( self.blocks_coeffs_np[ npos ] ):
            index                           = self.blocks_start[ npos ] + i
            coeffs_i[ self.ncx[ index ] ]   = self.cheby_coeffs[ index ]

        
        val = np.polynomial.chebyshev.chebval( xeu, coeffs_i ).sum( )

        return val
    
    def get_fit_fcn( self ) -> RefLevel:
        return self._fit_fcn
    

class ChebyFitHP2D:

    def __init__( self, data_fipath: str, fit_props: FitProperties ) -> None:
        # Define class attribues
        self._fit_props     = fit_props
        self._data_fipath   = data_fipath

        # Start fit
        self._initialize( )

    def __call__( self, xe: float, ye: float ) -> float:
        return self.evaluate( xe, ye )

    def _fit( self ) -> RefLevel:
        # Load database
        with h5py.File( self._data_fipath, "r" ) as fid:
            mu          = fid[ "mu" ][:]
            beta        = fid[ "beta" ][:]
            data_raw    = fid[ "fcn" ][:]
        
        # Interpolate database
        fit_function_raw        = sp.interpolate.RegularGridInterpolator( ( beta, mu ), data_raw )
        fit_function            = lambda x, y: fit_function_raw( ( x, y ) )

        # Create root FitRegion
        ref_level               = RefLevel( copy.copy( self._fit_props ) )

        # Configure fit additional settings
        is_square_ref           = True
        show_figs               = False

        # Create first fit to feed root refinement level
        ref_level.add_data( *fit_integral_2d( fit_function, self._fit_props, show_figs=show_figs ) )

        y_levels = np.array( [ -4.0, -3.0, -2.0, -1.0, np.log10( 0.9998 ) ] )
        for i in range( 39 ):
            for j in range( y_levels.shape[0]-1 ):
                fit_props_i         = copy.copy( ref_level.fit_props )
                fit_props_i.x_max   = ( i + 1 )
                fit_props_i.x_min   = i
                fit_props_i.y_max   = y_levels[j+1]
                fit_props_i.y_min   = y_levels[j]
                
                ref_level_i         = RefLevel( fit_props_i, parent=ref_level )

                # Create first fit to feed root refinement level
                ref_level_i.add_data( *fit_integral_2d( fit_function, fit_props_i, show_figs=show_figs ) )

                # Start adaptive refinement loop
                if not ref_level_i.check_tolerances( ) and ref_level_i.level < ref_level_i.fit_props.max_ref_level:
                    fit_residual_2D_adaptive_interface( fit_function, ref_level_i, is_square_ref=is_square_ref, show_figs=show_figs )

                # Storage child data
                ref_level.child.append( ref_level_i )

        # Set starting position
        ref_level.set_start_index( 0 )

        # Plot results summary
        ref_level.show_summary( folder_path )

        print( ref_level )

        return ref_level
    
    def _get_coeffs( self, xe: float, ye: float ) -> list:
        # Find interval
        xnpos = int( ( xe - self.x_min ) / self.dx )
        if xnpos > ( self.x_min_vec.shape[0] - 1 ): 
            xnpos = ( self.x_min_vec.shape[0] - 1 )

        ynpos = int( ( ye - self.y_min ) / self.dy )
        if ynpos > ( self.y_min_vec.shape[0] - 1 ): 
            ynpos = ( self.y_min_vec.shape[0] - 1 )

        npos = xnpos * self.intervals_np + ynpos
         
        # Calculate interpolated value
        n           = int( self.blocks_max_cheby_order[ npos ] + 1 )
        coeffs_i    = np.zeros( ( n, n ) )
        for i in range( self.blocks_coeffs_np[ npos ] ):
            index                                               = self.blocks_start[ npos ] + i
            coeffs_i[ self.ncx[ index ], self.ncy[ index ] ]    = self.cheby_coeffs[ index ]

        return coeffs_i

    def _initialize( self ) -> None:
        # Fit to target function by using adaptive interface
        ref_level = self._fit( )

        # Recompose the adaptive fit to have a working interpolation scheme
        self._recompose_fit( ref_level )

    def _recompose_fit( self, ref_level: RefLevel ) -> None:
        # Get cumulative data
        self.cum_coeffs      = ref_level.get_cheby_coeffs( )
        self.cheby_coeffs    = self.cum_coeffs[:, 0]
        self.ncx             = self.cum_coeffs[:, 1]
        self.ncy             = self.cum_coeffs[:, 2]
        self.max_cheby_order = self.cum_coeffs[:, 1:].max( )

        # Get intervals size
        self.max_ref_level   = ref_level.get_max_level( )
        self.x_min           = ref_level.fit_props.x_min
        self.x_max           = ref_level.fit_props.x_max
        self.y_min           = ref_level.fit_props.y_min
        self.y_max           = ref_level.fit_props.y_max

        # Get log scales
        self.x_log_scale     = ref_level.fit_props.x_log_scale
        self.y_log_scale     = ref_level.fit_props.y_log_scale
        self.fcn_log_scale   = ref_level.fit_props.fcn_log_scale

        # Calculate hash table
        self.intervals_np   = 2**(self.max_ref_level)
        self.dx             = ( self.x_max - self.x_min ) / self.intervals_np
        self.x              = np.linspace( self.x_min, self.x_max, self.intervals_np+1 )
        self.xm             = ( self.x[1:] + self.x[:-1] ) / 2.0

        self.dy             = ( self.y_max - self.y_min ) / self.intervals_np
        self.y              = np.linspace( self.y_min, self.y_max, self.intervals_np+1 )
        self.ym             = ( self.y[1:] + self.y[:-1] ) / 2.0

        self.blocks_np              = self.xm.shape[0] * self.ym.shape[0]
        self.blocks_start           = np.zeros( ( self.blocks_np, ), dtype=int )
        self.blocks_coeffs_np       = np.zeros( ( self.blocks_np, ), dtype=int )
        self.blocks_max_cheby_order = np.zeros( ( self.blocks_np, ), dtype=int )
        self.x_min_vec              = np.zeros( ( self.blocks_np, ), dtype=float )
        self.x_max_vec              = np.zeros( ( self.blocks_np, ), dtype=float )
        self.y_min_vec              = np.zeros( ( self.blocks_np, ), dtype=float )
        self.y_max_vec              = np.zeros( ( self.blocks_np, ), dtype=float )
        count                       = 0
        for i in range( self.xm.shape[0] ):
            xmi = self.xm[i]
            for j in range( self.ym.shape[0] ):
                ymi = self.ym[j]
                print( i, j, self.xm[i], self.ym[i] )
                bs, bc, _, _, mco, _, x_min_i, x_max_i, y_min_i, y_max_i, _, _      = ref_level.get_start_index( np.array( [ xmi, ymi ] ) )
                self.blocks_start[count]                                            = bs
                self.blocks_coeffs_np[count]                                        = bc
                self.blocks_max_cheby_order[count]                                  = mco
                self.x_min_vec[count]                                               = x_min_i
                self.x_max_vec[count]                                               = x_max_i
                self.y_min_vec[count]                                               = y_min_i
                self.y_max_vec[count]                                               = y_max_i
                count                                                               += 1

    def evaluate( self, xe: float, ye: float ) -> float:
        # Get patch properties
        coeffs_i, npos = self._get_coeffs( xe, ye )

        # Convert Xe to local coordinates [-1, 1]
        xeu = map_to_unit( self.x_min_vec[ npos ], self.x_max_vec[ npos ], xe )
        yeu = map_to_unit( self.y_min_vec[ npos ], self.y_max_vec[ npos ], ye )

        # Calculate polynomial value
        val = np.polynomial.chebyshev.chebval2d( xeu, yeu, coeffs_i ).sum( )

        return val
    
    def evaluate_dx( self, xe: float, ye: float ) -> float:
        # Get patch properties
        coeffs_i, npos = self._get_coeffs( xe, ye )
        
        # Convert Xe to local coordinates [-1, 1]
        xeu = map_to_unit( self.x_min_vec[ npos ], self.x_max_vec[ npos ], xe )
        yeu = map_to_unit( self.y_min_vec[ npos ], self.y_max_vec[ npos ], ye )

        # Derivative w.r.t x
        coeffs_i = np.polynomial.chebyshev.chebder( coeffs_i, m=1, axis=0 )

        # Calculate polynomial value
        val = np.polynomial.chebyshev.chebval2d( xeu, yeu, coeffs_i ).sum( )

        return val
    
    def evaluate_dy( self, xe: float, ye: float ) -> float:
        # Get patch properties
        coeffs_i, npos = self._get_coeffs( xe, ye )
        
        # Convert Xe to local coordinates [-1, 1]
        xeu = map_to_unit( self.x_min_vec[ npos ], self.x_max_vec[ npos ], xe )
        yeu = map_to_unit( self.y_min_vec[ npos ], self.y_max_vec[ npos ], ye )

        # Derivative w.r.t x
        coeffs_i = np.polynomial.chebyshev.chebder( coeffs_i, m=1, axis=1 )

        # Calculate polynomial value
        val = np.polynomial.chebyshev.chebval2d( xeu, yeu, coeffs_i ).sum( )

        return val
    

class ChebyFitHPAlpha:

    def __init__( self, magee_fit ) -> None:
        # Define class attributes
        self._mageefit2     = magee_fit
        self._interp_fcn    = [ ]

        # Initialize class
        self._initialize( )

    def __call__( self, x: float ) -> np.ndarray:
        return self.evaluate( x )
    
    def _initialize( self ) -> None:
        # Loop over Magee interpolators to regenerate with Chebyshev polynomials
        for mi in self._mageefit2.get_interpolators( ):
            self._interp_fcn.append( ChebyFitHP( mi ) )
    
    def evaluate( self, x: float ) -> np.ndarray:
        result = np.zeros( ( len( self._interp_fcn, ) ) )
        for i in range( result.shape[0] ):
            result[i] = self._interp_fcn[i]( x )

        return result
    
    def get_interpolators( self ) -> list:
        return self._interp_fcn


class SplineFit:

    def __init__( self, x: np.ndarray, y:np.ndarray, show_fit=False ) -> None:
        # Define class attributes
        self._interpf   = None
        self._x         = x
        self._y         = y
        self._show_fit  = show_fit

        # Initialize class
        self._initialize( )

    def __call__( self, x: np.ndarray ) -> np.ndarray:
        return self.evaluate( x )

    def _initialize( self ) -> None:
        # Generate interpolation object
        self._interpf = sp.interpolate.CubicSpline( self._x, self._y )

        # Show interpolation quality
        if self._show_fit:
            xx = np.linspace( self._x[0], self._x[-1], 100 )

            fig = plt.figure( )
            ax  = fig.add_subplot( 111 )

            ax.plot( xx, self._interpf( xx ), "--", label="Interpolated" )
            ax.plot( self._x, self._y, "o", color="orange", label="Original" )
            ax.legend( )

            plt.show( )

    def evaluate( self, x: np.ndarray ) -> np.ndarray:
        return self._interpf( x )
    

class MageeFit:

    def __init__( self, coeffs: np.ndarray ) -> None:
        # Define class attributes
        self._coeffs        = coeffs
        self._interp_fcn    = [ ]

        # Initialize class by fitting each of the coefficients using chebyshev expansion
        self._initialize( )

    def __call__( self, x: float ) -> np.ndarray:
        return self.evaluate( x )

    def _initialize( self ) -> None:
        for i in range( self._coeffs.shape[0] ):
            self._interp_fcn.append( ChebyFit( self._coeffs[i, :], show_fit=True ) )

    def evaluate( self, x: float ) -> np.ndarray:
        results = np.zeros( ( self._coeffs.shape[0], ) )
        for i in range( self._coeffs.shape[0] ):
            results[i] = self._interp_fcn[i]( x )

        return results
    

class MageeFit2:

    def __init__( self, mu: np.ndarray, coeffs: np.ndarray ) -> None:
        # Define class attributes
        self._coeffs        = coeffs
        self._interp_fcn    = [ ]
        self._mu            = mu

        # Initialize class by fitting each of the coefficients using chebyshev expansion
        self._initialize( )

    def __call__( self, x: float ) -> np.ndarray:
        return self.evaluate( x )

    def _initialize( self ) -> None:
        for i in range( self._coeffs.shape[0] ):
            self._interp_fcn.append( SplineFit( self._mu, self._coeffs[i, :], show_fit=True ) )

    def evaluate( self, x: float ) -> np.ndarray:
        results = np.zeros( ( self._coeffs.shape[0], ) )
        for i in range( self._coeffs.shape[0] ):
            results[i] = self._interp_fcn[i]( x )

        return results
    
    def get_interpolators( self ) -> list:
        return self._interp_fcn


def add_hover_tool( X, Y, Z, cnf ) -> None:
    # Create hover tool
    cursor = mplcursors.cursor( cnf, hover=True )

    # Add annotation
    @cursor.connect( "add" )
    def on_add( sel ) -> None:
        i, j = sel.index  # indices of the hovered cell
        print( i, j )
        x_val, y_val = X[i, j], Y[i, j]
        z_val = Z[i, j]
        sel.annotation.set_text(f"x={x_val:.2f}\ny={y_val:.2f}\nz={z_val:.2f}")


def calculate_residual_magee( fit_props: FitProperties, target_fcn, interp_raw, data_raw, beta: np.ndarray, mu: np.ndarray ) -> None:
    # Loop over magee interpolation points to fit coefficients
    for i in range( fit_props.magee_cheby_np ):
        mui = np.ones_like( beta ) * fit_props.magee_mui[i]
        print( mui )
        _, fit_props.magee_coeffs[:, i] = fit_magee_fcn( target_fcn, beta, fit_props.magee_mui[i], interp_raw( ( beta, mui ) ), fit_props.alpha_shift, show_recon=False )

    # magee_fit = MageeFit( fit_props.magee_coeffs )
    if fit_props.magee_cheby_np > 1:
        mu_fit              = fit_props.magee_mui_log
        magee_coeffs_fit    = fit_props.magee_coeffs

    else:
        np_fit              = 4
        mu_fit              = np.linspace( fit_props.y_min, fit_props.y_max, np_fit )
        magee_coeffs_fit    = fit_props.magee_coeffs * np.ones( ( np_fit, ) )

    magee_fit     = MageeFit2( mu_fit, magee_coeffs_fit )
    magee_chfhpa  = ChebyFitHPAlpha( magee_fit )

    data = data_raw.copy( )
    G0   = np.zeros( ( data.shape[0], ) )
    for i in range( mu.shape[0] ):
        print( "I: ", i )
        if mu[i] < 10**fit_props.magee_cheby_intv[1]:
            print( "Re-fitting..." )
            # coeffs_i    = magee_fit( map_to_unit( fit_props.magee_cheby_intv[0], fit_props.magee_cheby_intv[1], np.log10( mu[i] ) ) )
            # coeffs_i    = magee_fit( np.log10( mu[i] ) )
            coeffs_i    = magee_chfhpa( np.log10( mu[i] ) )
        
        G0          = eval_magee_fcn( target_fcn, coeffs_i, beta, mu[i], fit_props.alpha_shift )
        data[:, i]  = data_raw[:, i] - G0
        # data[:, i]  = data_raw[:, i] - target_fcn( beta, mu[i] )


    fit_props.magee_norm_f  = np.max( np.abs( data ) ) * 10
    # data_norm               = data / fit_props.magee_norm_f
    # B, M                    = np.meshgrid( beta, mu, indexing="ij" )
    # data                    = data_raw * np.exp( - 1/( mu )**0.3 )
    data_norm               = data.copy( )
    print( "norm_f: ", fit_props.magee_norm_f )

    return data, data_norm, magee_chfhpa


def magee_dalpha( tb: np.ndarray, alpha: float ) -> np.ndarray:
    g       = (
                    jv( 0.25, tb**2.0 / 8.0 + alpha )
                    *
                    jv( -0.25, tb**2.0 / 8.0 + alpha )
                    +
                    jv( 0.75, tb**2.0 / 8.0 + alpha )
                    *
                    jv( -0.75, tb**2.0 / 8.0 + alpha )
                )
    pos     = np.isnan( g )
    g[pos]  = 0.0
    pos     = tb < 1e-1
    g[pos]  = 0.0

    return g


def eval_magee_fcn( fcn_fit, coeffs: np.ndarray, beta: np.ndarray, mu: float, alpha: np.ndarray, show_recon=False ) -> np.ndarray:
    # Build system matrix
    A = np.zeros( ( beta.shape[0], alpha.shape[0] ) )
    for i in range( alpha.shape[0] ):
        A[:, i] = fcn_fit( beta, mu, alpha=alpha[i] )

    # Compute function value
    recon = A.dot( coeffs )

    if show_recon:
        fig, ax = plt.subplots( 2, 2 )

        cnf = ax[0, 0].matshow( np.log10( np.abs( A ) ) )
        plt.colorbar( cnf, ax=ax[0, 0] )
        cnf = ax[0, 1].matshow( np.log10( np.abs( coeffs ) ) )
        plt.colorbar( cnf, ax=ax[0, 1] )
        cnf = ax[1, 0].matshow( np.log10( np.abs( recon ) ) )
        plt.colorbar( cnf, ax=ax[1, 0] )

        plt.show( )

    return recon


def eval_magee_fcn_mu( fcn_fit, coeffs: np.ndarray, beta: float, mu: np.ndarray, alpha: np.ndarray, show_recon=False ) -> np.ndarray:
    # Build system matrix
    A = np.zeros( ( mu.shape[0], alpha.shape[0] ) )
    for i in range( alpha.shape[0] ):
        A[:, i] = fcn_fit( beta, mu, alpha=alpha[i] )

    # Compute function value
    recon = A.dot( coeffs )

    if show_recon:
        fig, ax = plt.subplots( 2, 2 )

        cnf = ax[0, 0].matshow( np.log10( np.abs( A ) ) )
        plt.colorbar( cnf, ax=ax[0, 0] )
        cnf = ax[0, 1].matshow( np.log10( np.abs( coeffs ) ) )
        plt.colorbar( cnf, ax=ax[0, 1] )
        cnf = ax[1, 0].matshow( np.log10( np.abs( recon ) ) )
        plt.colorbar( cnf, ax=ax[1, 0] )

        plt.show( )

    return recon


def eval_magee_fcn_v2( coeffs: np.ndarray, beta: np.ndarray, mu: float ) -> np.ndarray:
    # Build system matrix
    A           = np.zeros( ( beta.shape[0], 3 ) )
    lt          = - np.pi * beta**5.0 / 64 / np.sqrt( 2 ) * np.exp( -beta**2.0 * mu / 4.0 )
    nf          = np.max( np.abs( lt ) )
    
    if np.abs( nf ) < 1e-6:
        A[:, 0] = 0.0
        A[:, 1] = 0.0
        A[:, 2] = 0.0
    
    else:
        A[:, 0] = lt * magee_dalpha( beta, -1e-4 ) / nf
        A[:, 1] = lt * magee_dalpha( beta, +0.0  ) / nf
        A[:, 2] = lt * magee_dalpha( beta, +1e-4 ) / nf

    # Compute function value
    recon = ( A.dot( coeffs ) ) * nf

    # fig, ax = plt.subplots( 2, 2 )

    # cnf = ax[0, 0].matshow( np.log10( np.abs( A ) ) )
    # plt.colorbar( cnf, ax=ax[0, 0] )
    # cnf = ax[0, 1].matshow( np.log10( np.abs( coeffs ) ) )
    # plt.colorbar( cnf, ax=ax[0, 1] )
    # cnf = ax[1, 0].matshow( np.log10( np.abs( recon ) ) )
    # plt.colorbar( cnf, ax=ax[1, 0] )

    # plt.show( )

    return recon


def eval_magee_fcn_v3( coeffs: np.ndarray, beta: np.ndarray, mu: float ) -> np.ndarray:
    # Build system matrix
    A           = np.zeros( ( beta.shape[0], 3 ) )
    lt          = - np.pi * beta**7.0 / 256 / np.sqrt( 2 ) * np.exp( -beta**2.0 * mu / 4.0 )
    nf          = np.max( np.abs( lt ) )
    
    if np.abs( nf ) < 1e-6:
        A[:, 0] = 0.0
        A[:, 1] = 0.0
        A[:, 2] = 0.0
    
    else:
        A[:, 0] = lt * magee_dalpha( beta, -1e-4 ) / nf
        A[:, 1] = lt * magee_dalpha( beta, +0.0  ) / nf
        A[:, 2] = lt * magee_dalpha( beta, +1e-4 ) / nf

    # Compute function value
    recon = ( A.dot( coeffs ) ) * nf

    # fig, ax = plt.subplots( 2, 2 )

    # cnf = ax[0, 0].matshow( np.log10( np.abs( A ) ) )
    # plt.colorbar( cnf, ax=ax[0, 0] )
    # cnf = ax[0, 1].matshow( np.log10( np.abs( coeffs ) ) )
    # plt.colorbar( cnf, ax=ax[0, 1] )
    # cnf = ax[1, 0].matshow( np.log10( np.abs( recon ) ) )
    # plt.colorbar( cnf, ax=ax[1, 0] )

    # plt.show( )

    return recon


def fit_magee_fcn( fcn_fit, beta: np.ndarray, mu: float, f: np.ndarray, alpha: np.ndarray, show_recon=False ) -> None:
    # Build system matrix
    A = np.zeros( ( beta.shape[0], alpha.shape[0] ) )
    for i in range( alpha.shape[0] ):
        A[:, i] = fcn_fit( beta, mu, alpha=alpha[i] )

    # Weights: use radial weight r (if appropriate) ---
    # avoid exact zero weight at r=0 if using r: set small positive weight for first sample
    w       = np.linspace( 0, 1.0, beta.shape[0] )
    w[0]    = 1e-6
    W_sqrt  = np.sqrt(w)

    A_w     = A * W_sqrt[:, None]
    f_w     = f * W_sqrt

    # Solve weighted least-squares (with optional Tikhonov regularization) ---
    # plain least-squares:
    coeffs, residuals, rank, s = np.linalg.lstsq(A_w, f_w)

    # optional regularized solver (lambda parameter):
    # lam = 1e-6
    # ATA = A_w.T @ A_w + lam * np.eye(N+1)
    # ATb = A_w.T @ f_w
    # coeffs = np.linalg.solve(ATA, ATb)

    # --- reconstruct and diagnostics ---
    recon = A @ coeffs

    if show_recon:
        fig, ax = plt.subplots( 2, 1, figsize= ( 8, 8 ) )
        ax[0].plot( beta, f, label='original' )
        ax[0].plot( beta, recon, '--', label=f'recon' )
        ax[0].plot( beta, f - recon, '--', label=f'Difference' )
        ax[0].legend()
        ax[0].set_xlabel( "beta" )
        ax[0].set_ylabel( 'Fit' )
        ax[0].grid(True)

        ax[1].semilogy(np.abs(coeffs), 'o-')
        ax[1].set_xlabel('order n')
        ax[1].set_ylabel('|a_n|')
        ax[1].set_title('Coefficients by order')
        ax[1].grid(True)

        plt.show()

    return recon, coeffs


def fit_magee_fcn_v2( fcn_fit, beta: np.ndarray, mu: float, f: np.ndarray, N: int ) -> None:
    # print( "mu: ", mu )
    NN = N
    if mu > 0.3:
        NN = 1

    # --- build design matrix A_{j,n} = J_n(k r_j) ---
    A           = np.zeros((beta.shape[0], 3))
    lt          = - np.pi * beta**5.0 / 64 / np.sqrt( 2 ) * np.exp( -beta**2.0 * mu / 4.0 )
    nf          = np.max( np.abs( lt ) )
    A[:, 0]     = lt * magee_dalpha( beta, -1e-4 ) / nf
    A[:, 1]     = lt * magee_dalpha( beta, +0.0  ) / nf
    A[:, 2]     = lt * magee_dalpha( beta, +1e-4 ) / nf

    # --- weights: use radial weight r (if appropriate) ---
    # avoid exact zero weight at r=0 if using r: set small positive weight for first sample
    w       = np.linspace( 0, 1.0, beta.shape[0] )
    w[0]    = 1e-6
    W_sqrt  = np.sqrt(w)

    A_w     = A * W_sqrt[:, None]
    f_w     = f * W_sqrt / nf

    # --- solve weighted least-squares (with optional Tikhonov regularization) ---
    # plain least-squares:
    _coeffs, residuals, rank, s = np.linalg.lstsq(A_w, f_w)

    coeffs          = _coeffs.copy( )

    # print( "coeffs: ", coeffs )

    # optional regularized solver (lambda parameter):
    # lam = 1e-6
    # ATA = A_w.T @ A_w + lam * np.eye(N+1)
    # ATb = A_w.T @ f_w
    # coeffs = np.linalg.solve(ATA, ATb)

    # --- reconstruct and diagnostics ---
    recon = ( A @ _coeffs ) * nf

    # fig, ax = plt.subplots( 2, 1, figsize= ( 8, 8 ) )
    # ax[0].plot( beta, f, label='original' )
    # ax[0].plot( beta, recon, '--', label=f'recon N={N}' )
    # ax[0].plot( beta, f - recon, '--', label=f'Difference' )
    # ax[0].legend()
    # ax[0].set_xlabel( "beta" )
    # ax[0].set_ylabel( 'Fit' )
    # ax[0].grid(True)

    # ax[1].semilogy(np.abs(coeffs), 'o-')
    # ax[1].set_xlabel('order n')
    # ax[1].set_ylabel('|a_n|')
    # ax[1].set_title('Coefficients by order')
    # ax[1].grid(True)

    # plt.show()

    return recon, coeffs


def fit_magee_fcn_v3( fcn_fit, beta: np.ndarray, mu: float, f: np.ndarray, N: int ) -> None:
    # print( "mu: ", mu )
    NN = N
    if mu > 0.3:
        NN = 1

    # --- build design matrix A_{j,n} = J_n(k r_j) ---
    A           = np.zeros((beta.shape[0], 3))
    lt          = - np.pi * beta**7.0 / 256 / np.sqrt( 2 ) * np.exp( -beta**2.0 * mu / 4.0 )
    nf          = np.max( np.abs( lt ) )
    A[:, 0]     = lt * magee_dalpha( beta, -1e-4 ) / nf
    A[:, 1]     = lt * magee_dalpha( beta, +0.0  ) / nf
    A[:, 2]     = lt * magee_dalpha( beta, +1e-4 ) / nf

    # --- weights: use radial weight r (if appropriate) ---
    # avoid exact zero weight at r=0 if using r: set small positive weight for first sample
    w       = np.linspace( 0, 1.0, beta.shape[0] )
    w[0]    = 1e-6
    W_sqrt  = np.sqrt(w)

    A_w     = A * W_sqrt[:, None]
    f_w     = f * W_sqrt / nf

    # --- solve weighted least-squares (with optional Tikhonov regularization) ---
    # plain least-squares:
    _coeffs, residuals, rank, s = np.linalg.lstsq(A_w, f_w)

    coeffs          = _coeffs.copy( )

    # print( "coeffs: ", coeffs )

    # optional regularized solver (lambda parameter):
    # lam = 1e-6
    # ATA = A_w.T @ A_w + lam * np.eye(N+1)
    # ATb = A_w.T @ f_w
    # coeffs = np.linalg.solve(ATA, ATb)

    # --- reconstruct and diagnostics ---
    recon = ( A @ _coeffs ) * nf

    # fig, ax = plt.subplots( 2, 1, figsize= ( 8, 8 ) )
    # ax[0].plot( beta, f, label='original' )
    # ax[0].plot( beta, recon, '--', label=f'recon N={N}' )
    # ax[0].plot( beta, f - recon, '--', label=f'Difference' )
    # ax[0].legend()
    # ax[0].set_xlabel( "beta" )
    # ax[0].set_ylabel( 'Fit' )
    # ax[0].grid(True)

    # ax[1].semilogy(np.abs(coeffs), 'o-')
    # ax[1].set_xlabel('order n')
    # ax[1].set_ylabel('|a_n|')
    # ax[1].set_title('Coefficients by order')
    # ax[1].grid(True)

    # plt.show()

    return recon, coeffs


def fit_magee_fcn_v4( fcn_fit, beta: np.ndarray, mu: float, f: np.ndarray, N: int ) -> None:
    # print( "mu: ", mu )
    NN = N
    if mu > 0.3:
        NN = 1

    # --- build design matrix A_{j,n} = J_n(k r_j) ---
    A           = np.zeros((beta.shape[0], 3))
    lt          = - np.pi * beta**7.0 / 256 / np.sqrt( 2 ) * np.exp( -beta**2.0 * mu / 4.0 )
    nf          = np.max( np.abs( lt ) )
    A[:, 0]     = lt * magee_dalpha( beta, -1e-4 ) / nf
    A[:, 1]     = lt * magee_dalpha( beta, +0.0  ) / nf
    A[:, 2]     = lt * magee_dalpha( beta, +1e-4 ) / nf

    # --- weights: use radial weight r (if appropriate) ---
    # avoid exact zero weight at r=0 if using r: set small positive weight for first sample
    w       = np.linspace( 0, 1.0, beta.shape[0] )
    w[0]    = 1e-6
    W_sqrt  = np.sqrt(w)

    A_w     = A * W_sqrt[:, None]
    f_w     = f * W_sqrt / nf

    # --- solve weighted least-squares (with optional Tikhonov regularization) ---
    # plain least-squares:
    _coeffs, residuals, rank, s = np.linalg.lstsq(A_w, f_w)

    coeffs          = _coeffs.copy( )

    # print( "coeffs: ", coeffs )

    # optional regularized solver (lambda parameter):
    # lam = 1e-6
    # ATA = A_w.T @ A_w + lam * np.eye(N+1)
    # ATb = A_w.T @ f_w
    # coeffs = np.linalg.solve(ATA, ATb)

    # --- reconstruct and diagnostics ---
    recon = ( A @ _coeffs ) * nf

    # fig, ax = plt.subplots( 2, 1, figsize= ( 8, 8 ) )
    # ax[0].plot( beta, f, label='original' )
    # ax[0].plot( beta, recon, '--', label=f'recon N={N}' )
    # ax[0].plot( beta, f - recon, '--', label=f'Difference' )
    # ax[0].legend()
    # ax[0].set_xlabel( "beta" )
    # ax[0].set_ylabel( 'Fit' )
    # ax[0].grid(True)

    # ax[1].semilogy(np.abs(coeffs), 'o-')
    # ax[1].set_xlabel('order n')
    # ax[1].set_ylabel('|a_n|')
    # ax[1].set_title('Coefficients by order')
    # ax[1].grid(True)

    # plt.show()

    return recon, coeffs


def fit_residual_dGdt( folder_path: str, show_figs=False, is_square_ref=True, show_summary_fig=False ) -> None:
    # Define fit properties
    # Define boundary values
    fit_props               = FitProperties( )
    fit_props.dims          = 2
    fit_props.region_name   = "dGdt"
    fit_props.cheby_order_x = 15
    fit_props.cheby_order_y = 15
    fit_props.fcn_log_scale = False
    fit_props.x_hpatch_np   = 30
    fit_props.x_log_scale   = False
    fit_props.x_max         = 30.0
    fit_props.x_min         = 0.0
    fit_props.y_hpatch_np   = 4
    fit_props.y_log_scale   = True
    fit_props.y_max         = np.log10( 0.9998 )
    fit_props.y_min         = -4.0
    fit_props.cheby_abs_tol = 1E-4
    fit_props.cheby_rel_tol = 1E-4
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.max_ref_level = 10
    # fit_props.alpha_shift   = np.array( [ -1e-4, 0.0, 1e-4 ] )
    fit_props.alpha_shift   = np.array( [ 0.0 ] )
    # fit_props.alpha_shift   = np.array( [ 1e-3, 1e-2, 1e-1 ] )

    fit_props.num_x         = fit_props.cheby_order_x + 1
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = fit_props.cheby_order_y + 1
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1

    fit_props.generate_fitting_matrix( )

    # Load database
    fipath = r"D:\sergio\0050_OASIS_SM\TimeDomain\Gt.h5"
    with h5py.File( fipath, "r" ) as fid:
        mu          = fid[ "mu" ][:]
        beta        = fid[ "beta" ][:]
        data_raw    = fid[ "fcn" ][:]

    # Calculate residual Magee function
    target_fcn                  = dGdt_Magee
    interp_raw                  = sp.interpolate.RegularGridInterpolator( ( beta, mu ), data_raw )
    fit_props.magee_cheby_intv  = np.log10( np.array( [ mu[0], mu[-1] ] ) )
    fit_props.magee_cheby_np    = 1
    cheby_roots                 = sp.special.roots_chebyt( fit_props.magee_cheby_np )[0]
    fit_props.magee_mui_log     = np.array( [ fit_props.magee_cheby_intv[0] ] ) if fit_props.magee_cheby_np < 2 else map_to_interval( fit_props.magee_cheby_intv[0], fit_props.magee_cheby_intv[1], cheby_roots )
    fit_props.magee_mui         = 10**fit_props.magee_mui_log
    fit_props.magee_coeffs      = np.zeros( ( fit_props.alpha_shift.shape[0], fit_props.magee_cheby_np ) )

    data, data_norm, chfhpa     = calculate_residual_magee( fit_props, target_fcn, interp_raw, data_raw, beta, mu )
    B, M                        = np.meshgrid( beta, mu, indexing="ij" )

    # Plot residual function
    fig, ax = plt.subplots( 2, 2 )

    cnf = ax[0, 0].contourf( B, np.log10( M ), data_raw, levels=50 )
    plt.colorbar( cnf, ax=ax[0, 0] )

    cnf = ax[0, 1].contourf( B, np.log10( M ), data, levels=50 )
    plt.colorbar( cnf, ax=ax[0, 1] )

    fig, ax = plt.subplots( 2, 2, subplot_kw=dict(projection='3d') )

    cnf = ax[0, 0].plot_surface( B, np.log10( M ), data_raw, cmap=plt.cm.jet )
    plt.colorbar( cnf, ax=ax[0, 0] )

    cnf = ax[0, 1].plot_surface( B, np.log10( M ), data, cmap=plt.cm.jet)
    plt.colorbar( cnf, ax=ax[0, 1] )

    plt.show( )
    # raise ValueError( "Stop by user" )

    # Interpolate database
    fit_function_raw        = sp.interpolate.RegularGridInterpolator( ( beta, mu ), data_norm )
    fit_function            = lambda x, y: fit_function_raw( ( x, y ) )

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_2d( fit_function, fit_props, show_figs=show_figs ) )

    y_levels = np.array( [ -4.0, -3.0, -2.0, -1.0, np.log10( 0.9998 ) ] )
    for i in range( fit_props.x_hpatch_np ):
        for j in range( fit_props.y_hpatch_np ):
            fit_props_i         = copy.copy( ref_level.fit_props )
            fit_props_i.x_max   = ( i + 1 )
            fit_props_i.x_min   = i
            fit_props_i.y_max   = y_levels[j+1]
            fit_props_i.y_min   = y_levels[j]
            
            ref_level_i         = RefLevel( fit_props_i, parent=ref_level )

            # Create first fit to feed root refinement level
            ref_level_i.add_data( *fit_integral_2d( fit_function, fit_props_i, show_figs=show_figs ) )

            # Start adaptive refinement loop
            if not ref_level_i.check_tolerances( ) and ref_level_i.level < ref_level_i.fit_props.max_ref_level:
                fit_residual_2D_adaptive_interface( fit_function, ref_level_i, is_square_ref=is_square_ref, show_figs=show_figs )

            # Storage child data
            ref_level.child.append( ref_level_i )

    # Set starting position
    ref_level.set_start_index( 0 )

    # Plot results summary
    ref_level.show_summary( folder_path )

    print( ref_level )

    return ref_level, chfhpa


def fit_residual_dGdtx( folder_path: str, show_figs=False, is_square_ref=True, show_summary_fig=False ) -> None:
    # Define boundary values
    fit_props               = FitProperties( )
    fit_props.dims          = 2
    fit_props.region_name   = "dGdtx"
    fit_props.cheby_order_x = 30
    fit_props.cheby_order_y = 15
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = False
    fit_props.x_hpatch_np   = 30
    fit_props.x_max         = 30.0
    fit_props.x_min         = 0.0
    fit_props.y_hpatch_np   = 4
    fit_props.y_log_scale   = True
    fit_props.y_max         = np.log10( 0.9998 )
    fit_props.y_min         = -4.0
    fit_props.cheby_abs_tol = 1E-2
    fit_props.cheby_rel_tol = 1E-4
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.max_ref_level = 10
    # fit_props.alpha_shift   = np.array( [ -5e-4, 0.0, 5e-4 ] )
    # fit_props.alpha_shift   = np.array( [ -8e-3, 0.0, 8e-3 ] )
    fit_props.alpha_shift   = np.array( [ 1e-3, 1e-2, 1e-1 ] )
    # fit_props.alpha_shift   = np.array( [ 0.0 ] )

    fit_props.num_x         = fit_props.cheby_order_x + 1
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = fit_props.cheby_order_y + 1
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1

    fit_props.generate_fitting_matrix( )

    # Load database
    fipath = r"D:\sergio\0050_OASIS_SM\TimeDomain\Gtx.h5"
    with h5py.File( fipath, "r" ) as fid:
        mu          = fid[ "mu" ][:]
        beta        = fid[ "beta" ][:]
        data_raw    = fid[ "fcn" ][:]

    # Calculate residual Magee function
    target_fcn                  = dGdtx_Magee
    interp_raw                  = sp.interpolate.RegularGridInterpolator( ( beta, mu ), data_raw )
    fit_props.magee_cheby_intv  = np.log10( np.array( [ mu[0], mu[-1] ] ) )
    fit_props.magee_cheby_np    = 1
    cheby_roots                 = sp.special.roots_chebyt( fit_props.magee_cheby_np )[0]
    fit_props.magee_mui_log     = np.array( [ fit_props.magee_cheby_intv[0] ] ) if fit_props.magee_cheby_np < 2 else map_to_interval( fit_props.magee_cheby_intv[0], fit_props.magee_cheby_intv[1], cheby_roots )
    fit_props.magee_mui         = 10**fit_props.magee_mui_log
    fit_props.magee_coeffs      = np.zeros( ( fit_props.alpha_shift.shape[0], fit_props.magee_cheby_np ) )

    data, data_norm, chfhpa     = calculate_residual_magee( fit_props, target_fcn, interp_raw, data_raw, beta, mu )
    B, M                        = np.meshgrid( beta, mu, indexing="ij" )

    # Plot residual function
    fig, ax = plt.subplots( 2, 2 )

    cnf = ax[0, 0].contourf( B, np.log10( M ), data_raw, levels=50 )
    plt.colorbar( cnf, ax=ax[0, 0] )

    cnf = ax[0, 1].contourf( B, np.log10( M ), data, levels=50 )
    plt.colorbar( cnf, ax=ax[0, 1] )

    fig, ax = plt.subplots( 2, 2, subplot_kw=dict(projection='3d') )

    cnf = ax[0, 0].plot_surface( B, np.log10( M ), data_raw, cmap=plt.cm.jet )
    plt.colorbar( cnf, ax=ax[0, 0] )

    cnf = ax[0, 1].plot_surface( B, np.log10( M ), data, cmap=plt.cm.jet)
    plt.colorbar( cnf, ax=ax[0, 1] )

    plt.show( )

    # Interpolate database
    fit_function_raw        = sp.interpolate.RegularGridInterpolator( ( beta, mu ), data_norm )
    fit_function            = lambda x, y: fit_function_raw( ( x, y ) )

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_2d( fit_function, fit_props, show_figs=show_figs ) )

    y_levels = np.array( [ -4.0, -3.0, -2.0, -1.0, np.log10( 0.9998 ) ] )
    for i in range( fit_props.x_hpatch_np ):
        for j in range( fit_props.y_hpatch_np ):
            fit_props_i         = copy.copy( ref_level.fit_props )
            fit_props_i.x_max   = ( i + 1 )
            fit_props_i.x_min   = i
            fit_props_i.y_max   = y_levels[j+1]
            fit_props_i.y_min   = y_levels[j]
            
            ref_level_i         = RefLevel( fit_props_i, parent=ref_level )

            # Create first fit to feed root refinement level
            ref_level_i.add_data( *fit_integral_2d( fit_function, fit_props_i, show_figs=show_figs ) )

            # Start adaptive refinement loop
            if not ref_level_i.check_tolerances( ) and ref_level_i.level < ref_level_i.fit_props.max_ref_level:
                fit_residual_2D_adaptive_interface( fit_function, ref_level_i, is_square_ref=is_square_ref, show_figs=show_figs )

            # Storage child data
            ref_level.child.append( ref_level_i )

    # Set starting position
    ref_level.set_start_index( 0 )

    # Plot results summary
    ref_level.show_summary( folder_path )

    print( ref_level )

    return ref_level, chfhpa


def fit_residual_dGdtxx( folder_path: str, show_figs=False, is_square_ref=True, show_summary_fig=False ) -> None:
    # Define boundary values
    fit_props               = FitProperties( )
    fit_props.dims          = 2
    fit_props.region_name   = "dGdtxx"
    fit_props.cheby_order_x = 40
    fit_props.cheby_order_y = 15
    fit_props.fcn_log_scale = False
    fit_props.x_hpatch_np   = 30
    fit_props.x_log_scale   = False
    fit_props.x_max         = 30.0
    fit_props.x_min         = 0.0
    fit_props.y_hpatch_np   = 4
    fit_props.y_log_scale   = True
    fit_props.y_max         = np.log10( 0.9998 )
    fit_props.y_min         = -4.0
    fit_props.cheby_abs_tol = 1E1
    fit_props.cheby_rel_tol = 1E-4
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.max_ref_level = 10
    # fit_props.alpha_shift   = np.array( [ 0.0, 0.01, 0.05, 0.1 ] )
    fit_props.alpha_shift   = np.array( [ 1e-3, 1e-2, 1e-1 ] )

    fit_props.num_x         = fit_props.cheby_order_x + 1
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = fit_props.cheby_order_y + 1
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1

    fit_props.generate_fitting_matrix( )

    # Load database
    fipath = r"D:\sergio\0050_OASIS_SM\TimeDomain\Gtxx.h5"
    with h5py.File( fipath, "r" ) as fid:
        mu          = fid[ "mu" ][:]
        beta        = fid[ "beta" ][:]
        data_raw    = fid[ "fcn" ][:]

    # Calculate residual Magee function
    target_fcn                  = dGdtxx_Magee
    interp_raw                  = sp.interpolate.RegularGridInterpolator( ( beta, mu ), data_raw )
    fit_props.magee_cheby_intv  = np.log10( np.array( [ mu[0], mu[-1] ] ) )
    fit_props.magee_cheby_np    = 1
    cheby_roots                 = sp.special.roots_chebyt( fit_props.magee_cheby_np )[0]
    fit_props.magee_mui_log     = np.array( [ fit_props.magee_cheby_intv[0] ] ) if fit_props.magee_cheby_np < 2 else map_to_interval( fit_props.magee_cheby_intv[0], fit_props.magee_cheby_intv[1], cheby_roots )
    fit_props.magee_mui         = 10**fit_props.magee_mui_log
    fit_props.magee_coeffs      = np.zeros( ( fit_props.alpha_shift.shape[0], fit_props.magee_cheby_np ) )

    data, data_norm, chfhpa     = calculate_residual_magee( fit_props, target_fcn, interp_raw, data_raw, beta, mu )
    B, M                        = np.meshgrid( beta, mu, indexing="ij" )

    # Plot residual function
    fig, ax = plt.subplots( 2, 2 )

    cnf = ax[0, 0].contourf( B, np.log10( M ), data_raw, levels=50 )
    plt.colorbar( cnf, ax=ax[0, 0] )

    cnf = ax[0, 1].contourf( B, np.log10( M ), data, levels=50 )
    plt.colorbar( cnf, ax=ax[0, 1] )

    cnf = ax[1, 0].contourf( B, np.log10( M ), data_norm, levels=50 )
    plt.colorbar( cnf, ax=ax[1, 0] )


    fig, ax = plt.subplots( 2, 2, subplot_kw=dict(projection='3d') )

    cnf = ax[0, 0].plot_surface( B, np.log10( M ), data_raw, cmap=plt.cm.jet )
    plt.colorbar( cnf, ax=ax[0, 0] )

    cnf = ax[0, 1].plot_surface( B, np.log10( M ), data, cmap=plt.cm.jet)
    plt.colorbar( cnf, ax=ax[0, 1] )

    plt.show( )

    # Interpolate database
    fit_function_raw        = sp.interpolate.RegularGridInterpolator( ( beta, mu ), data_norm )
    fit_function            = lambda x, y: fit_function_raw( ( x, y ) )

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_2d( fit_function, fit_props, show_figs=show_figs ) )

    y_levels = np.array( [ -4.0, -3.0, -2.0, -1.0, np.log10( 0.9998 ) ] )
    for i in range( fit_props.x_hpatch_np ):
        for j in range( fit_props.y_hpatch_np ):
            fit_props_i         = copy.copy( ref_level.fit_props )
            fit_props_i.x_max   = ( i + 1 )
            fit_props_i.x_min   = i
            fit_props_i.y_max   = y_levels[ j + 1 ]
            fit_props_i.y_min   = y_levels[ j ]
            
            ref_level_i         = RefLevel( fit_props_i, parent=ref_level )

            # Create first fit to feed root refinement level
            ref_level_i.add_data( *fit_integral_2d( fit_function, fit_props_i, show_figs=show_figs ) )

            # Start adaptive refinement loop
            if not ref_level_i.check_tolerances( ) and ref_level_i.level < ref_level_i.fit_props.max_ref_level:
                fit_residual_2D_adaptive_interface( fit_function, ref_level_i, is_square_ref=is_square_ref, show_figs=show_figs )

            # Storage child data
            ref_level.child.append( ref_level_i )

    # Set starting position
    ref_level.set_start_index( 0 )

    # Plot results summary
    ref_level.show_summary( folder_path )

    print( ref_level )

    return ref_level, chfhpa


def fit_residual_dGdtt( folder_path: str, show_figs=False, is_square_ref=True, show_summary_fig=False ) -> None:
    # Define fit properties
    # Define boundary values
    fit_props               = FitProperties( )
    fit_props.dims          = 2
    fit_props.region_name   = "dGdtt"
    fit_props.cheby_order_x = 40
    fit_props.cheby_order_y = 15
    fit_props.fcn_log_scale = False
    fit_props.x_hpatch_np   = 30
    fit_props.x_log_scale   = False
    fit_props.x_max         = 19.0
    fit_props.x_min         = 0.0
    fit_props.y_hpatch_np   = 4
    fit_props.y_log_scale   = True
    fit_props.y_max         = np.log10( 0.9998 )
    fit_props.y_min         = -4.0
    fit_props.cheby_abs_tol = 1E-3
    fit_props.cheby_rel_tol = 1E-4
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.max_ref_level = 10
    # fit_props.alpha_shift   = np.array( [ 0.0, 1e-4 ] )
    fit_props.alpha_shift   = np.array( [ 1e-3, 1e-2, 1e-1 ] )
    # fit_props.alpha_shift   = np.array( [ 0.0 ] )

    fit_props.num_x         = fit_props.cheby_order_x + 1
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = fit_props.cheby_order_y + 1
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1

    fit_props.generate_fitting_matrix( )

    # Load database
    fipath = r"D:\sergio\0050_OASIS_SM\TimeDomain\Gtt.h5"
    with h5py.File( fipath, "r" ) as fid:
        mu          = fid[ "mu" ][:]
        beta        = fid[ "beta" ][:]
        data_raw    = fid[ "fcn" ][:]

    # Calculate residual Magee function
    target_fcn                  = dGdtt_Magee
    interp_raw                  = sp.interpolate.RegularGridInterpolator( ( beta, mu ), data_raw )
    fit_props.magee_cheby_intv  = np.log10( np.array( [ mu[0], mu[-1] ] ) )
    fit_props.magee_cheby_np    = 1
    # cheby_roots                 = sp.special.roots_chebyt( fit_props.magee_cheby_np )[0]
    # cheby_roots                 = sp.special.roots_chebyt( fit_props.magee_cheby_np )[0]
    cheby_roots                 = np.linspace( -1.0, 1.0, fit_props.magee_cheby_np )
    fit_props.magee_mui_log     = np.array( [ fit_props.magee_cheby_intv[0] ] ) if fit_props.magee_cheby_np < 2 else map_to_interval( fit_props.magee_cheby_intv[0], fit_props.magee_cheby_intv[1], cheby_roots )
    fit_props.magee_mui         = 10**fit_props.magee_mui_log
    fit_props.magee_coeffs      = np.zeros( ( fit_props.alpha_shift.shape[0], fit_props.magee_cheby_np ) )

    data, data_norm, chfhpa     = calculate_residual_magee( fit_props, target_fcn, interp_raw, data_raw, beta, mu )
    B, M                        = np.meshgrid( beta, mu, indexing="ij" )

    # Plot residual function
    fig, ax                     = plt.subplots( 2, 2 )
    fig.suptitle( f"Normalization factor: {fit_props.magee_norm_f}" )

    cnf = ax[0, 0].contourf( B, np.log10( M ), data_raw, levels=50 )
    plt.colorbar( cnf, ax=ax[0, 0] )
    # add_hover_tool( B, np.log10( M ), data_raw, cnf )

    cnf = ax[0, 1].contourf( B, np.log10( M ), data, levels=50 )
    plt.colorbar( cnf, ax=ax[0, 1] )
    # add_hover_tool( B, np.log10( M ), data, cnf )

    mplcursors.cursor( )

    fig, ax = plt.subplots( 2, 2, subplot_kw=dict(projection='3d') )

    cnf = ax[0, 0].plot_surface( B, np.log10( M ), data_raw, cmap=plt.cm.jet )
    plt.colorbar( cnf, ax=ax[0, 0] )

    cnf = ax[0, 1].plot_surface( B, np.log10( M ), data, cmap=plt.cm.jet)
    plt.colorbar( cnf, ax=ax[0, 1] )

    plt.show( )

    # Interpolate database
    fit_function_raw        = sp.interpolate.RegularGridInterpolator( ( beta, mu ), data_norm )
    fit_function            = lambda x, y: fit_function_raw( ( x, y ) )

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_2d( fit_function, fit_props, show_figs=show_figs ) )

    y_levels = np.array( [ -4.0, -3.0, -2.0, -1.0, np.log10( 0.9998 ) ] )
    for i in range( fit_props.x_hpatch_np ):
        for j in range( fit_props.y_hpatch_np ):
            fit_props_i         = copy.copy( ref_level.fit_props )
            fit_props_i.x_max   = ( i + 1 )
            fit_props_i.x_min   = i
            fit_props_i.y_max   = y_levels[j+1]
            fit_props_i.y_min   = y_levels[j]
            
            ref_level_i         = RefLevel( fit_props_i, parent=ref_level )

            # Create first fit to feed root refinement level
            ref_level_i.add_data( *fit_integral_2d( fit_function, fit_props_i, show_figs=show_figs ) )

            # Start adaptive refinement loop
            if not ref_level_i.check_tolerances( ) and ref_level_i.level < ref_level_i.fit_props.max_ref_level:
                fit_residual_2D_adaptive_interface( fit_function, ref_level_i, is_square_ref=is_square_ref, show_figs=show_figs )

            # Storage child data
            ref_level.child.append( ref_level_i )

    # Set starting position
    ref_level.set_start_index( 0 )

    # Plot results summary
    ref_level.show_summary( folder_path )

    print( ref_level )

    return ref_level, chfhpa


def fit_residual_dGdttx( folder_path: str, show_figs=False, is_square_ref=True, show_summary_fig=False ) -> None:
    # Define boundary values
    fit_props               = FitProperties( )
    fit_props.dims          = 2
    fit_props.region_name   = "dGdttx"
    fit_props.cheby_order_x = 40
    fit_props.cheby_order_y = 15
    fit_props.fcn_log_scale = False
    fit_props.x_hpatch_np   = 30
    fit_props.x_log_scale   = False
    fit_props.x_max         = 30.0
    fit_props.x_min         = 0.0
    fit_props.y_hpatch_np   = 4
    fit_props.y_log_scale   = True
    fit_props.y_max         = np.log10( 0.999 )
    fit_props.y_min         = -4.0
    fit_props.cheby_abs_tol = 1E0
    fit_props.cheby_rel_tol = 1E-4
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.max_ref_level = 10
    # fit_props.alpha_shift   = np.array( [ 1e-4, 0.0, 1e-2, 1e-1 ] )
    fit_props.alpha_shift   = np.array( [ 1e-3, 1e-2, 1e-1 ] )

    fit_props.num_x         = fit_props.cheby_order_x + 1
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = fit_props.cheby_order_y + 1
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1

    fit_props.generate_fitting_matrix( )

    # Load database
    fipath = r"D:\sergio\0050_OASIS_SM\TimeDomain\Gttx.h5"
    with h5py.File( fipath, "r" ) as fid:
        mu          = fid[ "mu" ][:]
        beta        = fid[ "beta" ][:]
        data_raw    = fid[ "fcn" ][:]

    # Calculate residual Magee function
    target_fcn                  = dGdttx_Magee
    interp_raw                  = sp.interpolate.RegularGridInterpolator( ( beta, mu ), data_raw )
    fit_props.magee_cheby_intv  = np.log10( np.array( [ mu[0], mu[-1] ] ) )
    fit_props.magee_cheby_np    = 100
    cheby_roots                 = sp.special.roots_chebyt( fit_props.magee_cheby_np )[0]
    fit_props.magee_mui_log     = np.array( [ fit_props.magee_cheby_intv[0] ] ) if fit_props.magee_cheby_np < 2 else map_to_interval( fit_props.magee_cheby_intv[0], fit_props.magee_cheby_intv[1], cheby_roots )
    fit_props.magee_mui         = 10**fit_props.magee_mui_log
    fit_props.magee_coeffs      = np.zeros( ( fit_props.alpha_shift.shape[0], fit_props.magee_cheby_np ) )

    data, data_norm, chfhpa     = calculate_residual_magee( fit_props, target_fcn, interp_raw, data_raw, beta, mu )
    B, M                        = np.meshgrid( beta, mu, indexing="ij" )

    # Plot residual function
    fig, ax = plt.subplots( 2, 2 )

    cnf = ax[0, 0].contourf( B, np.log10( M ), data_raw, levels=50 )
    plt.colorbar( cnf, ax=ax[0, 0] )

    cnf = ax[0, 1].contourf( B, np.log10( M ), data, levels=50 )
    plt.colorbar( cnf, ax=ax[0, 1] )

    fig, ax = plt.subplots( 2, 2, subplot_kw=dict(projection='3d') )

    cnf = ax[0, 0].plot_surface( B, np.log10( M ), data_raw, cmap=plt.cm.jet )
    plt.colorbar( cnf, ax=ax[0, 0] )

    cnf = ax[0, 1].plot_surface( B, np.log10( M ), data, cmap=plt.cm.jet)
    plt.colorbar( cnf, ax=ax[0, 1] )

    plt.show( )

    # Interpolate database
    fit_function_raw        = sp.interpolate.RegularGridInterpolator( ( beta, mu ), data_norm )
    fit_function            = lambda x, y: fit_function_raw( ( x, y ) )

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_2d( fit_function, fit_props, show_figs=show_figs ) )

    y_levels = np.array( [ -4.0, -3.0, -2.0, -1.0, np.log10( 0.9998 ) ] )
    for i in range( fit_props.x_hpatch_np ):
        for j in range( fit_props.y_hpatch_np ):
            fit_props_i         = copy.copy( ref_level.fit_props )
            fit_props_i.x_max   = ( i + 1 )
            fit_props_i.x_min   = i
            fit_props_i.y_max   = y_levels[ j + 1 ]
            fit_props_i.y_min   = y_levels[ j ]
            
            ref_level_i         = RefLevel( fit_props_i, parent=ref_level )

            # Create first fit to feed root refinement level
            ref_level_i.add_data( *fit_integral_2d( fit_function, fit_props_i, show_figs=show_figs ) )

            # Start adaptive refinement loop
            if not ref_level_i.check_tolerances( ) and ref_level_i.level < ref_level_i.fit_props.max_ref_level:
                fit_residual_2D_adaptive_interface( fit_function, ref_level_i, is_square_ref=is_square_ref, show_figs=show_figs )

            # Storage child data
            ref_level.child.append( ref_level_i )

    # Set starting position
    ref_level.set_start_index( 0 )

    # Plot results summary
    ref_level.show_summary( folder_path )

    print( ref_level )

    return ref_level, chfhpa


def fit_residual_dGdttxx( folder_path: str, show_figs=True, is_square_ref=True, show_summary_fig=False ) -> None:
    # Define boundary values
    fit_props               = FitProperties( )
    fit_props.dims          = 2
    fit_props.region_name   = "dGdttxx"
    fit_props.cheby_order_x = 30
    fit_props.cheby_order_y = 30
    fit_props.fcn_log_scale = False
    fit_props.x_log_scale   = False
    fit_props.x_hpatch_np   = 30
    fit_props.x_max         = 30.0
    fit_props.x_min         = 0.0
    fit_props.y_log_scale   = True
    fit_props.y_max         = np.log10( 0.999 )
    fit_props.y_hpatch_np   = 4
    fit_props.y_min         = -4.0
    fit_props.cheby_abs_tol = 1E+1
    fit_props.cheby_rel_tol = 1E-3
    fit_props.x_map_fcn     = lambda x: x
    fit_props.y_map_fcn     = lambda y: y
    fit_props.max_ref_level = 10
    # fit_props.alpha_shift   = np.array( [ 1e-4, 1e-2, 1e-1 ] )
    # fit_props.alpha_shift   = np.array( [ 1e-4, -3e-2, 5e-1 ] )
    # fit_props.alpha_shift   = np.array( [ 1e-4, -3e-2 ] )
    fit_props.alpha_shift   = np.array( [ 1e-3, 1e-2, 1e-1 ] )
    # fit_props.alpha_shift   = np.array( [ 0.0 ] )

    fit_props.num_x         = fit_props.cheby_order_x + 1
    fit_props.num_x_fit     = fit_props.cheby_order_x + 1
    fit_props.num_y         = fit_props.cheby_order_y + 1
    fit_props.num_y_fit     = fit_props.cheby_order_y + 1

    fit_props.generate_fitting_matrix( )

    # Load database
    fipath = r"D:\sergio\0050_OASIS_SM\TimeDomain\Gttxx.h5"
    with h5py.File( fipath, "r" ) as fid:
        mu          = fid[ "mu" ][:]
        beta        = fid[ "beta" ][:]
        data_raw    = fid[ "fcn" ][:]

    pos         = beta < fit_props.x_max+1e-6
    beta        = beta[ pos ]
    data_raw    = data_raw[ pos, : ]

    # Calculate residual Magee function
    target_fcn                  = dGdttxx_Magee
    interp_raw                  = sp.interpolate.RegularGridInterpolator( ( beta, mu ), data_raw )
    fit_props.magee_cheby_intv  = np.log10( np.array( [ mu[0], mu[-1] ] ) )
    fit_props.magee_cheby_np    = 100
    # cheby_roots                 = sp.special.roots_chebyt( fit_props.magee_cheby_np )[0]
    cheby_roots                 = np.linspace( -1.0, 1.0, fit_props.magee_cheby_np )
    fit_props.magee_mui_log     = np.array( [ fit_props.magee_cheby_intv[0] ] ) if fit_props.magee_cheby_np < 2 else map_to_interval( fit_props.magee_cheby_intv[0], fit_props.magee_cheby_intv[1], cheby_roots )
    fit_props.magee_mui         = 10**fit_props.magee_mui_log
    fit_props.magee_coeffs      = np.zeros( ( fit_props.alpha_shift.shape[0], fit_props.magee_cheby_np ) )
    
    data, data_norm, chfhpa     = calculate_residual_magee( fit_props, target_fcn, interp_raw, data_raw, beta, mu )
    B, M                        = np.meshgrid( beta, mu, indexing="ij" )

    # Plot residual function
    fig, ax = plt.subplots( 2, 2 )

    cnf = ax[0, 0].contourf( B, np.log10( M ), data_raw, levels=50 )
    plt.colorbar( cnf, ax=ax[0, 0] )

    cnf = ax[0, 1].contourf( B, np.log10( M ), data, levels=50 )
    plt.colorbar( cnf, ax=ax[0, 1] )

    cnf = ax[1, 0].contourf( B, np.log10( M ), np.log10( np.abs( data_raw ) ), levels=50 )
    plt.colorbar( cnf, ax=ax[1, 0] )

    cnf = ax[1, 1].contourf( B, np.log10( M ), np.log10( np.abs( data ) ), levels=50 )
    plt.colorbar( cnf, ax=ax[1, 1] )

    fig, ax = plt.subplots( 2, 2, subplot_kw=dict(projection='3d') )

    cnf = ax[0, 0].plot_surface( B, np.log10( M ), data_raw, cmap=plt.cm.jet )
    plt.colorbar( cnf, ax=ax[0, 0] )

    cnf = ax[0, 1].plot_surface( B, np.log10( M ), data, cmap=plt.cm.jet)
    plt.colorbar( cnf, ax=ax[0, 1] )

    plt.show( )

    # Interpolate database
    fit_function_raw        = sp.interpolate.RegularGridInterpolator( ( beta, mu ), data_norm )
    fit_function            = lambda x, y: fit_function_raw( ( x, y ) )

    # Create root FitRegion
    ref_level               = RefLevel( copy.copy( fit_props ) )

    # Create first fit to feed root refinement level
    ref_level.add_data( *fit_integral_2d( fit_function, fit_props, show_figs=show_figs ) )

    y_levels = np.array( [ -4.0, -3.0, -2.0, -1.0, np.log10( 0.9998 ) ] )
    for i in range( fit_props.x_hpatch_np ):
        for j in range( fit_props.y_hpatch_np ):
            fit_props_i         = copy.copy( ref_level.fit_props )
            fit_props_i.x_max   = ( i + 1 )
            fit_props_i.x_min   = i
            fit_props_i.y_max   = y_levels[ j + 1 ]
            fit_props_i.y_min   = y_levels[ j ]
            
            ref_level_i         = RefLevel( fit_props_i, parent=ref_level )

            # Create first fit to feed root refinement level
            ref_level_i.add_data( *fit_integral_2d( fit_function, fit_props_i, show_figs=show_figs ) )

            # Start adaptive refinement loop
            if not ref_level_i.check_tolerances( ) and ref_level_i.level < ref_level_i.fit_props.max_ref_level:
                fit_residual_2D_adaptive_interface( fit_function, ref_level_i, is_square_ref=is_square_ref, show_figs=show_figs )

            # Storage child data
            ref_level.child.append( ref_level_i )

    # Set starting position
    ref_level.set_start_index( 0 )

    # Plot results summary
    ref_level.show_summary( folder_path )

    print( ref_level )

    return ref_level, chfhpa


def generate_dGdt(folder_path: str, show_summary_fig=False, show_figs=False)->None:
    # Fit coefficients
    fit_region, chfhpa = fit_residual_dGdt(folder_path, show_summary_fig=show_summary_fig, show_figs=show_figs)

    # Write coefficients
    write_coeffs_module_adaptive_2d_only_header( fit_region, folder_path, "dGdt", is_time=True )

    for i,mi in enumerate( chfhpa.get_interpolators( ) ):
        write_coeffs_module_adaptive_1d_only_header( mi.get_fit_fcn( ), folder_path, f"dGdtA{i:d}" )


def generate_dGdtx(folder_path: str, show_summary_fig=False, show_figs=False)->None:
    # Fit coefficients
    fit_region, chfhpa = fit_residual_dGdtx(folder_path, show_summary_fig=show_summary_fig, show_figs=show_figs)

    # Write coefficients
    write_coeffs_module_adaptive_2d_only_header( fit_region, folder_path, "dGdtx", is_time=True )

    for i,mi in enumerate( chfhpa.get_interpolators( ) ):
        write_coeffs_module_adaptive_1d_only_header( mi.get_fit_fcn( ), folder_path, f"dGdtxA{i:d}" )


def generate_dGdtxx(folder_path: str, show_summary_fig=False, show_figs=False)->None:
    # Fit coefficients
    fit_region, chfhpa = fit_residual_dGdtxx(folder_path, show_summary_fig=show_summary_fig, show_figs=show_figs)

    # Write coefficients
    write_coeffs_module_adaptive_2d_only_header( fit_region, folder_path, "dGdtxx", is_time=True )

    for i,mi in enumerate( chfhpa.get_interpolators( ) ):
        write_coeffs_module_adaptive_1d_only_header( mi.get_fit_fcn( ), folder_path, f"dGdtxxA{i:d}" )


def generate_dGdtt(folder_path: str, show_summary_fig=False, show_figs=False)->None:
    # Fit coefficients
    fit_region, chfhpa = fit_residual_dGdtt(folder_path, show_summary_fig=show_summary_fig, show_figs=show_figs)

    # Write coefficients
    write_coeffs_module_adaptive_2d_only_header( fit_region, folder_path, "dGdtt", is_time=True )

    for i,mi in enumerate( chfhpa.get_interpolators( ) ):
        write_coeffs_module_adaptive_1d_only_header( mi.get_fit_fcn( ), folder_path, f"dGdttA{i:d}" )


def generate_dGdttx(folder_path: str, show_summary_fig=False, show_figs=False)->None:
    # Fit coefficients
    fit_region, chfhpa = fit_residual_dGdttx(folder_path, show_summary_fig=show_summary_fig, show_figs=show_figs)

    # Write coefficients
    write_coeffs_module_adaptive_2d_only_header( fit_region, folder_path, "dGdttx", is_time=True )

    for i,mi in enumerate( chfhpa.get_interpolators( ) ):
        write_coeffs_module_adaptive_1d_only_header( mi.get_fit_fcn( ), folder_path, f"dGdttxA{i:d}" )


def generate_dGdttxx(folder_path: str, show_summary_fig=False, show_figs=False)->None:
    # Fit coefficients
    fit_region, chfhpa = fit_residual_dGdttxx(folder_path, show_summary_fig=show_summary_fig, show_figs=show_figs)

    # Write coefficients
    write_coeffs_module_adaptive_2d_only_header( fit_region, folder_path, "dGdttxx", is_time=True )

    for i,mi in enumerate( chfhpa.get_interpolators( ) ):
        write_coeffs_module_adaptive_1d_only_header( mi.get_fit_fcn( ), folder_path, f"dGdttxxA{i:d}" )


def main(  ) -> None:
    # Load database
    fipath = r"D:\sergio\0050_OASIS_SM\TimeDomain\test_fcn_data.h5"
    with h5py.File( fipath, "r" ) as fid:
        mu      = fid[ "mu" ][:]
        beta    = fid[ "beta" ][:]
        data    = fid[ "fcn" ][:]

    B, M = np.meshgrid( beta, mu, indexing="ij" )

    # Interpolate database
    interp = sp.interpolate.RegularGridInterpolator( ( beta, mu ), data, method="cubic" )

    print( mu )
    print( beta )
    print( interp( ( 19.0, 0.5 ) ) )

    # Plot
    plt.contourf( B, M, interp( ( B, M ) ) - data )
    plt.colorbar( )
    plt.show( )


def G_Magee( beta: np.ndarray, mu ) -> np.ndarray:
    lt = np.pi * beta**2.0 / 8 / np.sqrt( 2 ) * np.exp( -beta**2.0 * mu / 4.0 )
    return lt * (
                    sp.special.jv( 0.25, beta**2.0 / 8.0 )
                    *
                    sp.special.jv( 0.75, beta**2.0 / 8.0 )
                    -
                    sp.special.jv( -0.25, beta**2.0 /8.0 )
                    *
                    sp.special.jv( -0.75, beta**2.0 / 8.0 )
                )


def dGdt_Magee( beta: np.ndarray, mu: float, alpha=0.0 ) -> np.ndarray:
    eta     = np.sqrt( 1 - mu**2.0 )
    eta     = eta if eta > 0.5 else 0.5
    lt      = np.pi * beta**3.0 / 16 / np.sqrt( 2 ) * np.exp( -beta**2.0 * mu / 4.0 )
    x2      = beta**2.0 / 8.0
    bt      = (
                    sp.special.jv( 0.25, x2 + alpha )
                    *
                    sp.special.jv( -0.25, x2 + alpha )
                    +
                    sp.special.jv( 0.75, x2 + alpha )
                    *
                    sp.special.jv( -0.75, x2 + alpha )
                )

    pos     = np.isnan( bt )
    bt[pos] = 0.0
    pos     = beta < 1e-1
    bt[pos] = 0.0

    return lt * bt


def dGdtx_Magee( beta: np.ndarray, mu: float, alpha=0.0 ) -> np.ndarray:
    return ( -beta**2.0 / 4.0 ) * dGdt_Magee( beta, mu, alpha=alpha )


def dGdtx_Magee2( beta: np.ndarray, mu: float, alpha=0.0 ) -> np.ndarray:
    eps = 1e-6
    return ( dGdt_Magee( beta, mu+eps, alpha=alpha ) - dGdt_Magee( beta, mu-eps, alpha=alpha ) ) / 2.0 / eps


def dGdtxx_Magee( beta: np.ndarray, mu: float, alpha=0.0 ) -> np.ndarray:
    return ( -beta**2.0 / 4.0 ) * dGdtx_Magee( beta, mu, alpha=alpha )


def dGdtt_Magee( beta: np.ndarray, mu: float, alpha=0.0 ) -> np.ndarray:
    # Calculate leading term
    lt      = np.pi / 16 / np.sqrt( 2 )
    expt    = np.exp( -beta**2.0 * mu / 4.0 )
    r       = 1.0
    x       = beta**2.0 / 8.0 + alpha
    pos     = x < 1e-6
    x[pos]  = 0.0

    # Calculate term 1
    a       = beta**3.0
    y0      = (
                    sp.special.jv( 0.25, x )
                    *
                    sp.special.jv( -0.25, x )
                    +
                    sp.special.jv( 0.75, x )
                    *
                    sp.special.jv( -0.75, x )
                )
    
    # plt.plot( beta, y0 )
    pos     = np.isnan( y0 )
    # print( y0[:5] )
    # print( pos )
    if pos.sum( ) > 0:
        num_pos = np.where( pos )[0][-1] + 1
        # y0[pos] = 0.0
        # pos     = beta < 1e-1
        # y0[pos] = 0.0
        # y0[0] = y0[3]
        # y0[1] = y0[3]
        y0[:num_pos] = y0[num_pos]
    # plt.plot( beta, y0 )
    # plt.show( )
    
    T1      = - 0.5 * a * y0 * beta * mu * expt

    # Calculate term 2
    
    a       = 3.0 * beta**2.0
    b       = beta**4.0 / 4.0
    c       = b / 2.0

    y0      = a * (
                        jv( 0.25, x ) * jv( -0.25, x )
                        +
                        jv( 0.75, x ) * jv( -0.75, x )
                )
    
    y1      = b * (
                        jv( -0.75, x ) * jv( -0.25, x )
                        -
                        jv( 0.75, x ) * jv( 0.25, x )
                    )
    
    y2      = c * (
                        - jv( 1.25, x ) * jv( -0.25, x )
                        +
                        jv( -1.25, x ) * jv( 0.25, x )
                    )
    
    y3      = c * (
                        - jv( 1.75, x ) * jv( -0.75, x )
                        +
                        jv( -1.75, x ) * jv( 0.75, x )
                    )
        
    yt      = ( y0 + y1 + y2 + y3 )


    # plt.plot( beta, yt )
    pos     = np.isnan( yt )
    yt[pos] = 0.0
    dt      = 5e-2
    pos     = beta < dt
    num_pos = np.where( pos )[0][-1]
    yt[pos] = yt[num_pos] / dt * beta[pos]
    # plt.plot( beta, yt )
    # plt.show( )
    
    T2      = yt * expt

    # plt.plot( lt * ( T1 + T2 ) )
    # plt.show( )

    return lt * ( T1 + T2 )


def dGdttx_Magee( beta: np.ndarray, mu: float, alpha=0.0 ) -> np.ndarray:
    return ( -beta**2.0 / 4.0 ) * dGdtt_Magee( beta, mu, alpha=alpha )


def dGdttxx_Magee( beta: np.ndarray, mu: float, alpha=0.0 ) -> np.ndarray:
    return ( -beta**2.0 / 4.0 ) * dGdttx_Magee( beta, mu, alpha=alpha )


def G0_Magee( beta: np.ndarray, mu ) -> np.ndarray:
    lt = np.pi * beta**3.0 / 16 / np.sqrt( 2 ) * np.exp( -beta**2.0 * mu / 4.0 )
    return lt * (
                    sp.special.jv( 0.25, beta**2.0 / 8.0 )
                    *
                    sp.special.jv( -0.25, beta**2.0 / 8.0 )
                    +
                    sp.special.jv( 0.75, beta**2.0 /8.0 )
                    *
                    sp.special.jv( -0.75, beta**2.0 / 8.0 )
                )


def G0_Magee_series( beta: np.ndarray, mu: float, n: int ) -> np.ndarray:
    lt  = np.pi * beta**3.0 / 16 / np.sqrt( 2 ) * np.exp( -beta**2.0 * mu / 4.0 )
    G   = lt * (
                    sp.special.jv( 0.25 * ( n + 1 ), beta**2.0 / 8.0 )
                    *
                    sp.special.jv( -0.25 * ( n + 1 ), beta**2.0 / 8.0 )
                    +
                    sp.special.jv( 0.75 * ( n + 1 ), beta**2.0 /8.0 )
                    *
                    sp.special.jv( -0.75 * ( n + 1 ), beta**2.0 / 8.0 )
                )
    pos     = beta < 1e-6
    G[pos]  = 0.0

    return G


def G0_Magee_dx( beta: np.ndarray, mu ) -> np.ndarray:
    lt      = - np.pi * beta**5.0 / 64 / np.sqrt( 2 ) * np.exp( -beta**2.0 * mu / 4.0 )
    G       = lt * (
                        sp.special.jv( 0.25, beta**2.0 / 8.0 )
                        *
                        sp.special.jv( -0.25, beta**2.0 / 8.0 )
                        +
                        sp.special.jv( 0.75, beta**2.0 /8.0 )
                        *
                        sp.special.jv( -0.75, beta**2.0 / 8.0 )
                    )

    return G


def G0_Magee_dx_series( beta: np.ndarray, mu: float, n: int ) -> np.ndarray:
    lt      = - np.pi * beta**5.0 / 64 / np.sqrt( 2 ) * np.exp( -beta**2.0 * mu / 4.0 )
    xarg    = np.sqrt( 1 - mu**2.0 ) * beta**2.0 / 8.0
    G       = lt * (
                        sp.special.jv( 0.25 * ( n + 1 ), xarg )
                        *
                        sp.special.jv( -0.25 * ( n + 1 ), xarg )
                        +
                        sp.special.jv( 0.75 * ( n + 1 ), xarg )
                        *
                        sp.special.jv( -0.75 * ( n + 1 ), xarg )
                    )
    
    pos     = beta < 1e-6
    G[pos]  = 0.0

    return G


def G0_Magee_dxx( beta: np.ndarray, mu ) -> np.ndarray:
    lt = np.pi * beta**7.0 / 256 / np.sqrt( 2 ) * np.exp( -beta**2.0 * mu / 4.0 )
    return lt * (
                    sp.special.jv( 0.25, beta**2.0 / 8.0 )
                    *
                    sp.special.jv( -0.25, beta**2.0 / 8.0 )
                    +
                    sp.special.jv( 0.75, beta**2.0 /8.0 )
                    *
                    sp.special.jv( -0.75, beta**2.0 / 8.0 )
                )


def G0_Magee_dxx_series( beta: np.ndarray, mu: float, n: int ) -> np.ndarray:
    lt = np.pi * beta**7.0 / 256 / np.sqrt( 2 ) * np.exp( -beta**2.0 * mu / 4.0 )
    return lt * (
                    sp.special.jv( 0.25 * ( n + 1 ), beta**2.0 / 8.0 )
                    *
                    sp.special.jv( -0.25 * ( n + 1 ), beta**2.0 / 8.0 )
                    +
                    sp.special.jv( 0.75 * ( n + 1 ), beta**2.0 /8.0 )
                    *
                    sp.special.jv( -0.75 * ( n + 1 ), beta**2.0 / 8.0 )
                )


def G0_Magee_dt( beta: np.ndarray, mu ) -> np.ndarray:
    # Calculate leading term
    lt  = np.pi / 16 / np.sqrt( 2 )

    # Calculate term 1
    a   = beta**3.0
    y0  = (
                sp.special.jv( 0.25, beta**2.0 / 8.0 )
                *
                sp.special.jv( -0.25, beta**2.0 / 8.0 )
                +
                sp.special.jv( 0.75, beta**2.0 /8.0 )
                *
                sp.special.jv( -0.75, beta**2.0 / 8.0 )
            )
    
    T1  = - 0.5 * a * y0 * beta * mu * np.exp( -beta**2.0 * mu / 4.0 )

    # Calculate term 2
    x   = beta**2.0 / 8.0
    a   = 3.0 * beta**2.0
    b   = beta**4.0 / 4.0
    c   = b / 2.0

    y0  = a * (
                    jv( 0.25, x ) * jv( -0.25, x )
                    +
                    jv( 0.75, x ) * jv( -0.75, x )
            )
    
    y1  = b * (
                    jv( -0.75, x ) * jv( -0.25, x )
                    -
                    jv( 0.75, x ) * jv( 0.25, x )
                )
    
    y2  = c * (
                    - jv( 1.25, x ) * jv( -0.25, x )
                    +
                    jv( -1.25, x ) * jv( 0.25, x )
                )
    
    y3  = c * (
                    - jv( 1.75, x ) * jv( -0.75, x )
                    +
                    jv( -1.75, x ) * jv( 0.75, x )
                )
    
    T2  = ( y0 + y1 + y2 + y3 ) * np.exp( -beta**2.0 * mu / 4.0 )

    return lt * ( T1 + T2 )


def map_to_interval( a: float, b: float, t: np.ndarray ) -> None:
    return 0.5 * ( ( 1 - t ) * a + ( t + 1 ) * b )

def map_to_unit( a: float, b: float, x: np.ndarray ) -> None:
    return ( 2.0 * x - ( a + b ) ) / ( b - a )




def plot_database( ) -> None:
    # x = np.linspace( 0, 20, 500 )
    # y = np.linspace( 0, 20, 500 )
    # X, Y = np.meshgrid( x, y )
    c = 1.0
    b = 6.0
    d = 5.0
    # F = np.exp( -c * ( ( X / b )**2.0 ) * np.exp( - c * ( Y * d ) **2.0 ) )
    # F = np.exp( -c * ( ( X / b )**2.0 ) ) * np.exp( - c * ( Y * d ) **2.0 )
    # F = np.exp( -c * ( ( X / b )**2.0 ) )
    # cnf = plt.contourf( X, Y, F, levels=50, cmap="jet" )
    # plt.colorbar( )
    # plt.show( )
    # raise ValueError( "Stop by user" )

    fipath = r"D:\sergio\0050_OASIS_SM\TimeDomain\Gtxx.h5"
    with h5py.File( fipath, "r" ) as fid:
        mu          = fid[ "mu" ][:]
        beta        = fid[ "beta" ][:]
        data_raw    = fid[ "fcn" ][:]

    pos_max     = updown_cross( beta, data_raw[:, 0] )
    beta_max    = beta[ pos_max ]
    fcn_theory  = G0_Magee_dx

    # fig, ax = plt.subplots( 3, 2 )

    # mu_colors = plt.cm.jet( np.linspace( 0, 1, mu.shape[0] ) )
    # for i in range( mu.shape[0] ):
    #     ax[0, 0].plot( beta, data_raw[:, i], color=mu_colors[i] )
    # ax[0, 0].set_yscale( "log" )

    # beta_colors = plt.cm.jet( np.linspace( 0, 1, beta_max.shape[0] ) )
    # for i in range( beta_max.shape[0] ):
    #     ax[0, 1].plot( mu, np.abs( data_raw[pos_max[i], :] ), color=beta_colors[i] )
    # ax[0, 1].set_yscale( "log" )

    # mu_colors = plt.cm.jet( np.linspace( 0, 1, mu_colors.shape[0] ) )
    # for i in range( mu.shape[0] ):
    #     dx = fcn_theory( beta, mu[i] )
    #     ax[1, 0].plot( beta, data_raw[:, i], color=mu_colors[i] )
    #     ax[1, 0].plot( beta, dx, "-x", color=mu_colors[-(i+1)] )
    # ax[1, 0].set_yscale( "log" )

    # beta_colors = plt.cm.jet( np.linspace( 0, 1, beta_max.shape[0] ) )
    # for i in range( beta_max.shape[0] ):
    #     dx = fcn_theory( beta_max[i], mu )
    #     ax[1, 1].plot( mu, np.abs( data_raw[pos_max[i], :] ), color=beta_colors[i] )
    #     ax[1, 1].plot( mu, np.abs( dx ), "-x", color=beta_colors[-(i+1)] )
    # ax[1, 1].set_yscale( "log" )
    
    # mu_colors = plt.cm.jet( np.linspace( 0, 1, mu_colors.shape[0] ) )
    # for i in range( mu.shape[0] ):
    #     dx = fcn_theory( beta, mu[i] )
    #     ax[2, 0].plot( beta, np.abs( data_raw[:, i] ) - np.abs( dx ), color=mu_colors[i] )
    # ax[2, 0].set_yscale( "log" )

    # beta_colors = plt.cm.jet( np.linspace( 0, 1, beta_max.shape[0] ) )
    # for i in range( beta_max.shape[0] ):
    #     dx = fcn_theory( beta_max[i], mu )
    #     ax[2, 1].plot( mu, np.abs( data_raw[pos_max[i], :] ) - np.abs( dx ), color=beta_colors[i] )
    # ax[2, 1].set_yscale( "log" )
    
    # plt.show( )
    # raise ValueError( "Stop by user" )

    # alpha                   = 1.00001
    # f                       = np.exp( ( (-1.0+(alpha-mu)) / ( alpha - mu )**12.0 ) )
    # plt.plot( mu, f )
    # plt.show( )
    # raise ValueError( "Stop by user" )

    # F = 1.0 - np.exp( -c * ( ( beta / b )**10) )
    # plt.plot( beta, F )
    # plt.show( )
    # raise ValueError( "Stop by user" )

    B, M = np.meshgrid( beta, mu, indexing="ij" )
    cnf = plt.contourf( B, M, data_raw )
    plt.colorbar( cnf )
    plt.show( )

    data = data_raw.copy( )
    for i in range( mu.shape[0] ):
        G0                      = G0_Magee_dxx( beta, mu[i] )
        G0[0]                   = 0.0

        # plt.plot( beta, data_raw[:, i] )
        # plt.plot( beta, G0 )
        # plt.plot( beta, data_raw[:, i] - G0 )
        # plt.show( )
        # pos                     = np.argmax( np.abs( data_raw[:, i] - G0 ) )
        # # G1                      = G0 * data_raw[pos, i] / G0[pos]
        # alpha                   = 1.00001
        # f                       = np.exp( ( (-1.0+(alpha-mu[i])) / ( alpha - mu[i])**8.0 ) )
        # # G1                      = G0 + G0 * f * ( data[pos, i] / G0[pos] - 1.0 )
        # G1                      = G0 * f * ( data[pos, i] / G0[pos] )

        # # G0                      = G0_Magee_1( beta, mu[i] )
        # # G0[0]                   = 0.0
        # # pos                     = np.argmax( np.abs( data[:, i] - G0 ) )
        # # alpha                   = 1.00001
        # # f                       = np.exp( ( (-1.0+(alpha-mu[i])) / ( alpha - mu[i])**4.0 ) )
        # # G1                      = G0 + G0 * f * ( data[pos, i] / G0[pos] - 1.0 )
        # # F                       = 1.0 - np.exp( -c * ( ( beta / b )**2.0 ) * np.exp( - c * ( mu[i] * d ) **2.0 ) )
        # F                       = 1.0 - np.exp( -c * ( ( beta / b )**10) )

        # data[:, i]              = data_raw[:, i] - G1 * F
        data[:, i]              = data_raw[:, i] - G0

    # data = data / np.max( np.abs( data ) )

    # Interpolate database
    interp  = sp.interpolate.RegularGridInterpolator( ( beta, mu ), data, method="cubic" )

    # Interpolate over reduced domain
    x_min   = 6.688
    x_max   = 6.719
    y_min   = 0.154
    y_max   = 0.1778
    N       = 300
    xi      = np.linspace( x_min, x_max, N )
    yi      = np.linspace( y_min, y_max, N )
    XI, YI  = np.meshgrid( xi, yi, indexing="ij" )
    FI      = interp( ( XI, YI ) )
    interp1 = sp.interpolate.RegularGridInterpolator( ( xi, yi ), FI, method="cubic", bounds_error=False )

    N2      = 4
    xi2     = np.linspace( x_min, x_max, N2 )
    yi2     = np.linspace( y_min, y_max, N2 )
    XI2,YI2 = np.meshgrid( xi2, yi2, indexing="ij" )
    FI2     = interp( ( XI2, YI2 ) )
    interp2 = sp.interpolate.RegularGridInterpolator( ( xi2, yi2 ), FI2, method="cubic", bounds_error=True )

    N3      = 5
    xi3     = np.linspace( x_min, x_max, N3 )
    yi3     = np.linspace( y_min, y_max, N3 )
    XI3,YI3 = np.meshgrid( xi3, yi3, indexing="ij" )

    FT1     = interp1( ( XI3, YI3 ) )
    FT2     = interp2( ( XI3, YI3 ) )

    print( FT1 )
    print( FT2 )

    # Plot graphs
    x = np.array( [ x_min, x_max, x_max, x_min ] )
    y = np.array( [ y_min, y_min, y_max, y_max ] )

    fig = plt.figure( )
    ax  = fig.add_subplot( 111 )
    cnf = ax.contourf( B, M, np.log10( np.abs( data ) ), cmap="jet", levels=50 )
    for i in range( x.shape[0] ):
        j = (i+1)%x.shape[0]
        ax.plot( [ x[i], x[j] ], [ y[i], y[j] ], color="red", linestyle="--" )
    
    ax.set_yscale( "log" )
    plt.colorbar( cnf, ax=ax )

    fig = plt.figure( )
    ax  = fig.add_subplot( 111, projection="3d" )
    cnf = ax.plot_surface( B, np.log10( M ), data, cmap=plt.cm.jet, antialiased=False )
    
    plt.colorbar( cnf, ax=ax )

    fig = plt.figure( )
    ax  = fig.add_subplot( 111, projection="3d" )
    ax.plot_surface( XI, YI, FI, cmap=cm.jet )

    fig = plt.figure( )
    ax  = fig.add_subplot( 111 )
    cnf = ax.contourf( XI3, YI3, np.abs( FT1-FT2 ) )
    ax.set_yscale( "log" )
    plt.colorbar( cnf, ax=ax )

    plt.show( )


def updown_cross( t: np.ndarray, f: np.ndarray, is_show=False ) -> np.ndarray:
    # Remove mean value
    f = f - f.mean( )

    # Cross along Y zero axis
    pos_up      = ( f[1:] > 0.0 ) & ( f[:-1] < 0.0 )
    pos_down    = ( f[1:] < 0.0 ) & ( f[:-1] > 0.0 )
    pos         = np.where( pos_up | pos_down )[0]

    # Loop over cross positions to get the periods
    pos_max     = np.zeros( ( pos.shape[0]-1, ), dtype=int )

    for i in range( pos.shape[0] - 1 ):
        # Define interval bounds
        pi          = np.argmax( np.abs( f[ pos[i]:pos[i+1] ] ) )
        pos_max[i]  = pos[i] + pi

    # Get semi-periods
    if is_show:
        plt.plot( t, f )
        plt.plot( t[pos_max], f[pos_max], "o" )
        plt.show( )

    return pos_max


if __name__ == "__main__":
    # t = np.linspace( 0, 100, 10000 )
    # dt = t[1] - t[0]
    # g = G_Magee( t, 0.0 )
    # g0 = G0_Magee( t, 0.0 )

    # plt.plot( t, g, label="g" )
    # plt.plot( t, g0, label="g0" )
    # plt.plot( ( t[:-1] + t[1:] ) / 2.0, np.diff( g )/dt, label="g0*" )
    # plt.legend( )
    # plt.show( )
    # main( )
    # plot_database( )
    folder_path = r"D:\sergio\developments\SeaMotions\aux_tools\0_databases\2_infinite_water_depth_time"
    # generate_dGdt( folder_path )
    # generate_dGdtx( folder_path )
    # generate_dGdtxx( folder_path )
    # generate_dGdtt( folder_path )
    # generate_dGdttx( folder_path )
    generate_dGdttxx( folder_path )