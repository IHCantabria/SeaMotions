
# Import general usage libraries
import copy
import h5py
import os

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
    """
    A class for performing Chebyshev polynomial interpolation on 1D data.
    
    This class fits data using Chebyshev polynomials by solving for coefficients
    at Chebyshev nodes. It provides a convenient interface for evaluating the
    fitted polynomial and optionally displaying the fit quality.
    
    Attributes
    ----------
    _coeffs : np.ndarray
        Array of Chebyshev polynomial coefficients.
    _vandermat : np.ndarray
        Vandermonde matrix used for solving the interpolation system.
    _rhs : np.ndarray
        Right-hand side values (function values at Chebyshev nodes).
    _show_fit : bool
        Flag indicating whether to display a plot showing the fit quality.
    
    Parameters
    ----------
    rhs : np.ndarray
        1D array of function values to be interpolated. The number of values
        determines the order of the Chebyshev polynomial (order = len(rhs) - 1).
    show_fit : bool, optional
        If True, displays diagnostic plots showing the original data, fitted
        polynomial, differences, and coefficient magnitudes. Default is False.
    
    Examples
    --------
    >>> rhs = np.array([1.0, 2.0, 1.5, 0.5])
    >>> cheby = ChebyFit(rhs, show_fit=True)
    >>> values = cheby(np.array([-1, 0, 1]))  # Evaluate at specific points
    """

    def __init__( self, rhs: np.ndarray, show_fit=False ) -> None:
        # Define class properties
        self._coeffs    = np.ndarray
        self._vandermat = np.ndarray
        self._rhs       = rhs
        self._show_fit  = show_fit

        # Initialize class by calculating the interpolation coefficients
        self._initialize( )

    def __call__( self, x: np.ndarray ) -> np.ndarray:
        """
        Evaluate the Chebyshev polynomial at given x values.
        
        Parameters
        ----------
        x : np.ndarray
            Array of x-coordinates at which to evaluate the polynomial.
            Values should typically be in the range [-1, 1] for standard
            Chebyshev interpolation.
        
        Returns
        -------
        np.ndarray
            Evaluated polynomial values at the given x-coordinates.
        """
        return self.evaluate( x )

    def _initialize( self ) -> None:
        """
        Initialize the Chebyshev interpolation by computing coefficients.
        
        This method computes Chebyshev nodes, constructs the Vandermonde matrix,
        solves for the Chebyshev coefficients, and optionally displays diagnostic
        plots showing the fit quality and coefficient decay.
        """
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
        """
        Evaluate the Chebyshev polynomial at given x values.
        
        Parameters
        ----------
        x : np.ndarray
            Array of x-coordinates at which to evaluate the polynomial.
        
        Returns
        -------
        np.ndarray
            Evaluated polynomial values at the given x-coordinates.
        """
        return np.polynomial.chebyshev.chebval( x, self._coeffs )
    

class ChebyFitHP:
    """
    A class for HP refinement using 1D Chebyshev polynomial fitting with adaptive refinement.
    
    This class takes a spline fit object and converts it to an adaptive Chebyshev
    polynomial representation using hierarchical refinement. It provides efficient
    evaluation through a hash table-based block structure.
    
    Attributes
    ----------
    _fit_fcn : RefLevel
        The adaptive refinement level structure containing Chebyshev coefficients.
    _spf : SplineFit
        The input spline fit object to be converted.
    cum_coeffs : np.ndarray
        Cumulative array of all Chebyshev coefficients and their orders.
    cheby_coeffs : np.ndarray
        Array of Chebyshev polynomial coefficients.
    ncx : np.ndarray
        Array of Chebyshev polynomial orders for each coefficient.
    max_cheby_order : int
        Maximum Chebyshev polynomial order used in the fit.
    max_ref_level : int
        Maximum refinement level achieved.
    x_min : float
        Minimum x value of the interpolation domain.
    x_max : float
        Maximum x value of the interpolation domain.
    x_log_scale : bool
        Whether x-axis uses logarithmic scaling.
    fcn_log_scale : bool
        Whether function values use logarithmic scaling.
    intervals_np : int
        Number of intervals in the hash table (2^max_ref_level).
    dx : float
        Width of each interval in the hash table.
    x : np.ndarray
        Array of interval boundaries.
    xm : np.ndarray
        Array of interval midpoints.
    blocks_np : int
        Total number of blocks in the hash table.
    blocks_start : np.ndarray
        Starting index in coefficient array for each block.
    blocks_coeffs_np : np.ndarray
        Number of coefficients for each block.
    blocks_max_cheby_order : np.ndarray
        Maximum Chebyshev order for each block.
    x_min_vec : np.ndarray
        Minimum x value for each block.
    x_max_vec : np.ndarray
        Maximum x value for each block.
    
    Parameters
    ----------
    sp_fit : SplineFit
        A SplineFit object to be converted to adaptive Chebyshev representation.
    
    Examples
    --------
    >>> x = np.linspace(0, 10, 100)
    >>> y = np.sin(x)
    >>> spline = SplineFit(x, y)
    >>> cheby_hp = ChebyFitHP(spline)
    >>> value = cheby_hp(5.0)  # Evaluate at x=5.0
    
    Notes
    -----
    The adaptive refinement process automatically subdivides regions where
    the fit error exceeds specified tolerances, providing high accuracy with
    minimal coefficient storage.
    """

    def __init__( self, sp_fit ) -> None:
        # Define class attribues
        self._fit_fcn   = RefLevel
        self._spf       = sp_fit

        # Start fit
        self._initialize( )

    def __call__( self, xe: float ) -> float:
        """
        Evaluate the Chebyshev fit at a given x value.
        
        Parameters
        ----------
        xe : float
            The x-coordinate at which to evaluate the fit.
        
        Returns
        -------
        float
            The interpolated value at xe.
        """
        return self.evaluate( xe )

    def _fit( self ) -> RefLevel:
        """
        Perform adaptive Chebyshev fitting on the spline function.
        
        Creates a FitProperties object with default settings for 1D fitting,
        performs initial fit, and applies adaptive refinement to achieve
        specified tolerances.
        
        Returns
        -------
        RefLevel
            The root refinement level containing the adaptive fit structure.
        """
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
        """
        Initialize the adaptive Chebyshev fit.
        
        Performs the fitting process and recomposes the result into an
        efficient hash table structure for fast evaluation.
        """
        # Fit to target function by using adaptive interface
        self._fit_fcn = self._fit( )

        # Recompose the adaptive fit to have a working interpolation scheme
        self._recompose_fit( self._fit_fcn )

    def _recompose_fit( self, ref_level: RefLevel ) -> None:
        """
        Recompose the adaptive fit into a hash table structure.
        
        Extracts coefficient data from the refinement level structure and
        organizes it into arrays for efficient block-based evaluation.
        
        Parameters
        ----------
        ref_level : RefLevel
            The root refinement level containing the fit data.
        """
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
        """
        Evaluate the Chebyshev fit at a given x value.
        
        Uses the hash table to quickly locate the appropriate block and
        evaluates the corresponding Chebyshev polynomial.
        
        Parameters
        ----------
        xe : float
            The x-coordinate at which to evaluate the fit.
        
        Returns
        -------
        float
            The interpolated value at xe.
        """
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
        """
        Get the underlying refinement level structure.
        
        Returns
        -------
        RefLevel
            The root RefLevel object containing the complete adaptive fit data.
        """
        return self._fit_fcn
    

class ChebyFitHP2D:
    """
    A class for high-precision 2D Chebyshev polynomial fitting using adaptive refinement.
    
    This class loads data from an HDF5 file and performs adaptive 2D Chebyshev polynomial
    fitting. It uses hierarchical refinement with custom subdivision strategies to achieve
    high accuracy. The resulting fit is organized into an efficient hash table structure
    for fast block-based evaluation and derivatives.
    
    Attributes
    ----------
    _fit_props : FitProperties
        Configuration object containing fitting parameters, tolerances, domain bounds,
        and scaling options.
    _data_fipath : str
        File path to the HDF5 database containing the raw data to be fitted.
    cum_coeffs : np.ndarray
        Cumulative array of all Chebyshev coefficients and their orders.
    cheby_coeffs : np.ndarray
        Array of Chebyshev polynomial coefficients.
    ncx : np.ndarray
        Array of x-direction Chebyshev polynomial orders for each coefficient.
    ncy : np.ndarray
        Array of y-direction Chebyshev polynomial orders for each coefficient.
    max_cheby_order : float
        Maximum Chebyshev polynomial order used across both directions.
    max_ref_level : int
        Maximum refinement level achieved during adaptive fitting.
    x_min : float
        Minimum x value of the interpolation domain.
    x_max : float
        Maximum x value of the interpolation domain.
    y_min : float
        Minimum y value of the interpolation domain.
    y_max : float
        Maximum y value of the interpolation domain.
    x_log_scale : bool
        Whether x-axis uses logarithmic scaling.
    y_log_scale : bool
        Whether y-axis uses logarithmic scaling.
    fcn_log_scale : bool
        Whether function values use logarithmic scaling.
    intervals_np : int
        Number of intervals in each direction (2^max_ref_level).
    dx : float
        Width of each interval in x-direction.
    dy : float
        Width of each interval in y-direction.
    x : np.ndarray
        Array of x-direction interval boundaries.
    xm : np.ndarray
        Array of x-direction interval midpoints.
    y : np.ndarray
        Array of y-direction interval boundaries.
    ym : np.ndarray
        Array of y-direction interval midpoints.
    blocks_np : int
        Total number of 2D blocks in the hash table.
    blocks_start : np.ndarray
        Starting index in coefficient array for each block.
    blocks_coeffs_np : np.ndarray
        Number of coefficients for each block.
    blocks_max_cheby_order : np.ndarray
        Maximum Chebyshev order for each block.
    x_min_vec : np.ndarray
        Minimum x value for each block.
    x_max_vec : np.ndarray
        Maximum x value for each block.
    y_min_vec : np.ndarray
        Minimum y value for each block.
    y_max_vec : np.ndarray
        Maximum y value for each block.
    
    Parameters
    ----------
    data_fipath : str
        Path to HDF5 file containing 'mu', 'beta', and 'fcn' datasets.
    fit_props : FitProperties
        Configuration object defining fitting behavior, domain, tolerances,
        and refinement parameters.
    
    Examples
    --------
    >>> fit_props = FitProperties()
    >>> fit_props.x_min = 0.0
    >>> fit_props.x_max = 30.0
    >>> fit_props.y_min = -4.0
    >>> fit_props.y_max = 0.0
    >>> cheby_2d = ChebyFitHP2D("data.h5", fit_props)
    >>> value = cheby_2d(5.0, -2.0)  # Evaluate at (x=5.0, y=-2.0)
    >>> dx_value = cheby_2d.evaluate_dx(5.0, -2.0)  # x-derivative
    >>> dy_value = cheby_2d.evaluate_dy(5.0, -2.0)  # y-derivative
    
    Notes
    -----
    The adaptive refinement process uses a hierarchical subdivision strategy
    with predefined y-levels to efficiently capture function behavior across
    different regions. The fitting uses Chebyshev polynomials evaluated at
    Chebyshev nodes to minimize interpolation error.
    
    The class provides efficient evaluation through a 2D hash table that quickly
    locates the appropriate polynomial block for any query point. Derivatives
    are computed analytically using Chebyshev derivative formulas.
    """

    def __init__( self, data_fipath: str, fit_props: FitProperties ) -> None:
        # Define class attribues
        self._fit_props     = fit_props
        self._data_fipath   = data_fipath

        # Start fit
        self._initialize( )

    def __call__( self, xe: float, ye: float ) -> float:
        """
        Evaluate the 2D Chebyshev fit at given (x, y) coordinates.
        
        Parameters
        ----------
        xe : float
            The x-coordinate at which to evaluate the fit.
        ye : float
            The y-coordinate at which to evaluate the fit.
        
        Returns
        -------
        float
            The interpolated value at (xe, ye).
        """
        return self.evaluate( xe, ye )

    def _fit( self ) -> RefLevel:
        """
        Perform adaptive 2D Chebyshev fitting on the loaded data.
        
        Loads data from HDF5 file, creates an interpolator, and performs
        hierarchical adaptive refinement using predefined subdivision levels
        in the y-direction and uniform subdivision in the x-direction.
        
        Returns
        -------
        RefLevel
            The root refinement level containing the complete adaptive fit structure.
        """
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
        """
        Get Chebyshev coefficients for the block containing point (xe, ye).
        
        Locates the appropriate block in the hash table and assembles the
        coefficient matrix for that block.
        
        Parameters
        ----------
        xe : float
            The x-coordinate of the query point.
        ye : float
            The y-coordinate of the query point.
        
        Returns
        -------
        list
            A list containing [coeffs_i, npos] where coeffs_i is a 2D array
            of Chebyshev coefficients and npos is the block index.
        """
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
        """
        Initialize the adaptive 2D Chebyshev fit.
        
        Performs the fitting process and recomposes the result into an
        efficient 2D hash table structure for fast evaluation.
        """
        # Fit to target function by using adaptive interface
        ref_level = self._fit( )

        # Recompose the adaptive fit to have a working interpolation scheme
        self._recompose_fit( ref_level )

    def _recompose_fit( self, ref_level: RefLevel ) -> None:
        """
        Recompose the adaptive fit into a 2D hash table structure.
        
        Extracts coefficient data from the refinement level structure and
        organizes it into arrays for efficient block-based evaluation in 2D.
        Creates a regular grid of blocks based on the maximum refinement level.
        
        Parameters
        ----------
        ref_level : RefLevel
            The root refinement level containing the fit data.
        """
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
        """
        Evaluate the 2D Chebyshev fit at given coordinates.
        
        Uses the hash table to locate the appropriate block and evaluates
        the corresponding 2D Chebyshev polynomial.
        
        Parameters
        ----------
        xe : float
            The x-coordinate at which to evaluate the fit.
        ye : float
            The y-coordinate at which to evaluate the fit.
        
        Returns
        -------
        float
            The interpolated value at (xe, ye).
        """
        # Get patch properties
        coeffs_i, npos = self._get_coeffs( xe, ye )

        # Convert Xe to local coordinates [-1, 1]
        xeu = map_to_unit( self.x_min_vec[ npos ], self.x_max_vec[ npos ], xe )
        yeu = map_to_unit( self.y_min_vec[ npos ], self.y_max_vec[ npos ], ye )

        # Calculate polynomial value
        val = np.polynomial.chebyshev.chebval2d( xeu, yeu, coeffs_i ).sum( )

        return val
    
    def evaluate_dx( self, xe: float, ye: float ) -> float:
        """
        Evaluate the x-derivative of the 2D Chebyshev fit.
        
        Computes the partial derivative with respect to x by differentiating
        the Chebyshev polynomial coefficients analytically.
        
        Parameters
        ----------
        xe : float
            The x-coordinate at which to evaluate the derivative.
        ye : float
            The y-coordinate at which to evaluate the derivative.
        
        Returns
        -------
        float
            The x-derivative value at (xe, ye).
        """
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
        """
        Evaluate the y-derivative of the 2D Chebyshev fit.
        
        Computes the partial derivative with respect to y by differentiating
        the Chebyshev polynomial coefficients analytically.
        
        Parameters
        ----------
        xe : float
            The x-coordinate at which to evaluate the derivative.
        ye : float
            The y-coordinate at which to evaluate the derivative.
        
        Returns
        -------
        float
            The y-derivative value at (xe, ye).
        """
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
    """
    A class to wrap up Chebyshev 1D interpolators and to have a convinient 
    way to use them through a easy interface. The idea is to receive a G0CoeffsFit
    instance and to reinterpolate the function using 1D Chebyshev polynomials in this 
    way the errors associated to the interpolation process will be the same as 
    the ones using on the high performance c++ routines.
    
    Attributes
    ----------
    _G0fit : G0CoeffsFit
        The input G0CoeffsFit containing the original interpolators by using high resolution spline.
    _interp_fcn : list
        List of ChebyFitHP objects, one for each converted interpolator.
    
    Parameters
    ----------
    G0_fit : G0CoeffsFit
        A G0 fit object that provides interpolators via get_interpolators() method.
    
    Examples
    --------
    >>> G0_fit = G0CoeffsFit( mu, coeffs )  # Some G0 fit object
    >>> cheby_alpha = ChebyFitHPAlpha( G0_fit )
    >>> alpha_values = cheby_alpha( 0.5 )  # Evaluate at x=0.5
    """

    def __init__( self, G0_fit ) -> None:
        # Define class attributes
        self._G0fit     = G0_fit
        self._interp_fcn    = [ ]

        # Initialize class
        self._initialize( )

    def __call__( self, x: float ) -> np.ndarray:
        """
        Evaluate all Chebyshev interpolators at a given x value.
        
        Parameters
        ----------
        x : float
            The x-coordinate at which to evaluate all interpolators.
        
        Returns
        -------
        np.ndarray
            Array of evaluated values from all interpolators.
        """
        return self.evaluate( x )
    
    def _initialize( self ) -> None:
        """
        Initialize the Chebyshev interpolators from the G0 fit.
        
        Converts each interpolator from the G0 fit object into a ChebyFitHP
        object for high-precision Chebyshev polynomial-based interpolation.
        """
        # Loop over G0 interpolators to regenerate with Chebyshev polynomials
        for mi in self._G0fit.get_interpolators( ):
            self._interp_fcn.append( ChebyFitHP( mi ) )
    
    def evaluate( self, x: float ) -> np.ndarray:
        """
        Evaluate all Chebyshev interpolators at a given x value.
        
        Parameters
        ----------
        x : float
            The x-coordinate at which to evaluate all interpolators.
        
        Returns
        -------
        np.ndarray
            1D array containing the evaluated values from all Chebyshev
            interpolators at the specified x-coordinate.
        """
        result = np.zeros( ( len( self._interp_fcn, ) ) )
        for i in range( result.shape[0] ):
            result[i] = self._interp_fcn[i]( x )

        return result
    
    def get_interpolators( self ) -> list:
        """
        Get the list of Chebyshev interpolation functions.
        
        Returns
        -------
        list
            List of ChebyFitHP objects used for interpolation.
        """
        return self._interp_fcn


class SplineFit:
    """
    A class for performing cubic spline interpolation on 1D data. 
    Just a convinient wrapper to scipy.interpolate.CubicSpline to 
    manage data and to show fit quality.

    
    Attributes
    ----------
    _interpf : scipy.interpolate.CubicSpline
        The cubic spline interpolation object.
    _x : np.ndarray
        1D array of x-coordinates of the data points.
    _y : np.ndarray
        1D array of y-coordinates of the data points.
    _show_fit : bool
        Flag indicating whether to display a plot comparing original and 
        interpolated data.
    
    Parameters
    ----------
    x : np.ndarray
        1D array of x-coordinates where data is defined.
    y : np.ndarray
        1D array of y-values corresponding to the x-coordinates.
    show_fit : bool, optional
        If True, displays a plot showing the original data points and the
        interpolated curve. Default is False.
    
    Examples
    --------
    >>> x = np.array([0, 1, 2, 3, 4])
    >>> y = np.array([0, 1, 4, 9, 16])
    >>> spline = SplineFit(x, y, show_fit=True)
    >>> interpolated_value = spline(2.5)
    """

    def __init__( self, x: np.ndarray, y:np.ndarray, show_fit=False ) -> None:
        # Define class attributes
        self._interpf   = None
        self._x         = x
        self._y         = y
        self._show_fit  = show_fit

        # Initialize class
        self._initialize( )

    def __call__( self, x: np.ndarray ) -> np.ndarray:
        """
        Evaluate the spline interpolation at given x values.
        
        Parameters
        ----------
        x : np.ndarray
            Array of x-coordinates at which to evaluate the interpolation.
        
        Returns
        -------
        np.ndarray
            Interpolated y-values at the given x-coordinates.
        """
        return self.evaluate( x )

    def _initialize( self ) -> None:
        """
        Initialize the cubic spline interpolator.
        
        Creates the scipy CubicSpline object from the input data and optionally
        displays a comparison plot if show_fit is True.
        """
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
        """
        Evaluate the spline interpolation at given x values.
        
        Parameters
        ----------
        x : np.ndarray
            Array of x-coordinates at which to evaluate the interpolation.
        
        Returns
        -------
        np.ndarray
            Interpolated y-values at the given x-coordinates.
        """
        return self._interpf( x )
    

class G0CoeffsFit:
    """
    A class for fitting and interpolating G0 fit coefficients using cubic splines.
    
    This class takes an array of coefficients and fits each row independently
    using cubic spline interpolation as a function of the mu parameter. Each row 
    of coefficients corresponds to a different phase angles used in the fit.
    
    Attributes
    ----------
    _coeffs : np.ndarray
        Array of coefficients to be fitted, shape (n_shifts, n_mu_points).
    _interp_fcn : list
        List of SplineFit objects, one for each coefficient row.
    _mu : np.ndarray
        Array of mu values corresponding to the coefficient columns.
    
    Parameters
    ----------
    mu : np.ndarray
        1D array of mu values at which coefficients are defined.
    coeffs : np.ndarray
        2D array of coefficients, where each row represents a different
        set of coefficient for phase shift to be fitted as a function of mu.
    
    Examples
    --------
    >>> mu = np.linspace(0, 1, 10)
    >>> coeffs = np.random.rand(5, 10)
    >>> fit = G0CoeffsFit(mu, coeffs)
    >>> result = fit(0.5)  # Evaluate all coefficients at mu=0.5
    """

    def __init__( self, mu: np.ndarray, coeffs: np.ndarray ) -> None:
        # Define class attributes
        self._coeffs        = coeffs
        self._interp_fcn    = [ ]
        self._mu            = mu

        # Initialize class by fitting each of the coefficients using chebyshev expansion
        self._initialize( )

    def __call__( self, x: float ) -> np.ndarray:
        """
        Evaluate all fitted coefficients at a given mu value.
        
        Parameters
        ----------
        x : float
            The mu value at which to evaluate the coefficients.
        
        Returns
        -------
        np.ndarray
            Array of evaluated coefficient values at the given mu.
        """
        return self.evaluate( x )

    def _initialize( self ) -> None:
        """
        Initialize the interpolation by creating a SplineFit for each coefficient row.
        
        This method creates cubic spline interpolators for each row of the
        coefficient matrix.
        """
        for i in range( self._coeffs.shape[0] ):
            self._interp_fcn.append( SplineFit( self._mu, self._coeffs[i, :], show_fit=True ) )

    def evaluate( self, x: float ) -> np.ndarray:
        """
        Evaluate all fitted coefficients at a given mu value.
        
        Parameters
        ----------
        x : float
            The mu value at which to evaluate the coefficients.
        
        Returns
        -------
        np.ndarray
            1D array containing the evaluated values of all coefficients
            at the specified mu value.
        """
        results = np.zeros( ( self._coeffs.shape[0], ) )
        for i in range( self._coeffs.shape[0] ):
            results[i] = self._interp_fcn[i]( x )

        return results
    
    def get_interpolators( self ) -> list:
        """
        Get the list of interpolation functions.
        
        Returns
        -------
        list
            List of SplineFit objects used for interpolating each coefficient.
        """
        return self._interp_fcn


def calculate_residual_G0( fit_props: FitProperties, target_fcn, interp_raw, data_raw, beta: np.ndarray, mu: np.ndarray ) -> None:
    """
    Calculate the residual between raw data and G0.
    
    This function fits G0 coefficients at specified mu values, creates interpolators
    for the coefficients, and subtracts the G0 approximation from the raw data to
    obtain a residual function.
    
    Parameters
    ----------
    fit_props : FitProperties
        Configuration object containing G0 fitting parameters including:
        - G0_cheby_np: Number of Chebyshev interpolation points
        - G0_mui: Array of mu values for G0 fitting
        - G0_mui_log: Log10 of mu values
        - G0_coeffs: Array to store fitted coefficients
        - G0_cheby_intv: Interval for Chebyshev interpolation in log10(mu)
        - alpha_shift: Array of phase shift angles for fitting
        - y_min, y_max: Domain bounds in y-direction
    target_fcn : callable
        Target function for G0 evaluation with signature fcn(beta, mu, alpha).
    interp_raw : scipy.interpolate.RegularGridInterpolator
        Interpolator for the raw data as function of (beta, mu).
    data_raw : np.ndarray
        2D array of raw data values, shape (n_beta, n_mu).
    beta : np.ndarray
        1D array of beta values (spatial coordinate).
    mu : np.ndarray
        1D array of mu values (typically cos(theta) or similar parameter).
    
    Returns
    -------
    data : np.ndarray
        2D array of residual data after subtracting G0 approximation,
        shape (n_beta, n_mu).
    G0_chfhpa : ChebyFitHPAlpha
        High-precision Chebyshev fit object for G0 coefficient interpolation
        as a function of log10(mu).
    
    Notes
    -----
    The function performs the following steps:
    1. Fits G0 coefficients at discrete mu points using least-squares
    2. Creates spline interpolators for coefficients as function of mu
    3. Converts spline interpolators to high-precision Chebyshev form
    4. Evaluates G0 approximation at all data points
    5. Subtracts G0 from raw data to obtain residual
    6. Computes normalization factor based on maximum residual magnitude
    
    The G0 approximation is only applied for mu values within the specified
    Chebyshev interpolation interval (G0_cheby_intv).
    
    Examples
    --------
    >>> fit_props = FitProperties()
    >>> define_G0_fit_properties(fit_props, 10, mu)
    >>> data_resid, g0_fit = calculate_residual_G0(
    ...     fit_props, dGdt_G0, interp_raw, data_raw, beta, mu
    ... )
    """
    # Loop over G0 interpolation points to fit coefficients
    for i in range( fit_props.G0_cheby_np ):
        mui = np.ones_like( beta ) * fit_props.G0_mui[i]
        print( mui )
        _, fit_props.G0_coeffs[:, i] = fit_G0_fcn( target_fcn, beta, fit_props.G0_mui[i], interp_raw( ( beta, mui ) ), fit_props.alpha_shift, show_recon=False )

    if fit_props.G0_cheby_np > 1:
        mu_fit        = fit_props.G0_mui_log
        G0_coeffs_fit = fit_props.G0_coeffs

    else:
        np_fit        = 4
        mu_fit        = np.linspace( fit_props.y_min, fit_props.y_max, np_fit )
        G0_coeffs_fit = fit_props.G0_coeffs * np.ones( ( np_fit, ) )

    G0_fit     = G0CoeffsFit( mu_fit, G0_coeffs_fit )
    G0_chfhpa  = ChebyFitHPAlpha( G0_fit )

    data = data_raw.copy( )
    G0   = np.zeros( ( data.shape[0], ) )
    for i in range( mu.shape[0] ):
        print( "I: ", i )
        if mu[i] < 10**fit_props.G0_cheby_intv[1]:
            print( "Re-fitting..." )
            coeffs_i    = G0_chfhpa( np.log10( mu[i] ) )
        
        G0          = eval_G0_fcn( target_fcn, coeffs_i, beta, mu[i], fit_props.alpha_shift )
        data[:, i]  = data_raw[:, i] - G0

    fit_props.G0_norm_f  = np.max( np.abs( data ) ) * 10
    print( "norm_f: ", fit_props.G0_norm_f )

    return data, G0_chfhpa


def define_G0_fit_properties( fit_props: FitProperties, n: int, mu: np.ndarray, points_dist="linear" ) -> None:
    """
    Define G0 fitting properties.
    
    This function configures the FitProperties object with parameters needed for
    G0 coefficient fitting, including the interpolation interval, number of points,
    and distribution of fitting points in the mu domain.
    
    Parameters
    ----------
    fit_props : FitProperties
        Configuration object to be updated with G0 fitting parameters.
        The function modifies the following attributes:
        - G0_cheby_intv: Log10 interval bounds [log10(mu_min), log10(mu_max)]
        - G0_cheby_np: Number of interpolation points
        - G0_mui_log: Log10 of mu values at fitting points
        - G0_mui: Actual mu values at fitting points (10^G0_mui_log)
        - G0_coeffs: Initialized coefficient array (zeros)
    n : int
        Number of interpolation points to use for G0 fitting. If np < 2,
        a single point fit is performed.
    mu : np.ndarray
        1D array of mu values defining the domain bounds. The first and last
        values are used to set the interpolation interval.
    points_dist : str, optional
        Distribution type for fitting points. Options:
        - "linear": Uniformly spaced points in [-1, 1] interval
        - "chebyshev": Chebyshev nodes (roots of Chebyshev polynomial)
        Default is "linear".
    
    Returns
    -------
    None
        The function modifies fit_props in-place.
    
    Raises
    ------
    ValueError
        If points_dist is not "linear" or "chebyshev".
    
    Notes
    -----
    The function uses logarithmic spacing in mu to better capture behavior
    across multiple orders of magnitude. The Chebyshev node distribution
    provides better interpolation accuracy near the boundaries of the interval.
    
    For np < 2, the function sets up a minimal fit with a single point,
    which is useful for constant G0 approximations.
    
    Examples
    --------
    >>> fit_props = FitProperties()
    >>> mu = np.logspace(-4, 0, 100)
    >>> define_G0_fit_properties(fit_props, 10, mu, points_dist="chebyshev")
    >>> print(fit_props.G0_cheby_np)
    10
    >>> print(fit_props.G0_cheby_intv)
    [-4.0, 0.0]
    """
    # Define G0 fit points
    if points_dist == "linear":
        cheby_roots = np.linspace( -1.0, 1.0, n )
    elif points_dist == "chebyshev":
        cheby_roots = sp.special.roots_chebyt( n )[0]
    else:
        raise ValueError( f"Unknown points distribution: {points_dist}" )
    
    # Define G0 Chebyshev fit properties
    fit_props.G0_cheby_intv     = np.log10( np.array( [ mu[0], mu[-1] ] ) )
    fit_props.G0_cheby_np       = n
    fit_props.G0_mui_log        = np.array( [ fit_props.G0_cheby_intv[0] ] ) if fit_props.G0_cheby_np < 2 else map_to_interval( fit_props.G0_cheby_intv[0], fit_props.G0_cheby_intv[1], cheby_roots )
    fit_props.G0_mui            = 10**fit_props.G0_mui_log
    fit_props.G0_coeffs         = np.zeros( ( fit_props.alpha_shift.shape[0], fit_props.G0_cheby_np ) )


def eval_G0_fcn( fcn_fit, coeffs: np.ndarray, beta: np.ndarray, mu: float, alpha: np.ndarray, show_recon=False ) -> np.ndarray:
    """
    Evaluate the G0 function using fitted coefficients.
    
    This function reconstructs the G0 approximation by forming a linear combination
    of basis functions evaluated at different phase shift angles (alpha values).
    The reconstruction is: G0(beta, mu) = sum_i coeffs[i] * fcn_fit(beta, mu, alpha[i])
    
    Parameters
    ----------
    fcn_fit : callable
        Target function for G0 evaluation with signature fcn_fit(beta, mu, alpha).
        This is typically one of the G0 basis functions like dGdt_G0, dGdtx_G0, etc.
    coeffs : np.ndarray
        1D array of fitted coefficients for the linear combination, shape (n_alpha,).
    beta : np.ndarray
        1D array of beta values (temporal coordinate) at which to evaluate.
    mu : float
        Single mu value.
    alpha : np.ndarray
        1D array of phase shift angles used in the basis functions, shape (n_alpha,).
    show_recon : bool, optional
        If True, displays diagnostic plots showing the system matrix, coefficients,
        and reconstructed function. Default is False.
    
    Returns
    -------
    recon : np.ndarray
        1D array of reconstructed G0 values at the given beta values, shape (n_beta,).
    
    Notes
    -----
    The function builds a system matrix A where A[i,j] = fcn_fit(beta[i], mu, alpha[j]),
    then computes the reconstruction as recon = A @ coeffs.
    
    The visualization (when show_recon=True) displays:
    - Log10 of absolute values of the system matrix
    - Log10 of absolute coefficient values
    - Log10 of absolute reconstructed values
    
    Examples
    --------
    >>> beta = np.linspace(0, 10, 100)
    >>> mu = 0.5
    >>> alpha = np.array([0.0, 0.01, 0.1])
    >>> coeffs = np.array([1.0, 0.5, 0.2])
    >>> G0 = eval_G0_fcn(dGdt_G0, coeffs, beta, mu, alpha, show_recon=True)
    """

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


def fit_G0_fcn( fcn_fit, beta: np.ndarray, mu: float, f: np.ndarray, alpha: np.ndarray, show_recon=False ) -> None:
    """
    Fit G0 coefficients using weighted least-squares regression.
    
    This function performs a weighted least-squares fit to determine the coefficients
    that best approximate a target function as a linear combination of basis functions
    evaluated at different phase shift angles. The weighting increases linearly with
    beta to emphasize fitting accuracy at larger values.
    
    Parameters
    ----------
    fcn_fit : callable
        Basis function for G0 evaluation with signature fcn_fit(beta, mu, alpha).
        This is typically one of the G0 basis functions like dGdt_G0, dGdtx_G0, etc.
    beta : np.ndarray
        1D array of beta values (temporal coordinate) at which to fit, shape (n_beta,).
    mu : float
        Single mu value.
    f : np.ndarray
        1D array of target function values to fit at the beta points, shape (n_beta,).
    alpha : np.ndarray
        1D array of phase shift angles used in the basis functions, shape (n_alpha,).
    show_recon : bool, optional
        If True, displays diagnostic plots showing:
        - Original data vs. reconstructed fit
        - Residual (difference between original and fit)
        - Coefficient magnitudes on log scale
        Default is False.
    
    Returns
    -------
    recon : np.ndarray
        1D array of reconstructed (fitted) values at the beta points, shape (n_beta,).
    coeffs : np.ndarray
        1D array of fitted coefficients for the linear combination, shape (n_alpha,).
    
    Notes
    -----
    The function uses a linear weighting scheme w = linspace(0, 1, n_beta) to give
    more importance to larger beta values. A small positive weight (1e-6) is assigned
    to the first point to avoid numerical issues with zero weight.
    
    The weighted least-squares problem is solved as:
        min ||W^(1/2) * A * coeffs - W^(1/2) * f||^2
    where W is a diagonal weight matrix and A is the system matrix.
    
    For regularized solutions, Tikhonov regularization can be enabled by uncommenting
    the relevant code section and adjusting the lambda parameter.
    
    Examples
    --------
    >>> beta = np.linspace(0, 10, 100)
    >>> mu = 0.5
    >>> alpha = np.array([0.0, 0.01, 0.1])
    >>> f_target = some_target_function(beta, mu)
    >>> recon, coeffs = fit_G0_fcn(dGdt_G0, beta, mu, f_target, alpha, show_recon=True)
    >>> print(f"Fitted coefficients: {coeffs}")
    """

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


def fit_residual_dGdt( folder_path: str, show_figs=False, is_square_ref=True, show_summary_fig=False ) -> None:
    # Define fit properties
    fit_props                   = FitProperties( )
    fit_props.dims              = 2
    fit_props.region_name       = "dGdt"
    fit_props.cheby_order_x     = 15
    fit_props.cheby_order_y     = 15
    fit_props.fcn_log_scale     = False
    fit_props.x_hpatch_np       = 30
    fit_props.x_log_scale       = False
    fit_props.x_max             = 30.0
    fit_props.x_min             = 0.0
    fit_props.y_hpatch_np       = 4
    fit_props.y_log_scale       = True
    fit_props.y_max             = np.log10( 0.9998 )
    fit_props.y_min             = -4.0
    fit_props.cheby_abs_tol     = 1E-4
    fit_props.cheby_abs_tol_f   = 2
    fit_props.cheby_rel_tol     = 1E-4
    fit_props.cheby_rel_tol_f   = 2
    fit_props.x_map_fcn         = lambda x: x
    fit_props.y_map_fcn         = lambda y: y
    fit_props.max_ref_level     = 10
    fit_props.alpha_shift       = np.array( [ 0.0 ] )

    fit_props.num_x             = fit_props.cheby_order_x + 1
    fit_props.num_x_fit         = fit_props.cheby_order_x + 1
    fit_props.num_y             = fit_props.cheby_order_y + 1
    fit_props.num_y_fit         = fit_props.cheby_order_y + 1

    fit_props.generate_fitting_matrix( )

    # Load database
    fipath = os.path.join( get_integrals_database_fopath( ), "1_time_domain", "Gt.h5" )
    with h5py.File( fipath, "r" ) as fid:
        mu          = fid[ "mu" ][:]
        beta        = fid[ "beta" ][:]
        data_raw    = fid[ "fcn" ][:]

    # Define G0 interpolation properties
    target_fcn      = dGdt_G0
    interp_raw      = sp.interpolate.RegularGridInterpolator( ( beta, mu ), data_raw )
    
    define_G0_fit_properties( fit_props, 1, mu, points_dist="chebyshev" )

    # Calculate residual function in between G and G0*
    data, chfhpa    = calculate_residual_G0( fit_props, target_fcn, interp_raw, data_raw, beta, mu )

    # Plot residual function
    plot_residual_function( beta, mu, data_raw, data )

    # Interpolate database
    fit_function_raw        = sp.interpolate.RegularGridInterpolator( ( beta, mu ), data )
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
    fit_props                   = FitProperties( )
    fit_props.dims              = 2
    fit_props.region_name       = "dGdtx"
    fit_props.cheby_order_x     = 30
    fit_props.cheby_order_y     = 15
    fit_props.fcn_log_scale     = False
    fit_props.x_log_scale       = False
    fit_props.x_hpatch_np       = 30
    fit_props.x_max             = 30.0
    fit_props.x_min             = 0.0
    fit_props.y_hpatch_np       = 4
    fit_props.y_log_scale       = True
    fit_props.y_max             = np.log10( 0.9998 )
    fit_props.y_min             = -4.0
    fit_props.cheby_abs_tol     = 1E-2
    fit_props.cheby_abs_tol_f   = 2
    fit_props.cheby_rel_tol     = 1E-4
    fit_props.cheby_rel_tol_f   = 2
    fit_props.x_map_fcn         = lambda x: x
    fit_props.y_map_fcn         = lambda y: y
    fit_props.max_ref_level     = 10
    fit_props.alpha_shift       = np.array( [ 1e-3, 1e-2, 1e-1 ] )

    fit_props.num_x             = fit_props.cheby_order_x + 1
    fit_props.num_x_fit         = fit_props.cheby_order_x + 1
    fit_props.num_y             = fit_props.cheby_order_y + 1
    fit_props.num_y_fit         = fit_props.cheby_order_y + 1

    fit_props.generate_fitting_matrix( )

    # Load database
    fipath = os.path.join( get_integrals_database_fopath( ), "1_time_domain", "Gtx.h5" )
    with h5py.File( fipath, "r" ) as fid:
        mu          = fid[ "mu" ][:]
        beta        = fid[ "beta" ][:]
        data_raw    = fid[ "fcn" ][:]

    # Define G0 interpolation properties
    target_fcn      = dGdtx_G0
    interp_raw      = sp.interpolate.RegularGridInterpolator( ( beta, mu ), data_raw )
    
    define_G0_fit_properties( fit_props, 1, mu, points_dist="chebyshev" )

    # Calculate residual function in between G and G0*
    data, chfhpa    = calculate_residual_G0( fit_props, target_fcn, interp_raw, data_raw, beta, mu )
    
    # Plot residual function
    plot_residual_function( beta, mu, data_raw, data )

    # Interpolate database
    fit_function_raw        = sp.interpolate.RegularGridInterpolator( ( beta, mu ), data )
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
    fit_props                   = FitProperties( )
    fit_props.dims              = 2
    fit_props.region_name       = "dGdtxx"
    fit_props.cheby_order_x     = 40
    fit_props.cheby_order_y     = 15
    fit_props.fcn_log_scale     = False
    fit_props.x_hpatch_np       = 30
    fit_props.x_log_scale       = False
    fit_props.x_max             = 30.0
    fit_props.x_min             = 0.0
    fit_props.y_hpatch_np       = 4
    fit_props.y_log_scale       = True
    fit_props.y_max             = np.log10( 0.9998 )
    fit_props.y_min             = -4.0
    fit_props.cheby_abs_tol     = 1E1
    fit_props.cheby_abs_tol_f   = 2
    fit_props.cheby_rel_tol     = 1E-4
    fit_props.cheby_rel_tol_f   = 2
    fit_props.x_map_fcn         = lambda x: x
    fit_props.y_map_fcn         = lambda y: y
    fit_props.max_ref_level     = 10
    fit_props.alpha_shift       = np.array( [ 1e-3, 1e-2, 1e-1 ] )

    fit_props.num_x             = fit_props.cheby_order_x + 1
    fit_props.num_x_fit         = fit_props.cheby_order_x + 1
    fit_props.num_y             = fit_props.cheby_order_y + 1
    fit_props.num_y_fit         = fit_props.cheby_order_y + 1

    fit_props.generate_fitting_matrix( )

    # Load database
    fipath = os.path.join( get_integrals_database_fopath( ), "1_time_domain", "Gtxx.h5" )
    with h5py.File( fipath, "r" ) as fid:
        mu          = fid[ "mu" ][:]
        beta        = fid[ "beta" ][:]
        data_raw    = fid[ "fcn" ][:]

    # Define G0 interpolation properties
    target_fcn      = dGdtxx_G0
    interp_raw      = sp.interpolate.RegularGridInterpolator( ( beta, mu ), data_raw )

    define_G0_fit_properties( fit_props, 1, mu, points_dist="chebyshev" )

    # Calculate residual function in between G and G0*
    data, chfhpa    = calculate_residual_G0( fit_props, target_fcn, interp_raw, data_raw, beta, mu )

    # Plot residual function
    plot_residual_function( beta, mu, data_raw, data )

    # Interpolate database
    fit_function_raw        = sp.interpolate.RegularGridInterpolator( ( beta, mu ), data )
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
    fit_props                   = FitProperties( )
    fit_props.dims              = 2
    fit_props.region_name       = "dGdtt"
    fit_props.cheby_order_x     = 40
    fit_props.cheby_order_y     = 15
    fit_props.fcn_log_scale     = False
    fit_props.x_hpatch_np       = 30
    fit_props.x_log_scale       = False
    fit_props.x_max             = 19.0
    fit_props.x_min             = 0.0
    fit_props.y_hpatch_np       = 4
    fit_props.y_log_scale       = True
    fit_props.y_max             = np.log10( 0.9998 )
    fit_props.y_min             = -4.0
    fit_props.cheby_abs_tol     = 1E-3
    fit_props.cheby_abs_tol_f   = 2
    fit_props.cheby_rel_tol     = 1E-4
    fit_props.cheby_rel_tol_f   = 2
    fit_props.x_map_fcn         = lambda x: x
    fit_props.y_map_fcn         = lambda y: y
    fit_props.max_ref_level     = 10
    fit_props.alpha_shift       = np.array( [ 1e-3, 1e-2, 1e-1 ] )

    fit_props.num_x             = fit_props.cheby_order_x + 1
    fit_props.num_x_fit         = fit_props.cheby_order_x + 1
    fit_props.num_y             = fit_props.cheby_order_y + 1
    fit_props.num_y_fit         = fit_props.cheby_order_y + 1

    fit_props.generate_fitting_matrix( )

    # Load database
    fipath = os.path.join( get_integrals_database_fopath( ), "1_time_domain", "Gtt.h5" )
    with h5py.File( fipath, "r" ) as fid:
        mu          = fid[ "mu" ][:]
        beta        = fid[ "beta" ][:]
        data_raw    = fid[ "fcn" ][:]

    # Define G0 interpolation properties
    target_fcn      = dGdtt_G0
    interp_raw      = sp.interpolate.RegularGridInterpolator( ( beta, mu ), data_raw )

    define_G0_fit_properties( fit_props, 1, mu, points_dist="linear" )

    # Calculate residual function in between G and G0*
    data, chfhpa    = calculate_residual_G0( fit_props, target_fcn, interp_raw, data_raw, beta, mu )

    # Plot residual function
    plot_residual_function( beta, mu, data_raw, data )

    # Interpolate database
    fit_function_raw        = sp.interpolate.RegularGridInterpolator( ( beta, mu ), data )
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
    fit_props                   = FitProperties( )
    fit_props.dims              = 2
    fit_props.region_name       = "dGdttx"
    fit_props.cheby_order_x     = 40
    fit_props.cheby_order_y     = 15
    fit_props.fcn_log_scale     = False
    fit_props.x_hpatch_np       = 30
    fit_props.x_log_scale       = False
    fit_props.x_max             = 30.0
    fit_props.x_min             = 0.0
    fit_props.y_hpatch_np       = 4
    fit_props.y_log_scale       = True
    fit_props.y_max             = np.log10( 0.999 )
    fit_props.y_min             = -4.0
    fit_props.cheby_abs_tol     = 1E0
    fit_props.cheby_abs_tol_f   = 2
    fit_props.cheby_rel_tol     = 1E-4
    fit_props.cheby_rel_tol_f   = 2
    fit_props.x_map_fcn         = lambda x: x
    fit_props.y_map_fcn         = lambda y: y
    fit_props.max_ref_level     = 10
    fit_props.alpha_shift       = np.array( [ 1e-3, 1e-2, 1e-1 ] )

    fit_props.num_x             = fit_props.cheby_order_x + 1
    fit_props.num_x_fit         = fit_props.cheby_order_x + 1
    fit_props.num_y             = fit_props.cheby_order_y + 1
    fit_props.num_y_fit         = fit_props.cheby_order_y + 1

    fit_props.generate_fitting_matrix( )

    # Load database
    fipath = os.path.join( get_integrals_database_fopath( ), "1_time_domain", "Gttx.h5" )
    with h5py.File( fipath, "r" ) as fid:
        mu          = fid[ "mu" ][:]
        beta        = fid[ "beta" ][:]
        data_raw    = fid[ "fcn" ][:]

    # Define G0 interpolation properties
    target_fcn      = dGdttx_G0
    interp_raw      = sp.interpolate.RegularGridInterpolator( ( beta, mu ), data_raw )

    define_G0_fit_properties( fit_props, 100, mu, points_dist="linear" )

    # Calculate residual function in between G and G0*
    data, chfhpa    = calculate_residual_G0( fit_props, target_fcn, interp_raw, data_raw, beta, mu )

    # Plot residual function
    plot_residual_function( beta, mu, data_raw, data )

    # Interpolate database
    fit_function_raw        = sp.interpolate.RegularGridInterpolator( ( beta, mu ), data )
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
    fit_props                   = FitProperties( )
    fit_props.dims              = 2
    fit_props.region_name       = "dGdttxx"
    fit_props.cheby_order_x     = 30
    fit_props.cheby_order_y     = 30
    fit_props.fcn_log_scale     = False
    fit_props.x_log_scale       = False
    fit_props.x_hpatch_np       = 30
    fit_props.x_max             = 30.0
    fit_props.x_min             = 0.0
    fit_props.y_log_scale       = True
    fit_props.y_max             = np.log10( 0.999 )
    fit_props.y_hpatch_np       = 4
    fit_props.y_min             = -4.0
    fit_props.cheby_abs_tol     = 1E+1
    fit_props.cheby_abs_tol_f   = 2
    fit_props.cheby_rel_tol     = 1E-3
    fit_props.cheby_rel_tol_f   = 2
    fit_props.x_map_fcn         = lambda x: x
    fit_props.y_map_fcn         = lambda y: y
    fit_props.max_ref_level     = 10
    fit_props.alpha_shift       = np.array( [ 1e-3, 1e-2, 1e-1 ] )

    fit_props.num_x             = fit_props.cheby_order_x + 1
    fit_props.num_x_fit         = fit_props.cheby_order_x + 1
    fit_props.num_y             = fit_props.cheby_order_y + 1
    fit_props.num_y_fit         = fit_props.cheby_order_y + 1

    fit_props.generate_fitting_matrix( )

    # Load database
    fipath = os.path.join( get_integrals_database_fopath( ), "1_time_domain", "Gttxx.h5" )
    with h5py.File( fipath, "r" ) as fid:
        mu          = fid[ "mu" ][:]
        beta        = fid[ "beta" ][:]
        data_raw    = fid[ "fcn" ][:]

    pos         = beta < fit_props.x_max+1e-6
    beta        = beta[ pos ]
    data_raw    = data_raw[ pos, : ]

    # Define G0 interpolation properties
    target_fcn      = dGdttxx_G0
    interp_raw      = sp.interpolate.RegularGridInterpolator( ( beta, mu ), data_raw )

    define_G0_fit_properties( fit_props, 100, mu, points_dist="linear" )

    # Calculate residual function in between G and G0*
    data, chfhpa     = calculate_residual_G0( fit_props, target_fcn, interp_raw, data_raw, beta, mu )

    # Plot residual function
    plot_residual_function( beta, mu, data_raw, data )

    # Interpolate database
    fit_function_raw        = sp.interpolate.RegularGridInterpolator( ( beta, mu ), data )
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


def generate_dGdt( show_summary_fig=False, show_figs=False )->None:
    # Define case type
    case_type   = "dGdt"

    # Define folder path to storage fit results
    fopath      = os.path.join( get_integrals_database_fit_fopath( ), "1_time_domain" )
    fopath_fit  = os.path.join( fopath, f"{case_type}_fit_results" )

    # Check if fit results folder exits
    if not os.path.isdir( fopath_fit ):
        os.makedirs( fopath_fit )

    # Fit coefficients
    fit_region, chfhpa = fit_residual_dGdt( fopath_fit, show_summary_fig=show_summary_fig, show_figs=show_figs )

    # Write coefficients
    write_coeffs_module_adaptive_2d_only_header( fit_region, fopath, case_type, is_time=True )

    for i,mi in enumerate( chfhpa.get_interpolators( ) ):
        write_coeffs_module_adaptive_1d_only_header( mi.get_fit_fcn( ), fopath, f"{case_type}A{i:d}" )


def generate_dGdtx( show_summary_fig=False, show_figs=False )->None:
    # Define case type
    case_type   = "dGdtx"

    # Define folder path to storage fit results
    fopath      = os.path.join( get_integrals_database_fit_fopath( ), "1_time_domain" )
    fopath_fit  = os.path.join( fopath, f"{case_type}_fit_results" )

    # Check if fit results folder exits
    if not os.path.isdir( fopath_fit ):
        os.makedirs( fopath_fit )

    # Fit coefficients
    fit_region, chfhpa = fit_residual_dGdtx( fopath_fit, show_summary_fig=show_summary_fig, show_figs=show_figs )

    # Write coefficients
    write_coeffs_module_adaptive_2d_only_header( fit_region, fopath, case_type, is_time=True )

    for i,mi in enumerate( chfhpa.get_interpolators( ) ):
        write_coeffs_module_adaptive_1d_only_header( mi.get_fit_fcn( ), fopath, f"{case_type}A{i:d}" )


def generate_dGdtxx( show_summary_fig=False, show_figs=False )->None:
    # Define case type
    case_type   = "dGdtxx"

    # Define folder path to storage fit results
    fopath      = os.path.join( get_integrals_database_fit_fopath( ), "1_time_domain" )
    fopath_fit  = os.path.join( fopath, f"{case_type}_fit_results" )

    # Check if fit results folder exits
    if not os.path.isdir( fopath_fit ):
        os.makedirs( fopath_fit )

    # Fit coefficients
    fit_region, chfhpa = fit_residual_dGdtxx( fopath_fit, show_summary_fig=show_summary_fig, show_figs=show_figs )

    # Write coefficients
    write_coeffs_module_adaptive_2d_only_header( fit_region, fopath, case_type, is_time=True )

    for i,mi in enumerate( chfhpa.get_interpolators( ) ):
        write_coeffs_module_adaptive_1d_only_header( mi.get_fit_fcn( ), fopath, f"{case_type}A{i:d}" )


def generate_dGdtt( show_summary_fig=False, show_figs=False )->None:
    # Define case type
    case_type   = "dGdtt"

    # Define folder path to storage fit results
    fopath      = os.path.join( get_integrals_database_fit_fopath( ), "1_time_domain" )
    fopath_fit  = os.path.join( fopath, f"{case_type}_fit_results" )

    # Check if fit results folder exits
    if not os.path.isdir( fopath_fit ):
        os.makedirs( fopath_fit )

    # Fit coefficients
    fit_region, chfhpa = fit_residual_dGdtt( fopath_fit, show_summary_fig=show_summary_fig, show_figs=show_figs )

    # Write coefficients
    write_coeffs_module_adaptive_2d_only_header( fit_region, fopath, case_type, is_time=True )

    for i,mi in enumerate( chfhpa.get_interpolators( ) ):
        write_coeffs_module_adaptive_1d_only_header( mi.get_fit_fcn( ), fopath, f"{case_type}A{i:d}" )


def generate_dGdttx( show_summary_fig=False, show_figs=False )->None:
    # Define case type
    case_type   = "dGdttx"

    # Define folder path to storage fit results
    fopath      = os.path.join( get_integrals_database_fit_fopath( ), "1_time_domain" )
    fopath_fit  = os.path.join( fopath, f"{case_type}_fit_results" )

    # Check if fit results folder exits
    if not os.path.isdir( fopath_fit ):
        os.makedirs( fopath_fit )

    # Check if fit results folder exits
    if not os.path.isdir( fopath_fit ):
        os.makedirs( fopath_fit )

    # Fit coefficients
    fit_region, chfhpa = fit_residual_dGdttx( fopath_fit, show_summary_fig=show_summary_fig, show_figs=show_figs )

    # Write coefficients
    write_coeffs_module_adaptive_2d_only_header( fit_region, fopath, case_type, is_time=True )

    for i,mi in enumerate( chfhpa.get_interpolators( ) ):
        write_coeffs_module_adaptive_1d_only_header( mi.get_fit_fcn( ), fopath, f"{case_type}A{i:d}" )


def generate_dGdttxx( show_summary_fig=False, show_figs=False )->None:
    # Define case type
    case_type   = "dGdttxx"

    # Define folder path to storage fit results
    fopath      = os.path.join( get_integrals_database_fit_fopath( ), "1_time_domain" )
    fopath_fit  = os.path.join( fopath, f"{case_type}_fit_results" )

    # Check if fit results folder exits
    if not os.path.isdir( fopath_fit ):
        os.makedirs( fopath_fit )

    # Fit coefficients
    fit_region, chfhpa = fit_residual_dGdttxx( fopath_fit, show_summary_fig=show_summary_fig, show_figs=show_figs )

    # Write coefficients
    write_coeffs_module_adaptive_2d_only_header( fit_region, fopath, case_type, is_time=True )

    for i,mi in enumerate( chfhpa.get_interpolators( ) ):
        write_coeffs_module_adaptive_1d_only_header( mi.get_fit_fcn( ), fopath, f"{case_type}A{i:d}" )


def get_root_fopath( ) -> str:
    """
    Get the absolute path to the project root directory.

    Computes the parent directory of this file's directory and returns its
    absolute path. Useful for building project-relative file paths.

    Returns
    -------
    str
        Absolute path to the project root directory.

    Examples
    --------
    >>> root = get_root_fopath()
    >>> isinstance(root, str)
    True
    """
    return os.path.abspath( os.path.join( os.path.dirname( __file__ ), "../" ) )


def get_integrals_database_fopath( ) -> str:
    """
    Get the absolute path to the integrals database folder.

    Uses the project root (via get_root_fopath) and appends 'aux_data/0_integrals_database'
    to build the absolute path to the folder that stores the raw integrals database.

    Returns
    -------
    str
        Absolute path to the integrals database directory.

    Examples
    --------
    >>> path = get_integrals_database_fopath()
    >>> os.path.basename(path)
    '0_integrals_database'
    """

    # Get repository root path
    root_path = get_root_fopath( )

    # Compose integrals database folder path
    db_fopath = os.path.join( root_path, "aux_data", "0_integrals_database" )

    return db_fopath


def get_integrals_database_fit_fopath( ) -> str:
    """
    Get the absolute path to the fitted integrals database folder.

    Uses the project root (via get_root_fopath) and appends 'aux_data/1_fit_integrals_database'
    to build the absolute path to the folder that stores the fitted (post-processed)
    integrals database, e.g., Chebyshev/HP-adaptive coefficients.

    Returns
    -------
    str
        Absolute path to the fitted integrals database directory.

    Examples
    --------
    >>> path = get_integrals_database_fit_fopath()
    >>> os.path.basename(path)
    '1_fit_integrals_database'
    """

    # Get repository root path
    root_path = get_root_fopath( )

    # Compose integrals database folder path
    db_fopath = os.path.join( root_path, "aux_data", "1_integrals_database_fit" )

    return db_fopath


def G_G0( beta: np.ndarray, mu ) -> np.ndarray:
    """
    Compute the G0 asymptotic Green's function for infinite water depth.
    
    This function calculates the asymptotic Green's function using Bessel functions
    of fractional order. The result represents the leading-order asymptotic behavior
    of the time-domain Green's function for wave-body interaction problems.
    
    Parameters
    ----------
    beta : np.ndarray
        1D array of temporal coordinate values (dimensionless time parameter).
    mu : float
        Ration in between the vertical coordinate and the radial distance (z/r).
    
    Returns
    -------
    np.ndarray
        1D array of G0 values evaluated at the given beta values, shape (n_beta,).
    
    Notes
    -----
    The function computes:
        G0 = (π * β² / (8√2)) * exp(-β² * μ / 4) * 
             [J_{0.25}(β²/8) * J_{0.75}(β²/8) - J_{-0.25}(β²/8) * J_{-0.75}(β²/8)]
    
    where J_ν denotes the Bessel function of the first kind of order ν.
    
    This formulation is based on: "Wehausen and Laitone, Surface Waves, (1960)".
    
    Examples
    --------
    >>> beta = np.linspace(0, 10, 100)
    >>> mu = 0.5
    >>> G0_vals = G_G0(beta, mu)
    >>> plt.plot(beta, G0_vals)
    >>> plt.xlabel('beta')
    >>> plt.ylabel('G0')
    >>> plt.show()
    
    See Also
    --------
    dGdt_G0 : Time derivative of G0
    dGdtx_G0 : Mixed derivative with respect to time and spatial coordinate
    """

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


def dGdt_G0( beta: np.ndarray, mu: float, alpha=0.0 ) -> np.ndarray:
    """
    Time derivative of G_G0 for infinite water depth.

    Computes ∂G0/∂t evaluated at the dimensionless time variable beta (β). This is the
    time derivative of G_G0(beta, mu) and is expressed using fractional-order Bessel
    functions. An optional phase shift alpha can be added inside the Bessel arguments
    to improve numerical robustness near small arguments.

    Parameters
    ----------
    beta : np.ndarray
        1D array of dimensionless time values β at which to evaluate the derivative.
    mu : float
        Parameter in [0, 1] related to vertical/radial geometry (e.g., z/r or cos(theta)).
    alpha : float, optional
        Small phase shift added to the Bessel arguments to stabilize the evaluation.
        Default is 0.0.

    Returns
    -------
    np.ndarray
        1D array of dG0/dt values evaluated at the provided beta points.

    Notes
    -----
    The implementation evaluates:
        dG0/dt = (π β^3 / (16√2)) exp(-β^2 μ / 4)
                 × [ J_{1/4}(x) J_{-1/4}(x) + J_{3/4}(x) J_{-3/4}(x) ],
    with x = β^2/8 + α and J_ν the Bessel function of the first kind.
    NaNs are zeroed and values for very small β (β < 1e-1) are set to zero to avoid
    numerical issues.

    Examples
    --------
    >>> beta = np.linspace(0, 10, 200)
    >>> mu = 0.5
    >>> dgdt = dGdt_G0(beta, mu)
    >>> dgdt_shifted = dGdt_G0(beta, mu, alpha=1e-3)
    """

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


def dGdtx_G0( beta: np.ndarray, mu: float, alpha=0.0 ) -> np.ndarray:
    """
    Mixed derivative ∂²G0/(∂μ ∂t) for infinite water depth.

    Computes the derivative with respect to μ of the time derivative of G_G0(β, μ).
    Since G_G0 depends on μ only through exp(-β² μ / 4), the μ-derivative is a
    multiplicative factor:
        ∂/∂μ [∂G0/∂t] = (-β²/4) * (∂G0/∂t).

    Parameters
    ----------
    beta : np.ndarray
        1D array of dimensionless time values β.
    mu : float
        Parameter in [0, 1] related to vertical/radial geometry (e.g., z/r or cos(theta)).
    alpha : float, optional
        Small phase shift forwarded to dGdt_G0 to stabilize Bessel evaluations. Default 0.0.

    Returns
    -------
    np.ndarray
        1D array of ∂²G0/(∂μ ∂t) evaluated at the provided beta points.

    Notes
    -----
    This routine leverages dGdt_G0 and applies the analytic factor (-β²/4) from
    differentiating the exp(-β² μ / 4) dependence with respect to μ.

    Examples
    --------
    >>> beta = np.linspace(0, 10, 200)
    >>> mu = 0.5
    >>> d2 = dGdtx_G0(beta, mu)
    """
    return ( -beta**2.0 / 4.0 ) * dGdt_G0( beta, mu, alpha=alpha )


def dGdtx_G0_num( beta: np.ndarray, mu: float, alpha=0.0 ) -> np.ndarray:
    eps = 1e-6
    return ( dGdt_G0( beta, mu+eps, alpha=alpha ) - dGdt_G0( beta, mu-eps, alpha=alpha ) ) / 2.0 / eps


def dGdtxx_G0( beta: np.ndarray, mu: float, alpha=0.0 ) -> np.ndarray:
    """
    Second derivative with respect to μ of the time derivative of G_G0.

    Computes ∂³G0/(∂μ² ∂t) for infinite water depth. Since G_G0 depends on μ only
    through exp(-β² μ / 4), each μ-derivative contributes a multiplicative factor
    (-β²/4). Therefore:
        dGdtxx_G0 = (-β²/4) * dGdtx_G0 = ((-β²/4)²) * dGdt_G0.

    Parameters
    ----------
    beta : np.ndarray
        1D array of dimensionless time values β.
    mu : float
        Geometry parameter in [0, 1] (e.g., z/r or cos(theta)).
    alpha : float, optional
        Phase shift forwarded to the underlying Bessel-based evaluations. Default 0.0.

    Returns
    -------
    np.ndarray
        1D array of ∂³G0/(∂μ² ∂t) evaluated at the provided beta points.

    Examples
    --------
    >>> beta = np.linspace(0, 10, 200)
    >>> mu = 0.5
    >>> d3 = dGdtxx_G0(beta, mu)
    """
    
    return ( -beta**2.0 / 4.0 ) * dGdtx_G0( beta, mu, alpha=alpha )


def dGdtt_G0( beta: np.ndarray, mu: float, alpha=0.0 ) -> np.ndarray:
    """
    Second time derivative of G_G0 for infinite water depth.

    Computes ∂²G0/∂t² evaluated at the dimensionless time variable beta (β). This is the
    second derivative in time of the asymptotic Green's function G_G0(β, μ). The
    implementation uses analytic expressions involving fractional-order Bessel functions
    and an optional phase shift α inside the Bessel arguments to improve robustness for
    small arguments.

    Parameters
    ----------
    beta : np.ndarray
        1D array of dimensionless time values β at which to evaluate the second derivative.
    mu : float
        Geometry parameter in [0, 1] related to vertical/radial ratio (e.g., z/r or cosθ).
    alpha : float, optional
        Small phase shift added to the Bessel arguments (x = β²/8 + α) to stabilize the
        evaluation near x ≈ 0. Default is 0.0.

    Returns
    -------
    np.ndarray
        1D array of ∂²G0/∂t² values evaluated at the provided beta points.

    Notes
    -----
    The computation factors the exponential dependence on μ and combines two terms:
      - T1: contribution from differentiating the exponential envelope
      - T2: contribution from differentiating the Bessel products
    with the overall scaling lt = π / (16√2) and exp factor exp(-β² μ / 4).

    Numerical guards are applied:
      - x < 1e-6 is set to 0
      - NaNs in intermediate Bessel products are patched
      - for very small β (< 5e-2), the result is linearly ramped to avoid singular behavior

    Examples
    --------
    >>> beta = np.linspace(0, 10, 200)
    >>> mu = 0.5
    >>> d2 = dGdtt_G0(beta, mu)
    """

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
    
    pos     = np.isnan( y0 )
    if pos.sum( ) > 0:
        num_pos = np.where( pos )[0][-1] + 1
        y0[:num_pos] = y0[num_pos]
    
    T1      = - 0.5 * a * y0 * beta * mu * expt

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

    pos     = np.isnan( yt )
    yt[pos] = 0.0
    dt      = 5e-2
    pos     = beta < dt
    num_pos = np.where( pos )[0][-1]
    yt[pos] = yt[num_pos] / dt * beta[pos]
    
    T2      = yt * expt

    return lt * ( T1 + T2 )


def dGdttx_G0( beta: np.ndarray, mu: float, alpha=0.0 ) -> np.ndarray:
    """
    Mixed derivative ∂³G0/(∂μ ∂t²) for infinite water depth.

    Computes the derivative with respect to μ of the second time derivative of
    G_G0(β, μ). Since G_G0 depends on μ only through exp(-β² μ / 4), the μ-derivative
    contributes an analytic multiplicative factor:
        ∂/∂μ [∂²G0/∂t²] = (-β²/4) * (∂²G0/∂t²).

    Parameters
    ----------
    beta : np.ndarray
        1D array of dimensionless time values β.
    mu : float
        Geometry parameter in [0, 1] (e.g., z/r or cos(theta)).
    alpha : float, optional
        Phase shift forwarded to dGdtt_G0 to stabilize Bessel evaluations. Default 0.0.

    Returns
    -------
    np.ndarray
        1D array of ∂³G0/(∂μ ∂t²) evaluated at the provided beta points.

    Notes
    -----
    Implemented as (-β²/4) * dGdtt_G0(beta, μ, α).

    Examples
    --------
    >>> beta = np.linspace(0, 10, 200)
    >>> mu = 0.5
    >>> d3_mix = dGdttx_G0(beta, mu)
    """
    return ( -beta**2.0 / 4.0 ) * dGdtt_G0( beta, mu, alpha=alpha )


def dGdttxx_G0( beta: np.ndarray, mu: float, alpha=0.0 ) -> np.ndarray:
    """
    Mixed derivative ∂⁴G0/(∂μ² ∂t²) for infinite water depth.

    Computes the second derivative with respect to μ of the second time derivative of
    G_G0(β, μ). Since G_G0 depends on μ only through exp(-β² μ / 4), each μ-derivative
    contributes a multiplicative factor (−β²/4). Therefore:
        dGdttxx_G0 = (−β²/4) * dGdttx_G0 = ((−β²/4)²) * dGdtt_G0.

    Parameters
    ----------
    beta : np.ndarray
        1D array of dimensionless time values β.
    mu : float
        Geometry parameter in [0, 1] (e.g., z/r or cos(theta)).
    alpha : float, optional
        Phase shift forwarded to the underlying Bessel-based evaluations. Default 0.0.

    Returns
    -------
    np.ndarray
        1D array of ∂⁴G0/(∂μ² ∂t²) evaluated at the provided beta points.

    Notes
    -----
    Implemented via the analytic μ-dependence factor as (−β²/4) times dGdttx_G0.
    """
    return ( -beta**2.0 / 4.0 ) * dGdttx_G0( beta, mu, alpha=alpha )


def map_to_interval( a: float, b: float, t: np.ndarray ) -> None:
    """
    Map points from the canonical interval [-1, 1] to a target interval [a, b].

    Applies the affine transform x = 0.5 * ((1 - t) * a + (t + 1) * b) elementwise
    to t, where t is assumed to lie in [-1, 1] (e.g., Chebyshev nodes).

    Parameters
    ----------
    a : float
        Lower bound of the target interval.
    b : float
        Upper bound of the target interval.
    t : np.ndarray
        Points in the canonical interval [-1, 1] to be mapped.

    Returns
    -------
    np.ndarray
        Points mapped to the interval [a, b] with the same shape as t.

    Examples
    --------
    >>> t = np.array([-1.0, 0.0, 1.0])
    >>> map_to_interval(0.0, 10.0, t)
    array([ 0.,  5., 10.])
    """
    return 0.5 * ( ( 1 - t ) * a + ( t + 1 ) * b )


def map_to_unit( a: float, b: float, x: np.ndarray ) -> None:
    """
    Map points from an interval [a, b] to the canonical interval [-1, 1].

    Applies the inverse affine transformation to map values from a given interval
    to the standard Chebyshev domain [-1, 1]. This is the inverse operation of
    map_to_interval.

    Parameters
    ----------
    a : float
        Lower bound of the source interval.
    b : float
        Upper bound of the source interval.
    x : np.ndarray
        Points in the interval [a, b] to be mapped to [-1, 1].

    Returns
    -------
    np.ndarray
        Points mapped to the canonical interval [-1, 1] with the same shape as x.

    Notes
    -----
    The transformation is: t = (2*x - (a + b)) / (b - a)
    
    This ensures that:
    - x = a maps to t = -1
    - x = b maps to t = 1
    - x = (a+b)/2 maps to t = 0

    Examples
    --------
    >>> x = np.array([0.0, 5.0, 10.0])
    >>> map_to_unit(0.0, 10.0, x)
    array([-1.,  0.,  1.])
    
    See Also
    --------
    map_to_interval : Inverse transformation from [-1, 1] to [a, b]
    """

    return ( 2.0 * x - ( a + b ) ) / ( b - a )


def plot_residual_function( beta: np.ndarray, mu: np.ndarray, data_raw: np.ndarray, data: np.ndarray ) -> None:
    # Generate grid for plotting
    B, M = np.meshgrid( beta, mu, indexing="ij" )

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


if __name__ == "__main__":
    # generate_dGdt( )
    generate_dGdtx( )
    # generate_dGdtxx( )
    # generate_dGdtt( )
    # generate_dGdttx( )
    # generate_dGdttxx( )