
# Import general usage libraries
import h5py
import os

# Import general usage scientific libraries
import numpy as np
import scipy as sp

# Import general usage scientific plotting libraries
import bokeh.plotting as blt
import matplotlib.pyplot as plt

# Import local modules
from fit_cheby import fit_integral_1d, FitProperties


class PolyInt:

    def __init__( self, fcn, x_min, x_max, deg=2 ) -> None:
        # Define class attributes
        self.deg        = deg
        self.fcn        = fcn
        self.poly       = np.array( [ ] )
        self.poly_int   = np.array( [ ] )
        self.x_max      = x_max
        self.x_min      = x_min

        # Fit function
        self.fit( )

    def calculate_fcn( self, t ) -> None:
        return np.polyval( self.poly, t )
    
    def calculate_fcn_primitive( self, t ) -> None:
        return np.polyval( self.poly_int, t )

    def fit( self ) -> None:
        # Fit polynomial to curve
        x           = np.linspace( self.x_min, self.x_max, 100 )
        y           = self.fcn( x )
        self.poly   = np.polyfit( x, y, self.deg )

        # Calcualte integrated polynomial
        self.poly_int       = np.zeros( ( self.poly.shape[0] + 1, ) )
        self.poly_int[:-1]  = self.poly / np.linspace( self.poly.shape[0], 1, self.poly.shape[0] )

    def get_fit_curve( self, dnp: int ) -> None:
        x   = np.linspace( self.x_min, self.x_max, dnp )
        fa  = self.calculate_fcn( x )

        return x, fa
    
    def get_primitive( self, dnp: int ) -> None:
        x   = np.linspace( self.x_min, self.x_max, dnp )
        fa  = self.calculate_fcn_primitive( x )

        return x, fa
    
    def get_primitive_end( self ) -> float:
        return self.calculate_fcn_primitive( self.x_max )

    def get_primitive_start( self ) -> float:
        return self.calculate_fcn_primitive( self.x_min )
    
    def plot_fit( self ) -> None:
        # Get fit curve
        x, fc = self.get_fit_curve( 100 )

        # Get original curve
        fo = self.fcn( x )

        # Plot functions
        fig = plt.figure( )
        ax0 = fig.add_subplot( 211 )
        ax1 = fig.add_subplot( 212 )

        ax0.plot( x, fo )
        ax0.plot( x, fc )

        ax1.plot( x, np.log10( np.abs( fo - fc ) ) )

        plt.show( )


class PolyIntCheby:

    def __init__( self, fcn, x_min, x_max ) -> None:
        # Define class attributes
        self.fcn        = fcn
        self.fit_props  = FitProperties( )
        self.x_max      = x_max
        self.x_min      = x_min

        # Fit function
        self.fit( )

    def calculate_fcn( self, t: np.ndarray ) -> np.ndarray:
        return self.c0 + self.c1 * t + self.c2 * t**2.0
    
    def calculate_fcn_primitive( self, t: np.ndarray ) -> np.ndarray:
        return self.c0 * t + self.c1 * t**2.0 + self.c2 * t**3.0

    def fit( self ) -> None:
        # Fit function using Chebyshev polynomials
        self.fit_props.x_min            = self.x_min
        self.fit_props.x_max            = self.x_max
        self.fit_props.cheby_order_x    = 30
        self.fit_props.cheby_tol        = 1e-6
        fcn                             = lambda x: self.fcn( x )
        self.cx, self.nx                = fit_integral_1d( fcn, self.fit_props, "Fcn Fit", show_stats=False, create_figs=False, show_figs=False )

        # Define collapsed coefficients
        self.c0                         = self.cx[0] - self.cx[2]
        self.c1                         = self.cx[1]
        self.c2                         = 2.0 * self.cx[2]

    def get_fit_curve( self, dnp: int ) -> None:
        t                       = np.linspace( -1.0, 1.0, dnp )
        x                       = self.fit_props.x_map_lin( t )
        fa                      = self.calculate_fcn( t )

        return x, fa
    
    def get_primitive( self, dnp: int ) -> None:
        t                       = np.linspace( -1.0, 1.0, dnp )
        x                       = self.fit_props.x_map_lin( t )
        fa                      = self.calculate_fcn_primitive( t )

        return x, fa
    
    def get_primitive_end( self ) -> float:
        return self.calculate_fcn_primitive( 1.0 )

    def get_primitive_start( self ) -> float:
        return self.calculate_fcn_primitive( -1.0 )
    
    def plot_fit( self ) -> None:
        # Get fit curve
        x, fc = self.get_fit_curve( 100 )

        # Get original curve
        fo = self.fcn( x )

        # Plot functions
        fig = plt.figure( )
        ax0 = fig.add_subplot( 211 )
        ax1 = fig.add_subplot( 212 )

        ax0.plot( x, fo )
        ax0.plot( x, fc )

        ax1.plot( x, np.log10( np.abs( fo - fc ) ) )

        plt.show( )


def fit_besselj(  ) -> None:
    fit_props               = FitProperties( )
    fit_props.x_min         = 1.0
    fit_props.x_max         = 1.1
    fit_props.cheby_order_x = 30
    fit_props.cheby_tol     = 1e-6
    fcn                     = lambda x: sp.special.jv(0, x)
    cx, nx                  = fit_integral_1d( fcn, fit_props, "BesselJ", create_figs=True, show_figs=True )

    t                       = np.linspace( -1.0, 1.0, 100 )
    x                       = fit_props.x_map_lin( t )
    c0                      = cx[0] - cx[2]
    c1                      = cx[1]
    c2                      = 2.0 * cx[2]
    fa                      = c0 + c1 * t + c2 * t**2.0

    fig = plt.figure( )
    ax0 = fig.add_subplot( 211 )
    ax1 = fig.add_subplot( 212 )
    ax0.plot( x, fcn( x ) )
    ax0.plot( x, fa )
    ax1.plot( x, fcn( x ) - fa )
    plt.show( )
    print( nx.shape[0] )


def fit_besselj_2(  ) -> None:
    dnp = 9
    F0 = np.zeros( ( dnp, ) )
    F1 = np.zeros( ( dnp, ) )
    for i in range( dnp ):
        # Fit ascending polynomial to sin(x)
        a = i*0.1
        b = (i+1)*0.1
        x = np.linspace( a, b, 100 )
        # y = np.sqrt( x )
        y = -sp.special.jv( 1, x )
        # y = np.sin( x ) * np.exp( -x ) * sp.special.jv( 0, x )
        # y = np.sqrt( x ) * np.sin( np.sqrt( x ) ) * np.exp( -x ) * sp.special.jv( 0, x )
        p = np.polyfit( x, y, 6 )
        yi = np.polyval( p, x )

        # Print fitted polynomial coefficients
        print( f"p[{i:d}]:", p )

        # Create integrated polynomial coefficients
        p0i      = np.zeros( ( p.shape[0]+1, ) )
        p0i[:-1] = p / np.linspace( p.shape[0], 1, p.shape[0] )

        # Evaluate polynomial integral
        f0 = np.polyval( p0i, a )
        f1 = np.polyval( p0i, b )

        F0[i] = f0
        F1[i] = f1

        print( f"Integral Value[{i:d}]: ", f1 - f0 )

    print( "Total integral Value: ", ( F1 - F0 ).sum( ) )
    print( "Total integral Value 2: ", F1[-1]-F0[0] )
    
    fig = plt.figure( )
    ax0 = fig.add_subplot( 211 )
    ax1 = fig.add_subplot( 212 )
    ax0.plot( F0, "-o" )
    ax0.plot( F1, "-o" )
    ax1.plot( np.log10( np.abs( F1[:-1] - F0[1:] ) ) )
    plt.show( )

    # # Plot fit
    # fig = plt.figure( )
    # ax0 = fig.add_subplot( 311 )
    # ax1 = fig.add_subplot( 312 )
    # ax2 = fig.add_subplot( 313 )
    # ax0.plot( x0, y0 )
    # ax0.plot( x0, y0i )
    # ax1.plot( x0, y0i - y0 )
    # ax2.plot( x0, sp.special.jv( 0, x0 ) )
    # ax2.plot( x0, np.polyval( p0i, x0 ) + 1 )

    # plt.show( )


def fit_exp(  ) -> None:
    # Fit ascending polynomial to sin(x)
    x = np.linspace( 0.0, 1.0, 100 )
    y = np.exp( -x )
    p = np.polyfit( x, y, 6 )
    yi = np.polyval( p, x )

    # Print fitted polynomial coefficients
    print( "p:", p )

    # Create integrated polynomial coefficients
    pi      = np.zeros( ( p.shape[0]+1, ) )
    pi[:-1] = p / np.linspace( p.shape[0], 1, p.shape[0] )

    # Evaluate polynomial integral
    F0      = np.polyval( pi, 0 )
    F1      = np.polyval( pi, 1.0 )

    print( "Integral Value: ", F1 - F0, - F0 + 1.0 )
    print( F0, F1 )

    # Plot fit
    fig = plt.figure( )
    ax0 = fig.add_subplot( 311 )
    ax1 = fig.add_subplot( 312 )
    ax2 = fig.add_subplot( 313 )
    ax0.plot( x, y )
    ax0.plot( x, yi )
    ax0.plot( x, np.polyval( pi, x ) )
    ax1.plot( x, yi - y )
    ax2.plot( x, -np.exp( -x ) )
    ax2.plot( x, np.polyval( pi, x ) - 1.0 )

    plt.show( )


def fit_exp_cos(  ) -> None:
    # Fit ascending polynomial to sin(x)
    x = np.linspace( 0.0, 1.0, 100 )
    y = np.exp( -x ) * np.cos( x )
    p = np.polyfit( x, y, 6 )
    yi = np.polyval( p, x )

    # Print fitted polynomial coefficients
    print( "p:", p )

    # Create integrated polynomial coefficients
    pi      = np.zeros( ( p.shape[0]+1, ) )
    pi[:-1] = p / np.linspace( p.shape[0], 1, p.shape[0] )

    # Evaluate polynomial integral
    F0      = np.polyval( pi, 0 )
    F1      = np.polyval( pi, 1.0 )

    print( "Integral Value: ", F1 - F0, - F0 + 0.5 )
    print( F0, F1 )

    # Plot fit
    fig = plt.figure( )
    ax0 = fig.add_subplot( 311 )
    ax1 = fig.add_subplot( 312 )
    ax2 = fig.add_subplot( 313 )
    ax0.plot( x, y )
    ax0.plot( x, yi )
    ax1.plot( x, yi - y )
    ax2.plot( x, np.exp( -x )*( np.sin( x ) - np.cos( x ) ) / 2.0 )
    ax2.plot( x, np.polyval( pi, x ) - 0.5 )

    plt.show( )


def fit_exp_series(  ) -> None:
    C = 0.1
    x = np.linspace( 0, 1.0, 10 )
    y = np.zeros( x.shape )
    yi = np.zeros( x.shape )
    yo = np.exp( C * x )
    yio = np.exp( C * x ) / C

    for i in range( x.shape[0] ):
        cc = 1
        nf = 1
        xf = 1.0
        y[i] = 1.0
        for j in range( 1, 100 ):
            cc *= C
            nf *= j
            xf *= x[i]
            y[i] +=  cc * xf / nf

    for i in range( x.shape[0] ):
        cc = 1
        nf = 1
        xf = x[i]
        yi[i] = xf
        for j in range( 1, 100 ):
            cc *= C
            nf *= j
            xf *= x[i]
            yi[i] +=  cc * xf / nf / ( i + 1 )
        print( yi[i] )

    print( yi[-1]-yi[0], yio[-1]-yio[0] )

    fig = plt.figure( )
    ax0 = fig.add_subplot( 211 )
    ax1 = fig.add_subplot( 212 )
    ax0.plot( x, y )
    ax0.plot( x, yo )
    ax1.plot( x, yi+10 )
    ax1.plot( x, yio )
    plt.show( )


def fit_exp_sin(  ) -> None:
    # Fit ascending polynomial to sin(x)
    x = np.linspace( 0.0, 1.0, 100 )
    y = np.exp( -x ) * np.sin( x )
    p = np.polyfit( x, y, 6 )
    yi = np.polyval( p, x )

    # Print fitted polynomial coefficients
    print( "p:", p )

    # Create integrated polynomial coefficients
    pi      = np.zeros( ( p.shape[0]+1, ) )
    pi[:-1] = p / np.linspace( p.shape[0], 1, p.shape[0] )

    # Evaluate polynomial integral
    F0      = np.polyval( pi, 0 )
    F1      = np.polyval( pi, 1.0 )

    print( "Integral Value: ", F1 - F0, - F0 + 0.5 )
    print( F0, F1 )

    # Plot fit
    fig = plt.figure( )
    ax0 = fig.add_subplot( 311 )
    ax1 = fig.add_subplot( 312 )
    ax2 = fig.add_subplot( 313 )
    ax0.plot( x, y )
    ax0.plot( x, yi )
    ax1.plot( x, yi - y )
    ax2.plot( x, np.exp( -x )*( -np.sin( x ) - np.cos( x ) ) / 2.0 )
    ax2.plot( x, np.polyval( pi, x ) - 0.5 )

    plt.show( )


def fit_sin_ascending(  ) -> None:
    # Fit ascending polynomial to sin(x)
    x = np.linspace( 0.0, np.pi/2.0, 100 )
    y = np.sin( x )
    p = np.polyfit( x, y, 6 )
    yi = np.polyval( p, x )

    x0 = np.linspace( 0.0, np.pi/3.0, 100 )
    y0 = np.sin( x0 )
    p0 = np.polyfit( x0, y0, 6 )
    y0i = np.polyval( p0, x0 )

    x1 = np.linspace( np.pi/3.0, 2.0*np.pi/3.0, 100 )
    y1 = np.sin( x1 )
    p1 = np.polyfit( x1, y1, 6 )
    y1i = np.polyval( p1, x1 )

    x2 = np.linspace( 2.0*np.pi/3.0, np.pi, 100 )
    y2 = np.sin( x2 )
    p2 = np.polyfit( x2, y2, 8 )
    y2i = np.polyval( p2, x2 )

    # Print fitted polynomial coefficients
    print( "p:", p )
    print( "p0:", p0 )
    print( "p1:", p1 )

    # Create integrated polynomial coefficients
    pi      = np.zeros( ( p.shape[0]+1, ) )
    pi[:-1] = p / np.linspace( p.shape[0], 1, p.shape[0] )

    p0i      = np.zeros( ( p0.shape[0]+1, ) )
    p0i[:-1] = p0 / np.linspace( p0.shape[0], 1, p0.shape[0] )

    p1i      = np.zeros( ( p1.shape[0]+1, ) )
    p1i[:-1] = p1 / np.linspace( p1.shape[0], 1, p1.shape[0] )

    p2i      = np.zeros( ( p2.shape[0]+1, ) )
    p2i[:-1] = p2 / np.linspace( p2.shape[0], 1, p2.shape[0] )

    # Evaluate polynomial integral
    F0      = np.polyval( pi, 0 )
    F1      = np.polyval( pi, np.pi/2.0 )

    F00      = np.polyval( p0i, 0 )
    F10      = np.polyval( p0i, np.pi/3.0 )

    F01      = np.polyval( p1i, np.pi/3.0 )
    F11      = np.polyval( p1i, 2.0*np.pi/3.0 )

    F02      = np.polyval( p2i, 2.0*np.pi/3.0 )
    F12      = np.polyval( p2i, np.pi )

    print( "Integral Value: ", F1 - F0 )
    print( "Integral Value 1: ", ( F11 - F01 ) + ( F10 - F00 ), F11 - F00 )
    print( "Integral Value 2: ", ( F12 - F02 ) + ( F11 - F01 ) + ( F10 - F00 ), F12 - F00 )
    print( F0, F1 )
    print( F00, F10, F01, F11, F02, F12 )

    # Plot fit
    fig = plt.figure( )
    ax0 = fig.add_subplot( 311 )
    ax1 = fig.add_subplot( 312 )
    ax2 = fig.add_subplot( 313 )
    ax0.plot( x, y )
    ax0.plot( x, yi )
    ax0.plot( x, np.polyval( pi, x ) )
    ax1.plot( x, yi - y )
    ax2.plot( x, -np.cos( x ) )
    ax2.plot( x, np.polyval( pi, x ) - 1.0 )

    fig = plt.figure( )
    ax0 = fig.add_subplot( 211 )
    ax1 = fig.add_subplot( 212 )
    ax0.plot( x2, y2 )
    ax0.plot( x2, y2i )
    ax1.plot( x2, y2i - y2 )

    plt.show( )


def fit_tanh(  ) -> None:
    fit_props               = FitProperties( )
    fit_props.x_min         = 1.0
    fit_props.x_max         = 1.5
    fit_props.cheby_order_x = 30
    fit_props.cheby_tol     = 1e-5
    cx, nx                  = fit_integral_1d( lambda x: np.sqrt( x * np.tanh( x ) ), fit_props, "Srqt(x·tanh(x))", create_figs=True, show_figs=True )
    # cx, nx                  = fit_integral_1d( lambda x: np.sqrt( x * np.tanh( x ) ) / ( 1 + np.exp( -4*x ) ), fit_props, "Srqt(x·tanh(x))", create_figs=True, show_figs=True )
    print( cx )
    print( nx )
    print( nx.shape[0] )


def fit_tanh_2(  ) -> None:
    x = np.linspace( 0.0, 1.0, 100 )
    y = np.sqrt( x * np.tanh( x ) )
    p = np.polyfit( x, y, 10 )
    yi = np.polyval( p, x )

    print( "p: ", p )

    fig = plt.figure( )
    ax0 = fig.add_subplot( 211 )
    ax1 = fig.add_subplot( 212 )

    ax0.plot( x, y )
    ax0.plot( x, yi )

    ax1.plot( x, yi - y )

    plt.show( )


def integrate_besselj_poly(  ) -> None:
    # Define function to fit and integrate
    fcn         = lambda x: -sp.special.jv( 1, x )
    fcnp        = lambda x: sp.special.jv( 0, x )

    # Define start interval
    a           = 0.0
    b           = 0.01
    poly_int_0  = PolyInt( fcn, a, b )

    # Define end interval
    c           = 1.99
    d           = 2.0
    poly_int_1  = PolyInt( fcn, c, d )

    # Integrate function by looping over the intervals
    dnp             = 200
    val             = 0.0
    F0              = np.zeros( ( dnp, ) )
    F1              = np.zeros( ( dnp, ) )
    int_values_o    = np.zeros( ( dnp, ) )
    int_values_i    = np.zeros( ( dnp, ) )
    fig             = plt.figure( )
    ax0             = fig.add_subplot( 211 )
    ax1             = fig.add_subplot( 212 )
    for i in range( 0, dnp ):
        ai              = i*0.01
        bi              = (i+1)*0.01
        poly_int_i      = PolyInt( fcn, ai, bi )
        val             = ( poly_int_i.get_primitive_end( ) - poly_int_i.get_primitive_start( ) )
        int_values_i[i] = val

        F0[i]           = poly_int_i.get_primitive_start( )
        F1[i]           = poly_int_i.get_primitive_end( )

        xp, f           = poly_int_i.get_fit_curve( 100 )
        fo              = fcn( xp )

        xp, fp          = poly_int_i.get_primitive( 100 )
        fpo             = fcnp( xp )
        int_values_o[i] = fpo[-1] - fpo[0]

        ax0.plot( xp, f, "b" )
        ax0.plot( xp, fo, "r" )

        ax1.plot( xp, fp, "b" )
        ax1.plot( xp, fpo, "r" )

        print( i, ai, bi, poly_int_i.get_primitive_start( ), poly_int_i.get_primitive_end( ), val )

    # Integrate function using gauss-knord quadarature
    value, value_err    = integrate_gauss( fcn, a, d, 0.5 )

    # Integrate using polynomials
    value_poly_ends     = ( poly_int_1.get_primitive_end( ) - poly_int_0.get_primitive_start( ) )
    
    print( value, value_err, value_poly_ends, int_values_o.sum( ), int_values_i.sum( ) )

    fig = plt.figure( )
    ax0 = fig.add_subplot( 311 )
    ax1 = fig.add_subplot( 312 )
    ax2 = fig.add_subplot( 313 )

    ax0.plot( F0 )
    ax0.plot( F1 )

    ax1.plot( F1[:-1] - F0[1:] )
    
    ax2.plot( int_values_o, label="Orig" )
    ax2.plot( int_values_i, label="Interp" )

    plt.show( )


def integrate_gauss( fcn, a: float, b: float, dx: float ) -> float:
    value       = 0.0
    value_err   = 0.0
    n_iter      = int( np.floor( ( b - a ) / dx ) )

    for i in range( n_iter ):
        ri  = sp.integrate.quad(
                                    fcn,
                                    a+i*dx,
                                    a+(i+1)*dx
                                )
        value       += ri[0]
        value_err   += ri[1]

    if np.abs( ( b - a ) - ( n_iter * dx - a ) ) > 1e-3:
        ri  = sp.integrate.quad(
                                    fcn,
                                    a+n_iter*dx,
                                    b
                                )
        value       += ri[0]
        value_err   += ri[1]

    return value, value_err


def ma_filter( x, w ) -> np.ndarray:
    return np.convolve(x, np.ones( w ), 'same') / w


def main( fipath: str ) -> None:
    # Read results from disk
    with h5py.File( fipath, "r" ) as fid:
        A           = fid[ "A" ][:]
        B           = fid[ "B" ][:]
        T           = fid[ "T" ][:]
        fcn_err     = fid[ "fcn_err" ][:]
        fcn_time    = fid[ "fcn_time" ][:]
        fcn_val     = fid[ "fcn_val" ][:]

    bp = -1

    # amp     = 160.0
    # mu      = 1.5
    # sigma   = 2.7
    # lam     = 4.0
    # eta     = 8.0
    # xi      = 20.0
    # y = np.exp( - ( ( T - mu )**2.0 / sigma ) ) * amp * np.sin( 2.5 * T * T ) + np.exp( - ( ( T - eta )**1.1 / xi ) ) / lam

    # B_mat, T_mat = np.meshgrid( B, T )
    # fig = plt.figure( )
    # ax  = fig.add_subplot( 111, projection="3d" )

    # ax.plot_surface( B_mat.T, T_mat.T, fcn_val[0, :, :] )

    # plt.show( )

    # f = fcn_val[0, bp, :]
    # fd = np.diff( f )
    # pos = np.where( fd[1:] * fd[:-1] < 0 )[0] + 1
    # tm = ( T[pos[:-1]] + T[pos[1:]] ) / 2.0
    # fm = ( f[pos[:-1]] + f[pos[1:]] ) / 2.0

    # # fa = ma_filter( fcn_val[0, bp, :], 100 )
    # plt.plot( T, fcn_val[0, bp, :], label="Ref" )
    # plt.plot( tm, fm, "o" )
    # # plt.plot( T, fa, label="MA" )
    # # plt.plot( T, y, label="Model" )
    # plt.show( )
    # raise ValueError( "Stop by user" )

    tools   = "pan,wheel_zoom,box_zoom,hover,reset"
    fig     = blt.figure( 
                                title        = f"A: {A[0]:0.3f}",
                                x_axis_label = "T [-]",
                                y_axis_label = "G [-]",
                                sizing_mode  = "stretch_both",
                                tools        = tools
                            )
    
    for i in range( B.shape[0] ):
        fig.line( 
                    T, 
                    fcn_val[0, i, :], 
                    color=( int( i/B.shape[0]*255 ), 0, 0),
                    legend_label=f"B: {B[i]:0.3f}"
                )
        
    # for i in range( 20 ):
    #     fig.line( 
    #                 B, 
    #                 fcn_val[0, :, i], 
    #                 color=( int( i/20*255 ), 0, 0),
    #                 legend_label=f"T: {T[i]:0.3f}"
    #             )
        
    fig.legend.click_policy="hide"
    blt.show( fig )
    raise ValueError( "Stop by user" )
    # udc( G )
    step_np         = 5
    intervals_np    = T.shape[0] // step_np
    irnp            = T.shape[0] % step_np

    if irnp > 0:
        intervals                   = np.zeros( ( intervals_np + 2, ), dtype=int )
        intervals[:intervals_np+1]  = np.linspace( 0, intervals_np, intervals_np+1, dtype=int ) * step_np
        intervals[intervals_np]     = intervals[ intervals_np - 1 ] + irnp
        intervals[intervals_np+1]   = intervals[intervals_np] + step_np
        intervals_np                += 2
    else:
        intervals                   = np.linspace( 0, intervals_np, intervals_np+1, dtype=int ) * step_np
        intervals[intervals_np]     = intervals[intervals_np-1] + step_np
        intervals_np                += 1
    
    count_coeffs    = 0
    vec_np          = 30 * ( intervals_np - 1 )
    cum_points      = np.zeros( ( intervals_np, ), dtype=int )
    T_intv_vec      = np.zeros( ( intervals_np - 1, 2 ), )
    nx_vec          = np.zeros( ( vec_np, ) )
    cx_vec          = np.zeros( ( vec_np, ) )
    for i in range( intervals_np - 1 ):
        n0 = i * step_np
        n1 = ( i + 1 ) * step_np
        Tv = T[n0:n1]
        Fv = fcn_val[0, bp, n0:n1]
        f  = sp.interpolate.interp1d( Tv, Fv, kind="cubic" )
        TT = np.linspace( Tv[0], Tv[-1], Tv.shape[0]*10 )

        # plt.plot( T, fcn_val[0, -1, :], label="Raw" )
        # plt.plot( TT, f( TT ), label="Interpolation" )
        # plt.show( )

        fit_props               = FitProperties( )
        fit_props.x_min         = Tv[0]
        fit_props.x_max         = Tv[-1]
        fit_props.cheby_order_x = 30
        fit_props.cheby_tol     = 1e-6
        cx, nx                  = fit_integral_1d( f, fit_props, "A: 0.01 - B: 1.99", create_figs=False, show_figs=False )
        c1                      = count_coeffs + nx.shape[0]
        nx_vec[count_coeffs:c1] = nx
        cx_vec[count_coeffs:c1] = cx
        count_coeffs            = c1
        cum_points[i+1]         = count_coeffs
        T_intv_vec[i, 0]        = Tv[0]
        T_intv_vec[i, 1]        = Tv[-1]
        print( T_intv_vec[i, :] )

    nx_vec = nx_vec[:count_coeffs]
    cx_vec = cx_vec[:count_coeffs]

    fopath = r"E:\sergio\0050_OASIS_SM\TimeDomain"
    finame = f"test_coeffs_{step_np:d}.h5"
    fipath = os.path.join( fopath, finame )
    with h5py.File( fipath, "w" ) as fid:
        fid.attrs[ "step_np" ] = step_np
        fid.create_dataset( "T", data=T )
        fid.create_dataset( "T_intv", data=T_intv_vec )
        fid.create_dataset( "intervals", data=intervals )
        fid.create_dataset( "cum_points", data=cum_points )
        fid.create_dataset( "nx", data=nx_vec )
        fid.create_dataset( "cx", data=cx_vec )
        fid.create_dataset( "fcn_val", data=fcn_val[0, bp, :] )
    
    # plt.plot( T, G )
    # plt.show( )


def udc( t: np.ndarray, data: np.ndarray ) -> np.ndarray:
    # Remove sample mean
    points  = np.linspace( 0, data.shape[0]-1, data.shape[0] )
    data    -= data.mean( )

    # Detect zero crossing points
    pos_zc  = np.where( ( data[:-1] * data[1:] ) < 0.0 )[0]
    
    plt.plot( points, data, "*" )
    plt.plot( points[ pos_zc ], data[ pos_zc ], "o", color="red" )
    plt.plot( points[ pos_zc + 1 ], data[ pos_zc + 1 ], "o", color="green" )
    plt.show( )

    # 


if __name__ == "__main__":
    import sys
    file_path = r"E:\sergio\0050_OASIS_SM\TimeDomain\database.h5"

    print( sys.argv )

    # main( file_path )

    # fit_besselj(  )
    # fit_besselj_2(  )
    # fit_tanh(  )
    # fit_tanh_2(  )

    # fit_exp( )
    # fit_exp_cos( )
    # fit_exp_sin( )
    # fit_besselj_2( )
    # fit_exp_series( )
    # fit_sin_ascending( )

    # integrate_besselj_poly( )

    # T = 10.0
    # x = np.linspace( 0, 5.0, 100000 )
    # y = x * np.sin( T * x )
    # y1 = np.sqrt( x * np.tanh( x ) ) * np.sin( T * np.sqrt( x * np.tanh( x ) ) )
    # y2 = x * T * x

    # z = 1 + np.exp( -4*x )
    # z1 = 1 + 1 + 1/4*x

    # fig = plt.figure( )
    # ax0 = fig.add_subplot( 211 )
    # ax1 = fig.add_subplot( 212 )
    # ax0.plot( x, y )
    # ax0.plot( x, y1 )
    # ax0.plot( x, y2 )
    # ax1.plot( x, np.log( np.abs( y1 - y ) ) )
    # ax1.plot( x, np.log( np.abs( y1 - y2 ) ) )

    # fig = plt.figure( )
    # ax0 = fig.add_subplot( 211 )
    # ax1 = fig.add_subplot( 212 )
    # ax0.plot( x, z )
    # ax0.plot( x, z1 )

    # ax1.plot( x, np.log( np.abs( z1 - z ) ) )

    # plt.show( )

    # x = np.linspace( 0, 100, 10000 )
    # for i in range( 50 ):
    #     plt.plot( x, sp.special.jv( i, x ) )
    # plt.show( )
