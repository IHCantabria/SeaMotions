
# Import general usage libraries
import h5py
import os

# Import general usage scientific libraries
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy.special import jv


class IntegralModel:

    def __init__( self, beta: np.ndarray, mu: np.ndarray, data: np.ndarray, name: str ) -> None:
        # Storage class attributes
        self.B          = np.ndarray
        self.B_dm       = np.ndarray
        self.beta       = beta
        self.data       = data
        self.data_db    = np.ndarray
        self.data_dm    = np.ndarray
        self.fidata     = None
        self.M          = np.ndarray
        self.M_dm       = np.ndarray
        self.mu         = mu
        self.name       = name
        
        # Initialize data
        self._initialize( )

    def _initialize( self ) -> None:
        # Create data interpolator
        self.fdata = sp.interpolate.RegularGridInterpolator( ( self.beta, self.mu ), self.data )

        # Derivate data w.r.t mu
        self.data_db = ( ( self.data[2:, :] - self.data[:-2, :] ).T / ( self.beta[2:] - self.beta[:-2] ) ).T

        # Derivate data w.r.t mu
        self.data_dm = ( self.data[:, 2:] - self.data[:, :-2] ) / ( self.mu[2:] - self.mu[:-2] )

        # Create 2D domain supports
        self.B, self.M          = np.meshgrid( self.beta, self.mu, indexing="ij" )
        self.B_db, self.M_db    = np.meshgrid( self.beta[1:-1], self.mu, indexing="ij" )
        self.B_dm, self.M_dm    = np.meshgrid( self.beta, self.mu[1:-1], indexing="ij" )


def compare_discretization( gt: IntegralModel, gtd: IntegralModel, title, show=False ) -> None:
    # Create new evaluation matrix
    N       = 110
    beta    = np.linspace( gtd.beta[0], gtd.beta[-1], N )
    mu      = np.linspace( gtd.mu[0], gtd.mu[-1], N )
    B, M    = np.meshgrid( beta, mu, indexing="ij" )

    # Create figure layout
    fig, ax = plt.subplots( 2, 2 )
    fig.suptitle( title )

    # Plot reference data
    ax[0, 0].set_title( "Ref.Data" )
    cnf = ax[0, 0].contourf( gt.B, gt.M, gt.data )
    plt.colorbar( cnf, ax=ax[0, 0] )

    # Plot Downsampled data
    ax[0, 1].set_title( "Downsampled data" )
    cnf = ax[0, 1].contourf( gtd.B, gtd.M, gtd.data )
    plt.colorbar( cnf, ax=ax[0, 1] )

    # Plot difference in between Reference and downsampled data
    ax[1, 0].set_title( "Difference: Ref-Down" )
    cnf = ax[1, 0].contourf( B, M, gt.fdata( ( B, M ) ) - gtd.fdata( ( B, M ) ) )
    plt.colorbar( cnf, ax=ax[1, 0] )

    # Plot relative difference in between Reference and downsampled data
    ax[1, 1].set_title( "Rel.Difference: Ref-Down" )
    cnf = ax[1, 1].contourf( B, M, ( gt.fdata( ( B, M ) ) - gtd.fdata( ( B, M ) ) ) / np.abs( gtd.fdata( ( B, M ) ) ) )
    plt.colorbar( cnf, ax=ax[1, 1] )


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


def dGdtt_Magee( beta: np.ndarray, mu: float, mu_coeffs: np.ndarray, alpha=0.0 ) -> np.ndarray:
    # Check if input array is a Ndarray
    isnda = isinstance( beta, np.ndarray )

    # Extract coefficients
    A = 0.0
    B = mu_coeffs[0] # 5.0
    C = mu_coeffs[1] # 1.2

    # Calculate leading term
    lt      = np.pi / 16 / np.sqrt( 2 )
    expt    = np.exp( -beta**2.0 * ( A * mu**3.0 + B * mu**2.0 + C * mu ) / 4.0 )
    x       = beta**2.0 / 8.0 + alpha

    eps     = 1e-6
    if isnda:
        pos     = beta < eps
        x[pos]  = 0.0
    else:
        x       = 0.0 if beta < eps else x

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
    
    if isnda:
        pos = np.isnan( y0 )
        if pos.sum( ) > 0:
            num_pos = np.where( pos )[0][-1] + 1
            y0[:num_pos] = y0[num_pos]

    else:
        y0 = 0.0 if np.isnan( y0 ) else y0
    
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
    if isnda:
        pos     = np.isnan( yt )
        yt[pos] = 0.0
        dt      = 5e-2
        pos     = beta < dt
        num_pos = np.where( pos )[0][-1]
        yt[pos] = yt[num_pos] / dt * beta[pos]
    T2      = yt * expt

    return lt * ( T1 + T2 )


def load_model( fopath: str, name: str, downsample=0 ) -> IntegralModel:
    # Load data
    with h5py.File( os.path.join( fopath, name + ".h5" ) ) as fid:
        mu      = fid[ "mu" ][:]
        beta    = fid[ "beta" ][:]
        data    = fid[ "fcn" ][:]

    if downsample > 0:
        # Downsample mu
        f       = int( downsample + 1 )
        N       = mu.shape[0] // f
        nposm   = np.linspace( 0, N-1, N, dtype=int ) * f
        mur     = mu[ nposm ]

        # Downsample beta
        f       = int( downsample + 1 )
        N       = beta.shape[0] // f
        nposb   = np.linspace( 0, N-1, N, dtype=int ) * f
        betar   = beta[ nposb ]

        # Downsample data
        datar   = data[ nposb, : ]
        # datar   = datar[ :, nposm ]

        # Create Integral model
        imodel = IntegralModel( betar, mu, datar, name )

    else:
        imodel = IntegralModel( beta, mu, data, name )

    return imodel


def main( fopath: str ) -> None:
    # Load models
    gt      = load_model( fopath, "Gt" )
    gtx     = load_model( fopath, "Gtx" )
    gtxx    = load_model( fopath, "Gtxx" )
    gtt     = load_model( fopath, "Gtt" )
    gttx    = load_model( fopath, "Gttx" )
    gttxx   = load_model( fopath, "Gttxx" )

    gtd     = load_model( fopath, "Gt", downsample=2 )

    show_differences_dm( gt, gtx, "Gtx" )
    show_differences_dm( gtx, gtxx, "Gtxx" )
    show_differences_db( gt, gtt, "Gtt" )
    show_differences_dm( gtt, gttx, "Gttx" )
    show_differences_dm( gttx, gttxx, "Gttxx" )
    compare_discretization( gt, gtd, "Gt" )

    fig, ax = plt.subplots( 3, 1, sharex=True )
    ax[0].plot( gt.beta, norm_data( gt.data[:, 0] ), label="Gt" )
    ax[0].plot( gtx.beta, norm_data( gtx.data[:, 0] ), label="Gtx" )
    ax[0].plot( gtxx.beta, norm_data( gtxx.data[:, 0] ), label="Gtxx" )
    ax[0].plot( gtt.beta, norm_data( gtt.data[:, 0] ), label="Gtt" )
    ax[0].plot( gttx.beta, norm_data( gttx.data[:, 0] ), label="Gttx" )
    ax[0].plot( gttxx.beta, norm_data( gttxx.data[:, 0] ), label="Gttxx" )
    ax[0].set_xlabel( "Beta [s]" )
    ax[0].legend( )

    wrap_fcn = wrap_data( 
                            gt.data[:, 0],
                            gtx.data[:, 0],
                            gtxx.data[:, 0],
                            gtt.data[:, 0],
                            gttx.data[:, 0],
                            gttxx.data[:, 0]
                        )
    ax[1].plot( gt.beta, wrap_fcn )

    wrap_fcn = wrap_data( 
                            norm_data( gt.data[:, 0] ),
                            norm_data( gtx.data[:, 0] ),
                            norm_data( gtxx.data[:, 0] ),
                            norm_data( gtt.data[:, 0] ),
                            norm_data( gttx.data[:, 0] ),
                            norm_data( gttxx.data[:, 0] )
                        )
    ax[2].plot( gt.beta, wrap_fcn )

    plt.show( )


def main2( fopath: str ) -> None:
    # Load models
    gt      = load_model( fopath, "Gt" )
    gtx     = load_model( fopath, "Gtx" )
    gtxx    = load_model( fopath, "Gtxx" )
    gtt     = load_model( fopath, "Gtt" )
    gttx    = load_model( fopath, "Gttx" )
    gttxx   = load_model( fopath, "Gttxx" )

    gtd     = load_model( fopath, "Gt", downsample=2 )

    # Generate approximate functions
    dgdt    = np.zeros_like( gt.data )
    dgdtt   = np.zeros_like( gtt.data )
    dgdtto  = np.zeros_like( gtt.data )
    for i in range( gtt.mu.shape[0] ):
        dgdt[:, i]  = dGdt_Magee( gt.beta, gt.mu[i] )
        dgdtt[:, i] = dGdtt_Magee( gtt.beta, gtt.mu[i], np.array( [ 0.0, 1.0 ] ) )
        dgdtto[:, i] = dGdtt_Magee( gtt.beta, gtt.mu[i], np.array( [ 0.0, 1.3 ] ) )

    # Derivate function dgdt to compare with dgdtt
    dgdttn  = ( ( dgdt[2:, :] - dgdt[:-2, :] ).T / ( gt.beta[2:] - gt.beta[:-2] ) ).T
    beta2   = gt.beta[1:-1]

    nb      = 1200

    
    A   = np.zeros_like( gtt.beta )
    B   = np.zeros_like( gtt.beta )
    res = np.zeros_like( gtt.beta )
    st  = 5
    for i in range( st, gtt.beta.shape[0] ):
        optifcn = lambda x: dGdtt_Magee( gtt.beta[i], gtt.mu, x, alpha=0.0 )
        resfcn  = lambda x: np.max( np.abs( gtt.data[i, :] - optifcn( x ) ) )
        optmr   = sp.optimize.minimize( resfcn, np.array( [ 0.0, 1.0 ] ) )
        print( "Step: ", i, " - ", optmr.x, resfcn( optmr.x ) )
        A[i]    = optmr.x[0]
        B[i]    = optmr.x[1]
        res[i]  = resfcn( optmr.x )

    A[:st] = A[st]
    B[:st] = B[st]

    intfcn_a = sp.interpolate.CubicSpline( gtt.beta, A )
    intfcn_b = sp.interpolate.CubicSpline( gtt.beta, B )

    for i in range( st, gtt.beta.shape[0] ):
        dgdtto[i, :] = dGdtt_Magee( gtt.beta[i], gtt.mu, np.array( [ intfcn_a( ( gtt.beta[i] ) ), intfcn_b( ( gtt.beta[i] ) ) ] ) )
        if i % 100:
            print( "ROF: ", i )

    fig, ax = plt.subplots( 3, 1, sharex=True )

    ax[0].plot( gtt.beta, A )
    ax[0].set_xlabel( "Beta" )
    ax[0].set_ylabel( "A" )

    ax[1].plot( gtt.beta, B )
    ax[1].set_xlabel( "Beta" )
    ax[1].set_ylabel( "B" )

    ax[2].plot( gtt.beta, res )
    ax[2].set_xlabel( "Beta" )
    ax[2].set_ylabel( "res" )

    plt.show( )

    # N  = 100
    # xa = np.linspace( 0.0, 25, N )
    # xb = np.linspace( 0.0, 5.0, N )
    # XA, XB = np.meshgrid( xa, xb, indexing="ij" )
    # RES = np.zeros( ( N, N ) )
    # for i in range( N ):
    #     for j in range( N ):
    #         RES[i, j] = np.sum( np.abs( gtt.data[nb, :] - optifcn( np.array( [ XA[i, j], XB[i, j] ] ) ) ) )

    # cnf = plt.contourf( XA, XB, np.log10( np.abs( RES ) ) )
    # plt.colorbar( )
    # plt.show( )

    # res = sp.optimize.minimize( resfcn, np.array( [ 0.0, 5.0, 1.2 ] ) )

    # print( "Coeffs:     ", res.x )
    # print( "Status:     ", res.status )
    # print( "Message:    ", res.message )

    # plt.plot( gtt.mu, optifcn( np.array( [ 5.0, 1.2 ] ) ) )
    # plt.show( )

    xo = np.array( [ 0.0, 1.3 ] )

    fig, ax = plt.subplots( 2, 1, sharex=True )
    fig.suptitle( f"Beta: {gt.beta[nb]}" )

    ax[0].plot( gt.mu, dgdt[nb, :], label="dgdt" )
    ax[0].plot( gtt.mu, dgdtt[nb, :], label="dgdtt" )
    ax[0].plot( gtt.mu, dgdttn[nb-1, :], label="dgdttn" )
    ax[0].plot( gtt.mu, optifcn( xo ), label="optifcn" )
    ax[0].plot( gt.mu, gt.data[nb, :], label="Gt" )
    ax[0].plot( gtt.mu, gtt.data[nb, :], label="Gtt" )
    ax[0].set_xscale( "log" )
    ax[0].legend( )

    ax[1].plot( gt.mu, gt.data[nb, :] - dgdt[nb, ], label="Gt-dgdt" )
    ax[1].plot( gtt.mu, gtt.data[nb, :] - dgdtt[nb, ], label="Gtt-dgdtt" )
    ax[1].plot( gtt.mu, gtt.data[nb, :] - optifcn( xo ), label="Gtt-optifcn" )
    ax[1].set_xscale( "log" )
    ax[1].legend( )

    fig, ax = plt.subplots( 2, 3, sharex=True, sharey=True )

    ax[0, 0].set_title( "Ref.Function" )
    cnf = ax[0, 0].contourf( gtt.B, gtt.M, gtt.data )
    plt.colorbar( cnf, ax=ax[0, 0] )

    ax[0, 1].set_title( "MageeStandard" )
    cnf = ax[0, 1].contourf( gtt.B, gtt.M, dgdtt )
    plt.colorbar( cnf, ax=ax[0, 1] )

    ax[0, 2].set_title( "MageeOpti" )
    cnf = ax[0, 2].contourf( gtt.B, gtt.M, dgdtto )
    plt.colorbar( cnf, ax=ax[0, 2] )

    ax[1, 1].set_title( "MageeStandard - Diff" )
    cnf = ax[1, 1].contourf( gtt.B, gtt.M, gtt.data - dgdtt, vmin=-1.0, vmax=1.0, levels=50 )
    plt.colorbar( cnf, ax=ax[1, 1] )

    ax[1, 2].set_title( "MageeOpti - Diff" )
    cnf = ax[1, 2].contourf( gtt.B, gtt.M, gtt.data - dgdtto, vmin=-1.0, vmax=1.0, levels=50 )
    plt.colorbar( cnf, ax=ax[1, 2] )

    plt.show( )

    # plt.plot( gtt.mu, gtt.data[nb, :] )
    # plt.plot( gtt.mu, dgdtt[nb, :] )
    # plt.show( )


def norm_data( fcn: np.ndarray ) -> np.ndarray:
    return fcn / np.max( np.abs( fcn ) )


def show_differences_db( gt: IntegralModel, gtd: IntegralModel, title: str, show=False ) -> None:
    # Create figure and canvas layout
    fig, ax = plt.subplots( 2, 2 )
    fig.suptitle( title )
    
    cnf = ax[0, 0].contourf( gtd.B, gtd.M, gtd.data )
    plt.colorbar( cnf, ax=ax[0, 0] )

    cnf = ax[0, 1].contourf( gt.B_db, gt.M_db, gt.data_db )
    plt.colorbar( cnf, ax=ax[0, 1] )

    cnf = ax[1, 0].contourf( gt.B_db, gt.M_db, gtd.data[1:-1, :]-gt.data_db )
    plt.colorbar( cnf, ax=ax[1, 0] )

    # Plot lines along beta direction
    colors  = plt.cm.tab10.colors
    N       = 10
    pos     = np.linspace( 0, gt.beta.shape[0]-3, N, dtype=int )
    for i in range( N ):
        ax[1, 1].plot( gtd.mu, gtd.data[pos[i]+1, :], color=colors[i] )
        ax[1, 1].plot( gt.mu, gt.data_db[pos[i]+0, :], "--", color=colors[i] )

    if show:
        plt.show( )


def show_differences_dm( gt: IntegralModel, gtd: IntegralModel, title: str, show=False ) -> None:
    # Create figure and canvas layout
    fig, ax = plt.subplots( 2, 2 )
    fig.suptitle( title )
    
    cnf = ax[0, 0].contourf( gtd.B, gtd.M, gtd.data )
    plt.colorbar( cnf, ax=ax[0, 0] )

    cnf = ax[0, 1].contourf( gt.B_dm, gt.M_dm, gt.data_dm )
    plt.colorbar( cnf, ax=ax[0, 1] )

    cnf = ax[1, 0].contourf( gt.B_dm, gt.M_dm, gtd.data[:, 1:-1]-gt.data_dm )
    plt.colorbar( cnf, ax=ax[1, 0] )

    # Plot lines along beta direction
    colors  = plt.cm.tab10.colors
    N       = 10
    pos     = np.linspace( 0, gt.mu.shape[0]-3, N, dtype=int )
    for i in range( N ):
        ax[1, 1].plot( gtd.beta, gtd.data[:, pos[i]+1], color=colors[i] )
        ax[1, 1].plot( gt.beta, gt.data_dm[:, pos[i]+0], "--", color=colors[i] )

    if show:
        plt.show( )


def wrap_data( *args ) -> np.ndarray:
    res = np.zeros_like( args[0] )
    for i in range( len( args ) ):
        pos         = np.abs( args[i][:] ) > np.abs( res )
        res[pos]    = args[i][pos]

    return res


if __name__ == "__main__":
    folder_path = r"D:\sergio\0050_OASIS_SM\TimeDomain"
    # main( folder_path )
    main2( folder_path )