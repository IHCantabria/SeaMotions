
# Import general usage libraries
import os

# Import general usage scientific libraries
import numpy as np

# Import general usage scientific plotting libraries
import matplotlib.pyplot as plt
from bokeh import plotting as bkl

# Import local modules
from base_integrals_v2 import fxy, fxy_dx, L1, L1_dA, L1_dB, L2, L3, L3_dA, L3_dB, M1, M1_dA, M1_dB, M2, M3, M3_dA, M3_dB


def calculate_L1( fopath: str ) -> None:
    N       = 10
    Amat    = np.linspace( 0, 1, N ) 
    Bmat    = np.linspace( 0, 1, N )
    Hmat    = 10 ** np.linspace( -2.999565922520681, 0, N )
    data    = np.zeros( ( N**3, 3 ) )
    count   = 0

    for i in range( N ):
        for j in range( N ):
            for k in range( N ):
                data[count, 0]  = L1( Amat[j], Bmat[k], Hmat[i] )[0].real
                data[count, 1]  = L1_dA( Amat[j], Bmat[k], Hmat[i] )[0].real
                data[count, 2]  = L1_dB( Amat[j], Bmat[k], Hmat[i] )[0].real
                count           += 1

    fipath = os.path.join( fopath, "L1.dat" )
    write_3d_file( fipath, Amat, Bmat, Hmat, data )


def calculate_L2( fopath: str ) -> None:
    N       = 10
    Hmat    = 10 ** np.linspace( -2.999565922520681, 0, N )
    data    = np.zeros( ( N, 1 ) )
    count   = 0

    for i in range( N ):
        data[count, 0]  = L2( Hmat[i] )[0].real
        count           += 1

    fipath = os.path.join( fopath, "L2.dat" )
    write_1d_file( fipath, Hmat, data )


def calculate_L3( fopath: str ) -> None:
    N       = 10
    Amat    = np.linspace( 0, 1, N ) 
    Bmat    = np.linspace( 0, 1, N )
    Hmat    = 10 ** np.linspace( 0, 3, N )
    data    = np.zeros( ( N**3, 3 ) )
    count   = 0

    for i in range( N ):
        for j in range( N ):
            for k in range( N ):
                data[count, 0]  = L3( Amat[j], Bmat[k], Hmat[i] )[0].real
                data[count, 1]  = L3_dA( Amat[j], Bmat[k], Hmat[i] )[0].real
                data[count, 2]  = L3_dB( Amat[j], Bmat[k], Hmat[i] )[0].real
                count           += 1

    fipath = os.path.join( fopath, "L3.dat" )
    write_3d_file( fipath, Amat, Bmat, Hmat, data )


def calculate_M1( fopath: str ) -> None:
    N       = 10
    Amat    = np.linspace( 0, 1, N ) 
    Bmat    = np.linspace( 1, 2, N )
    Hmat    = 10 ** np.linspace( -2.999565922520681, 0, N )
    data    = np.zeros( ( N**3, 3 ) )
    count   = 0

    for i in range( N ):
        for j in range( N ):
            for k in range( N ):
                data[count, 0]  = M1( Amat[j], Bmat[k], Hmat[i] )[0].real
                data[count, 1]  = M1_dA( Amat[j], Bmat[k], Hmat[i] )[0].real
                data[count, 2]  = M1_dB( Amat[j], Bmat[k], Hmat[i] )[0].real
                count           += 1

    fipath = os.path.join( fopath, "M1.dat" )
    write_3d_file( fipath, Amat, Bmat, Hmat, data )


def calculate_M2( fopath: str ) -> None:
    N       = 10
    Hmat    = 10 ** np.linspace( -2.999565922520681, 0, N )
    data    = np.zeros( ( N, 1 ) )
    count   = 0

    for i in range( N ):
        data[count, 0]  = M2( Hmat[i] )[0].real
        count           += 1

    fipath = os.path.join( fopath, "M2.dat" )
    write_1d_file( fipath, Hmat, data )


def calculate_M3( fopath: str ) -> None:
    N       = 10
    Amat    = np.linspace( 0, 1, N ) 
    Bmat    = np.linspace( 1, 2, N )
    Hmat    = 10 ** np.linspace( 0, 3, N )
    data    = np.zeros( ( N**3, 3 ) )
    count   = 0

    for i in range( N ):
        for j in range( N ):
            for k in range( N ):
                data[count, 0]  = M3( Amat[j], Bmat[k], Hmat[i] )[0].real
                data[count, 1]  = M3_dA( Amat[j], Bmat[k], Hmat[i] )[0].real
                data[count, 2]  = M3_dB( Amat[j], Bmat[k], Hmat[i] )[0].real
                count           += 1

    fipath = os.path.join( fopath, "M3.dat" )
    write_3d_file( fipath, Amat, Bmat, Hmat, data )


def calculate_R11( fopath: str ) -> None:
    N       = 10
    X       = 10**np.linspace( -5.0, 5.0, N ) 
    Y       = 10**np.linspace( -5.0, 5.0, N )
    data    = np.zeros( ( N**2, 2 ) )
    count   = 0

    for i in range( N ):
        for j in range( N ):
            data[count, 0]  = fxy( X[i], Y[j], only_int=True )
            data[count, 1]  = fxy_dx( X[i], Y[j], only_int=True )
            count           += 1

    fipath = os.path.join( fopath, "R11.dat" )
    write_2d_file( fipath, X, Y, data )


def main( fopath: str ) -> None:
    calculate_L1( fopath )
    calculate_L2( fopath )
    calculate_L3( fopath )
    calculate_M1( fopath )
    calculate_M2( fopath )
    calculate_M3( fopath )
    calculate_R11( fopath )


def plot_HBA_fig( fopath: str, name: str, A: np.ndarray, B: np.ndarray, H: np.ndarray, data: np.ndarray ) -> None:
    # Set storage path
    fipath = os.path.join( fopath, name + "_HBA.html" )
    bkl.output_file( fipath )

    # Create figure
    

def plot_integrals(  ) -> None:
    plot_L1( )


def plot_L1( fopath: str ) -> None:
    NA      = 10
    NB      = 10
    NH      = 100
    Amat    = np.linspace( 0, 1, NA ) 
    Bmat    = np.linspace( 0, 1, NB )
    Hmat    = 10 ** np.linspace( -2.999565922520681, 0, NH )
    G       = np.zeros( ( NA, NB, NH ) )
    G_dA    = np.zeros( ( NA, NB, NH ) )
    G_dB    = np.zeros( ( NA, NB, NH ) )
    count   = 0

    for i in range( NH ):
        for j in range( NA ):
            for k in range( NB ):
                G[j, k, i]      = L1( Amat[j], Bmat[k], Hmat[i] )[0].real
                G_dA[j, k, i]   = L1_dA( Amat[j], Bmat[k], Hmat[i] )[0].real
                G_dB[j, k, i]   = L1_dB( Amat[j], Bmat[k], Hmat[i] )[0].real
                count           += 1

    plot_HBA_fig( fopath, "L1", Amat, Bmat, H_mat, G )


def plot_L3( fopath: str ) -> None:
    name    = "L3"
    NA      = 3
    NB      = 3
    NH      = 100
    # Amat    = np.linspace( 0, 1, NA ) 
    # Bmat    = np.linspace( 0, 1, NB )
    # Hmat    = 10 ** np.linspace( 0, 3, NH )
    # Amat    = np.linspace( 0.8125, 0.875, NA )
    # Bmat    = np.linspace( 0.75, 0.8125, NB )
    # Hmat    = np.linspace( 10.75, 13.34, NH )
    # Amat    = np.linspace( 0.75, 0.875, NA )
    # Bmat    = np.linspace( 0.75, 0.875, NB )
    # Hmat    = np.linspace( 10.75, 13.34, NH )
    # Amat    = np.linspace( 0.75, 1.0, NA )
    # Bmat    = np.linspace( 0.75, 1.0, NB )
    # Hmat    = np.linspace( 0.3125e2, 1e3 , NH )
    Amat    = np.linspace( 0.25, 0.375, NA )
    Bmat    = np.linspace( 0.375, 0.5, NB )
    Hmat    = np.linspace( 2.371, 5.623 , NH )
    G       = np.zeros( ( NA, NB, NH ) )
    G_dA    = np.zeros( ( NA, NB, NH ) )
    G_dB    = np.zeros( ( NA, NB, NH ) )
    count   = 0

    for i in range( NH ):
        for j in range( NA ):
            for k in range( NB ):
                G[j, k, i]      = L3( Amat[j], Bmat[k], Hmat[i] )[0].real
                G_dA[j, k, i]   = L3_dA( Amat[j], Bmat[k], Hmat[i] )[0].real
                G_dB[j, k, i]   = L3_dB( Amat[j], Bmat[k], Hmat[i] )[0].real
                count           += 1

    Hmat = np.log10( Hmat )
    Hmat2 = ( Hmat[1:] + Hmat[:-1] ) / 2.0
    Hmat4 = ( Hmat2[1:] + Hmat2[:-1] ) / 2.0

    for i in range( NA ):
        fig = plt.figure( )
        fig.suptitle( f"A: {Amat[i]}" )

        ax0  = fig.add_subplot( 311 )
        ax0.set_xlabel( "H" )
        ax0.set_ylabel( f"{name}" )
        ax1  = fig.add_subplot( 312 )
        ax1.set_xlabel( "H" )
        ax1.set_ylabel( f"{name}_dA" )
        ax2  = fig.add_subplot( 313 )
        ax2.set_xlabel( "H" )
        ax2.set_ylabel( f"{name}_dB" )

        for j in range( NB ):
            ax0.plot( Hmat, G[i, j, :], label=f"B: {Bmat[j]:0.2f}" )
            # ax1.plot( Hmat, G_dA[i, j, :], label=f"B: {Bmat[j]:0.2f}" )
            # ax2.plot( Hmat, G_dB[i, j, :], label=f"B: {Bmat[j]:0.2f}" )
            ax1.plot( Hmat2, np.diff( G[i, j, :] ), label=f"B: {Bmat[j]:0.2f}" )
            ax2.plot( Hmat4, np.diff( np.diff( G[i, j, :] ) ), label=f"B: {Bmat[j]:0.2f}" )
        ax0.legend( )
        ax1.legend( )
        ax2.legend( )

    for i in range( NB ):
        fig = plt.figure( )
        fig.suptitle( f"B: {Bmat[i]}" )

        ax0  = fig.add_subplot( 311 )
        ax0.set_xlabel( "H" )
        ax0.set_ylabel( f"{name}" )
        ax1  = fig.add_subplot( 312 )
        ax1.set_xlabel( "H" )
        ax1.set_ylabel( f"{name}_dA" )
        ax2  = fig.add_subplot( 313 )
        ax2.set_xlabel( "H" )
        ax2.set_ylabel( f"{name}_dB" )

        for j in range( NA ):
            ax0.plot( Hmat, G[j, i, :], label=f"A: {Amat[j]:0.2f}" )
            # ax1.plot( Hmat, G_dA[j, i, :], label=f"A: {Amat[j]:0.2f}" )
            # ax2.plot( Hmat, G_dB[j, i, :], label=f"A: {Amat[j]:0.2f}" )
            ax1.plot( Hmat2, np.diff( G[j, i, :] ), label=f"A: {Amat[j]:0.2f}" )
            ax2.plot( Hmat4, np.diff( np.diff( G[j, i, :] ) ), label=f"A: {Amat[j]:0.2f}" )
        ax0.legend( )
        ax1.legend( )
        ax2.legend( )

    plt.show( )


def write_vector(file, vec, name):
    file.write(f'Num.{name}: {len(vec)}\n')
    format_string = ' '.join(['{:.15f}'] * len(vec)) + '\n'
    file.write(format_string.format(*vec))

def write_1d_file(fipath, H, fcn):
    with open(fipath, 'w') as file:
        write_vector(file, H, 'H')
        file.write('G\n')
        fcn = np.asarray(fcn).flatten()
        for value in fcn:
            file.write(f'{value:.15f}\n')

def write_2d_file(fipath, A, B, fcn):
    with open(fipath, 'w') as file:
        write_vector(file, A, 'A')
        write_vector(file, B, 'B')
        file.write('G    G_dX\n')
        fcn = np.asarray(fcn)
        for row in fcn:
            file.write('{:.15f} {:.15f}\n'.format(*row))

def write_3d_file(fipath, A, B, H, fcn):
    with open(fipath, 'w') as file:
        write_vector(file, A, 'A')
        write_vector(file, B, 'B')
        write_vector(file, H, 'H')
        file.write('G    G_dA    G_dB\n')
        fcn = np.asarray(fcn)
        for row in fcn:
            file.write('{:.15f} {:.15f} {:.15f}\n'.format(*row))


if __name__ == "__main__":
    this_path   = os.path.dirname( os.path.abspath( __file__ ) )
    folder_path = os.path.join( this_path, "..", "aux_data", "2_check_database_files", "0_frequency_domain" )

    main( folder_path )
