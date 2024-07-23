
# Import general usage libraries
import h5py

# Import general usage scientific libraries
from numpy import exp, pi, sqrt, zeros
from numpy import abs as np_abs
from numpy.linalg import solve as np_solve

# Import general usage plotting libraries
import matplotlib.pyplot as plt


def main(
            sm_fipath: str,
            aq_fipath: str
        ) -> None:
    # Get data from results file
    with h5py.File( sm_fipath, "r" ) as fid:
        # Get data from seamotions
        sm_added_mass   = fid[ "added_mass" ][:]
        sm_damping      = fid[ "damping_rad" ][:]
        sm_freqs        = fid[ "frequencies" ][:]
        sm_hydrostatic  = fid[ "hydstiffness" ][:]
        sm_mass         = fid[ "mass" ][:]
        sm_rao_mag      = fid[ "rao_mag" ][:]
        sm_wex_mag      = fid[ "wave_exciting_mag" ][:]
        sm_wex_pha      = fid[ "wave_exciting_pha" ][:]

    with h5py.File( aq_fipath, "r" ) as fid:
        # Get data from Aqwa
        aq_freqs        = fid[ "body_0" ][ "frequencies" ][:]
        aq_added_mass   = fid[ "body_0" ][ "added_mass" ][:]
        aq_damping      = fid[ "body_0" ][ "damping_radiation" ][:]
        aq_hydrostatic  = fid[ "body_0" ][ "hydstiffness" ][:]
        aq_mass         = fid[ "body_0" ][ "mass" ][:]
        aq_rao_mag      = fid[ "body_0" ][ "rao_movements_mag" ][:]
        aq_wex_mag      = fid[ "body_0" ][ "wave_exciting_mag" ][:]
        aq_wex_pha      = fid[ "body_0" ][ "wave_exciting_mag" ][:]

    # print( "AQ - MASS" )
    # print( aq_mass )
    # print( "SM - MASS" )
    # print( sm_mass )

    aq_rao  = zeros( ( aq_freqs.shape[0], 6 ) )
    sm_rao  = zeros( ( aq_freqs.shape[0], 6 ) )
    hd_pos  = 0
    for fn in range( aq_freqs.shape[0] ):
        w = 2 * pi * aq_freqs[fn]
        # print( "Wex: ", aq_wex_mag[1, fn, 0] )
        # print( "A: ", aq_added_mass[0, fn, 0, 0] )
        # print( "B: ", aq_damping[0, fn, 0, 0] )
        # rao = aq_wex_mag[1, fn, 0] / sqrt( ( w**2.0 * ( aq_mass[0, 0] + aq_added_mass[0, fn, 0, 0] ) )**2.0 + ( w * aq_damping[0, fn, 0, 0] )**2.0 )
        # print( "rao_aq: ", aq_rao_mag[1, fn, 0] )
        # print( "rao: ", rao )

        aq_hydromech_mat = -w**2.0*(aq_mass + aq_added_mass[0, fn, :, :]) - w*1j*aq_damping[0, fn, :, :] + aq_hydrostatic
        aq_wex_vec = aq_wex_mag[hd_pos, fn, :]*exp(-1j*aq_wex_pha[hd_pos, fn, :])
        aq_rao[fn, :] = np_abs( np_solve(aq_hydromech_mat, aq_wex_vec) )

        # print("AQ - SYSMAT")
        # print(np_abs(aq_hydromech_mat))
        # print( "AQ - RHS" )
        # print(np_abs(aq_wex_vec))
        # print( "AQ - RAO" )
        # print(aq_rao)

        sm_hydromech_mat = -w**2.0*(sm_mass + sm_added_mass[0, 0, fn, :, :]) - w*1j*sm_damping[0, 0, fn, :, :] + sm_hydrostatic
        sm_wex_vec = sm_wex_mag[0, hd_pos, fn, :]*exp(-1j*sm_wex_pha[0, hd_pos, fn, :])
        # print(sm_hydromech_mat.shape)
        # print(sm_wex_vec.shape)
        sm_rao[fn, :] = np_abs( np_solve(sm_hydromech_mat[0, :, :], sm_wex_vec) )


        # print( "SM - SYSMAT" )
        # print(np_abs(sm_hydromech_mat))
        # print( "SM - RHS" )
        # print(np_abs(sm_wex_vec))
        # print( "SM - RAO" )
        # print(sm_rao)

        # print( "AQ vs SM - SYSMAT" )
        # print( np_abs( aq_hydromech_mat - sm_hydromech_mat[0, :, :] ) )

    # Plot results
    fig = plt.figure( )
    ax0 = fig.add_subplot( 221 )
    ax1 = fig.add_subplot( 222 )
    ax2 = fig.add_subplot( 223 )
    ax3 = fig.add_subplot( 224 )

    ax0.set_title( "Added Mass" )
    ax1.set_title( "Damping" )
    ax2.set_title( "Wex Total" )
    ax3.set_title( "RAO" )

    DOF = 2

    ax0.plot( sm_freqs, sm_added_mass[0, 0, :, DOF, DOF], "-o", label="SM" )
    ax0.plot( aq_freqs, aq_added_mass[0, :, DOF, DOF], "-o", label="AQ" )
    ax0.legend(  )

    ax1.plot( sm_freqs, sm_damping[0, 0, :, DOF, DOF], "-o", label="SM" )
    ax1.plot( aq_freqs, aq_damping[0, :, DOF, DOF], "-o", label="AQ" )
    ax1.legend(  )

    ax2.plot( sm_freqs, sm_wex_mag[0, 0, :, DOF], "-o", label="SM" )
    ax2.plot( aq_freqs, aq_wex_mag[1, :, DOF], "-o", label="AQ" )
    ax2.legend(  )

    ax3.plot( sm_freqs, sm_rao_mag[0, 0, :, DOF], "-o", label="SM" )
    ax3.plot( aq_freqs, aq_rao_mag[1, :, DOF], "-o", label="AQ" )
    # ax3.plot( sm_freqs, sm_rao[:, DOF], "-o", label="SMpy" )
    # ax3.plot( aq_freqs, aq_rao[:, DOF], "-o", label="AQpy" )
    ax3.legend(  )

    plt.show( )


def main_2(
            sm_fipath: str,
            aq_fipath: str
        ) -> None:
    # Get data from results file
    with h5py.File( sm_fipath, "r" ) as fid:
        # Get data from seamotions
        sm_added_mass   = fid[ "added_mass" ][:]
        sm_damping      = fid[ "damping_rad" ][:]
        sm_freqs        = fid[ "frequencies" ][:]
        sm_hydrostatic  = fid[ "hydstiffness" ][:]
        sm_mass         = fid[ "mass" ][:]
        sm_rao_mag      = fid[ "rao_mag" ][:]
        sm_wex_mag      = fid[ "wave_exciting_mag" ][:]
        sm_wex_pha      = fid[ "wave_exciting_pha" ][:]

    with h5py.File( aq_fipath, "r" ) as fid:
        # Get data from Aqwa
        aq_freqs        = fid[ "frequencies" ][:]
        aq_added_mass   = fid[ "added_mass" ][:]
        aq_damping      = fid[ "damping_rad" ][:]
        aq_hydrostatic  = fid[ "hydstiffness" ][:]
        aq_mass         = fid[ "mass" ][:]
        aq_rao_mag      = fid[ "rao_mag" ][:]
        aq_wex_mag      = fid[ "wave_exciting_mag" ][:]
        aq_wex_pha      = fid[ "wave_exciting_mag" ][:]

    # Plot results
    fig = plt.figure( )
    ax0 = fig.add_subplot( 221 )
    ax1 = fig.add_subplot( 222 )
    ax2 = fig.add_subplot( 223 )
    ax3 = fig.add_subplot( 224 )

    ax0.set_title( "Added Mass" )
    ax1.set_title( "Damping" )
    ax2.set_title( "Wex Total" )
    ax3.set_title( "RAO" )

    DOF = 2

    ax0.plot( sm_freqs, sm_added_mass[0, 0, :, DOF, DOF], "-o", label="SM" )
    ax0.plot( aq_freqs, aq_added_mass[0, 0, :, DOF, DOF], "-o", label="AQ" )
    ax0.legend(  )

    ax1.plot( sm_freqs, sm_damping[0, 0, :, DOF, DOF], "-o", label="SM" )
    ax1.plot( aq_freqs, aq_damping[0, 0, :, DOF, DOF], "-o", label="AQ" )
    ax1.legend(  )

    ax2.plot( sm_freqs, sm_wex_mag[0, 0, :, DOF], "-o", label="SM" )
    ax2.plot( aq_freqs, aq_wex_mag[0, 0, :, DOF], "-o", label="AQ" )
    ax2.legend(  )

    ax3.plot( sm_freqs, sm_rao_mag[0, 0, :, DOF], "-o", label="SM" )
    ax3.plot( aq_freqs, aq_rao_mag[0, 0, :, DOF], "-o", label="AQ" )
    ax3.legend(  )

    plt.show( )


if __name__ == "__main__":
    sm_file_path = r"E:\sergio\developments\SeaMotions\examples\_cube\results_noLID.hydb.h5"
    # sm_file_path = r"E:\sergio\developments\SeaMotions\examples\_cube\results_EL_500_Pr_1Em7.hydb.h5"
    aq_file_path = r"E:\sergio\developments\SeaMotions\examples\_cube\results_EL_500_Pr_1Em7.hydb.h5"
    # aq_file_path = r"E:\sergio\developments\SeaMotions\examples\_cube\ANALYSIS_noLID.hydb.h5"
    # main( sm_file_path, aq_file_path )
    main_2( sm_file_path, aq_file_path )