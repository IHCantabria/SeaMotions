
// Include local modules
#include "wave_elevation.hpp"

#include "../../containers/mpi_config.hpp"
#include "../../inout/input.hpp"
#include "../../math/euler_transforms.hpp"
#include "../../mesh/mesh_group.hpp"
#include "../../solvers/freq_domain/potential.hpp"


void    calculate_relative_wave_elevation_lin(
                                                Input*          input,
                                                MLGCmpx*        pot_gp,
                                                cuscomplex*     potpanel_total,
                                                cusfloat        ang_freq,
                                                int             ang_freq_num,
                                                cuscomplex*     raos,
                                                cuscomplex*     rel_wave_elevation
                                            )
{
    // Calculate wave elevation
    calculate_wave_elevation_lin(
                                    potpanel_total,
                                    pot_gp->field_points_np,
                                    ang_freq,
                                    input->grav_acc,
                                    rel_wave_elevation
                                );
    
    std::string     base_path( "E:/sergio/0050_OASIS_SM/_check_potentials/sm_freqs/wave_elevation/" );
    std::stringstream ss0; ss0 << "wave_elevation_ang_freq_" << ang_freq_num << ".dat";
    std::string     wave_elevation_fipath = base_path + ss0.str( );
    std::ofstream   wave_elevation_outf( wave_elevation_fipath );
    CHECK_FILE_UNIT_STATUS( wave_elevation_outf, wave_elevation_fipath );

    wave_elevation_outf << pot_gp->field_points_np << "  " << input->heads_np << std::endl;

    int idx = 0;
    for ( int i=0; i<input->heads_np; i++ )
    {
        for ( int j=0; j<pot_gp->field_points_np; j++ )
        {
            idx = i * pot_gp->field_points_np + j;
            wave_elevation_outf << rel_wave_elevation[idx].real( ) << "  ";
            wave_elevation_outf << rel_wave_elevation[idx].imag( ) << std::endl;
        }
    }

    wave_elevation_outf.close( );

    // std::string pot_wl_fipath( "E:/sergio/0050_OASIS_SM/wl_data.dat" );

    // std::ofstream of_pot_wl( pot_wl_fipath );

    // std::string space4( "    " );
    // for ( int i=0; i<pot_gp->sysmat_nrows; i++ )
    // {
    //     of_pot_wl << i+1 << space4;
    //     for ( int j=0; j<3; j++ )
    //     {
    //         of_pot_wl << pot_gp->field_points[3*i+j] << space4;
    //     }
    //     of_pot_wl << rel_wave_elevation[i].real( ) << space4 << rel_wave_elevation[i].imag( ) << space4;
    //     of_pot_wl << std::abs( rel_wave_elevation[i] ) << space4 << 57.3 * std::atan2( rel_wave_elevation[i].imag( ), rel_wave_elevation[i].real( ) ) << std::endl;

    // }

    // of_pot_wl.close( );

    // std::cout  << "PANELS POTENTIAL: " << std::endl;
    // for ( int i=0; i<pot_gp->field_points_np; i++ )
    // {
    //     std::cout << "FP[" << i << "]: " << potpanel_total[i] << " - ";
    //     std::cout << std::abs( potpanel_total[i] ) << " - " << 57.3 * std::atan2( potpanel_total[i].imag(), potpanel_total[i].real() ) << std::endl;
    // }
    // std::cout << std::endl;

    // std::cout << "WAVE ELEVATION: " << std::endl;
    // for ( int i=0; i<pot_gp->field_points_np; i++ )
    // {
    //     std::cout << "FP[" << i << "]: " << pot_gp->field_points[3*i] << " " << pot_gp->field_points[3*i+1] << " " << pot_gp->field_points[3*i+2];
    //     std::cout << " - Complex: " << rel_wave_elevation[i];
    //     std::cout << " - Mag/Pha: " << std::abs( rel_wave_elevation[i] ) << " - " << 57.3 * std::atan2( rel_wave_elevation[i].imag(), rel_wave_elevation[i].real() ) << std::endl;
    // }
    // std::cout << std::endl;

    // Calculate relative wave elevation
    int         fp_np           = pot_gp->field_points_np;
    int         we_index        = 0;
    cuscomplex  point_disp[3]   ;   clear_vector( 3, point_disp );
    cuscomplex  radius[3]       ;   clear_vector( 3, radius );
    cuscomplex  rao_trans[3]    ;   clear_vector( 3, rao_trans );
    cuscomplex  rao_rot[3]      ;   clear_vector( 3, rao_rot );
    int         rao_index       = 0;

    for ( int i=0; i<input->heads_np; i++ )
    {
        for ( int j=0; j<pot_gp->field_points_nb; j++ )
        {
            rao_index       =  (
                                    i * ( input->dofs_np * pot_gp->field_points_nb )
                                    +
                                    j * input->dofs_np
                                );

            for ( int k=pot_gp->field_points_cnp[j]; k<pot_gp->field_points_cnp[j+1]; k++ )
            {
                clear_vector( 3, point_disp );
                
                // Calculate movement of the kth point at the WL
                for ( int r=0; r<3; r++ )
                {
                    radius[r]       = cuscomplex( pot_gp->cog_to_field_points[3*k+r], 0.0 );
                    rao_trans[r]    = raos[rao_index+r] * input->wave_amplitude;
                    rao_rot[r]      = raos[rao_index+3+r] * input->wave_amplitude;
                }
                cross( 
                            rao_rot,
                            radius,
                            point_disp
                    );
                sv_add(
                            3,
                            point_disp,
                            rao_trans,
                            point_disp
                        );

                // Calculate relative water height by substracting the 
                // vertical motion at the required point to the wave elevation
                we_index                        = i * fp_np + k;
                std::cout << "WeIndex: " << we_index << " - " << pot_gp->field_points[3*we_index] << " " << pot_gp->field_points[3*we_index+1] << " " << pot_gp->field_points[3*we_index+2];
                std::cout << " - " << std::abs( rel_wave_elevation[we_index] ) << " - " << std::abs( point_disp[2] );
                std::cout << " - " << 57.3 * std::atan2( point_disp[2].imag( ), point_disp[2].real( ) ) << std::endl;
                rel_wave_elevation[we_index]    -= point_disp[2];
            }
        }
    }

    std::cout << "RELATIVE WAVE ELEVATION: " << std::endl;
    for ( int i=0; i<pot_gp->field_points_np; i++ )
    {
        std::cout << "FP[" << i << "]: " << pot_gp->field_points[3*i] << " " << pot_gp->field_points[3*i+1] << " " << pot_gp->field_points[3*i+2];
        std::cout << " - Complex: " << rel_wave_elevation[i];
        std::cout << " - Mag/Pha: " << std::abs( rel_wave_elevation[i] ) << " - " << 57.3 * std::atan2( rel_wave_elevation[i].imag(), rel_wave_elevation[i].real() ) << std::endl;
    }
    std::cout << std::endl;
 
}


void    calculate_wave_elevation_lin(
                                        cuscomplex*     pot_total,
                                        int             pot_total_np,
                                        cusfloat        ang_freq,
                                        cusfloat        grav_acc,
                                        cuscomplex*     wave_elevation
                                    )
{
    for ( int i=0; i<pot_total_np; i++ )
    {
        wave_elevation[i] =  cuscomplex( 0.0, 1.0 ) * pot_total[i] * ang_freq / grav_acc;
    }

}