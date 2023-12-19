
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
                                                cuscomplex*     raos,
                                                cuscomplex*     rel_wave_elevation
                                            )
{
    // Calculate wave elevation
    calculate_wave_elevation_lin(
                                    potpanel_total,
                                    input->heads_np * pot_gp->field_points_np,
                                    ang_freq,
                                    input->grav_acc,
                                    rel_wave_elevation
                                );
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
                rel_wave_elevation[we_index]    -= point_disp[2];
            }
        }
    }
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