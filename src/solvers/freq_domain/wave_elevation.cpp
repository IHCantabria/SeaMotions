
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
                                    pot_gp->field_points_np,
                                    ang_freq,
                                    input->grav_acc,
                                    rel_wave_elevation
                                );

    // Calculate relative wave elevation
    int         fp_np           = pot_gp->field_points_np;
    int         we_index        = 0;
    cusfloat    point_disp[3]   = { 0.0, 0.0, 0.0 };
    cusfloat    rao_trans[3]    = { 0.0, 0.0, 0.0 };
    cusfloat    rao_rot[3]      = { 0.0, 0.0, 0.0 };
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
            rao_trans[0]    = std::abs( raos[rao_index] );
            rao_trans[1]    = std::abs( raos[rao_index+1] );
            rao_trans[2]    = std::abs( raos[rao_index+2] );

            rao_rot[0]      = std::abs( raos[rao_index+3] );
            rao_rot[1]      = std::abs( raos[rao_index+4] );
            rao_rot[2]      = std::abs( raos[rao_index+5] );

            for ( int k=pot_gp->field_points_cnp[j]; k<pot_gp->field_points_cnp[j+1]; k++ )
            {
                // Calculate movement of the kth point at the WL
                euler_local_to_global_disp( 
                                                rao_trans,
                                                rao_rot,
                                                &(pot_gp->cog_to_field_points[3*k]),
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
        wave_elevation[i] = pot_total[i].imag( ) * ang_freq / grav_acc;
    }

}