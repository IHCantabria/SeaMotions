
// Include local modules
#include "qtf.hpp"

#include "froude_krylov.hpp"
#include "../../math/euler_transforms.hpp"
#include "../../math/math_interface.hpp"
#include "../../math/math_tools.hpp"
#include "../../waves.hpp"


void    calculate_second_order_force(
                                        Input*          input,
                                        MpiConfig*      mpi_config,
                                        MeshGroup*      mesh_gp,
                                        cuscomplex*     mdrift_rel_we_i,
                                        cuscomplex*     mdrift_rel_we_j,
                                        cuscomplex*     raos_i,
                                        cuscomplex*     raos_j,
                                        cuscomplex*     vel_x_i,
                                        cuscomplex*     vel_y_i,
                                        cuscomplex*     vel_z_i,
                                        cuscomplex*     vel_x_j,
                                        cuscomplex*     vel_y_j,
                                        cuscomplex*     vel_z_j,
                                        cusfloat        ang_freq_i,
                                        cusfloat        ang_freq_j,
                                        cuscomplex*     qtf_values
                                    )
{
    // Define aux variables to be used along the function
    int         idx0    = 0;
    int         idx1    = 0;
    cuscomplex  int_mod = 0.0;
    PanelGeom*  panel_k = nullptr;

    // Calculate second order force due to relative wave
    // elevation at the WL
    for ( int i=0; i<input->heads_np; i++ )
    {
        for ( int j=0; j<mesh_gp->meshes_np; j++ )
        {
            idx0 = i * ( input->dofs_np * input->bodies_np ) + j * input->bodies_np;
            for ( int k=mesh_gp->panels_wl_cnp[j]; k<mesh_gp->panels_wl_cnp[j+1]; j++ )
            {
                panel_k             = mesh_gp->panels_wl[k];
                idx1                = i * mesh_gp->panels_wl_tnp + k;
                int_mod             = mdrift_rel_we_i[idx1] * mdrift_rel_we_j[idx1] * panel_k->len_wl;

                for ( int r=0; r<input->dofs_np; r++ )
                {
                    qtf_values[idx0+r] += - 0.5 * input->grav_acc * input->water_density * int_mod * panel_k->normal_vec[r];
                }
            }
        }
    }

    // Calculate second order force due to the bernouilly contribution
    for ( int i=0; i<input->heads_np; i++ )
    {
        for ( int j=0; j<mesh_gp->meshes_np; j++ )
        {
            idx0 = i * ( input->dofs_np * input->bodies_np ) + j * input->bodies_np;
            for ( int k=mesh_gp->panels_cnp[j]; k<mesh_gp->panels_cnp[j+1]; j++ )
            {
                panel_k             = mesh_gp->panels[k];
                idx1                = i * mesh_gp->panels_tnp + k;
                int_mod             = (
                                            vel_x_i[idx1] * vel_x_j[idx1]
                                            +
                                            vel_y_i[idx1] * vel_y_j[idx1]
                                            +
                                            vel_z_i[idx1] * vel_z_j[idx1]
                                        ) * panel_k->area;

                for ( int r=0; r<input->dofs_np; r++ )
                {
                    qtf_values[idx0+r] += input->water_density * int_mod * panel_k->normal_vec[r];
                }
            }
        }
    }

    // Calculate second order force due to acceleration term
    cusfloat    cog_to_fp[3];   clear_vector( 3, cog_to_fp );
    cusfloat    point_disp[3];  clear_vector( 3, point_disp );
    cusfloat    rao_trans[3];   clear_vector( 3, rao_trans );
    cusfloat    rao_rot[3];     clear_vector( 3, rao_rot );

    for ( int i=0; i<input->heads_np; i++ )
    {
        for ( int j=0; j<mesh_gp->meshes_np; j++ )
        {
            idx0 = i * ( input->dofs_np * input->bodies_np ) + j * input->bodies_np;

            rao_trans[0]    = std::abs( raos_i[idx0] );
            rao_trans[1]    = std::abs( raos_i[idx0+1] );
            rao_trans[2]    = std::abs( raos_i[idx0+2] );

            rao_rot[0]      = std::abs( raos_i[idx0+3] );
            rao_rot[1]      = std::abs( raos_i[idx0+4] );
            rao_rot[2]      = std::abs( raos_i[idx0+5] );

            for ( int k=mesh_gp->panels_cnp[j]; k<mesh_gp->panels_cnp[j+1]; j++ )
            {
                // Get handle to the current panel for calculations
                panel_k             = mesh_gp->panels[k];

                // Define vector from cog to field point
                sv_sub( 3, panel_k->center, panel_k->body_cog, cog_to_fp );

                // Calculate first order displacement of the panel centre
                euler_local_to_global_disp( 
                                                rao_trans,
                                                rao_rot,
                                                cog_to_fp,
                                                point_disp
                                            );

                // Calculate point displacement
                idx1                = i * mesh_gp->panels_tnp + k;
                int_mod             = (
                                            point_disp[0] * vel_x_j[idx1]
                                            +
                                            point_disp[1] * vel_y_j[idx1]
                                            +
                                            point_disp[2] * vel_z_j[idx1]
                                        ) * panel_k->area * cuscomplex( 0.0, ang_freq_j );

                for ( int r=0; r<input->dofs_np; r++ )
                {
                    qtf_values[idx0+r] += input->water_density * int_mod * panel_k->normal_vec[r];
                }
            }
        }
    }

    // Calculate second order force due to momentum
    cusfloat    ang_2                           = pow2s( ang_freq_j );
    cuscomplex  hydro_force[input->dofs_np];    clear_vector( input->dofs_np, hydro_force );
    cuscomplex  mom_i[3];                       clear_vector( 3, mom_i );

    for ( int i=0; i<input->heads_np; i++ )
    {
        for ( int j=0; j<mesh_gp->meshes_np; j++ )
        {
            // Define chunk index
            idx0 = i * ( input->dofs_np * input->bodies_np ) + j * input->bodies_np;

            // Calculate total hydrodynamic force
            hydro_force[0] = - input->bodies[j]->mass * raos_j[idx0] * ang_2;
            hydro_force[1] = - input->bodies[j]->mass * raos_j[idx0+1] * ang_2;
            hydro_force[2] = - input->bodies[j]->mass * raos_j[idx0+2] * ang_2;
            hydro_force[3] = - (
                                    input->bodies[j]->inertia[0] * raos_j[idx0+3]
                                    +
                                    input->bodies[j]->inertia[1] * raos_j[idx0+4]
                                    +
                                    input->bodies[j]->inertia[2] * raos_j[idx0+5]
                                ) * ang_2;
            hydro_force[4] = - (
                                    input->bodies[j]->inertia[1] * raos_j[idx0+3]
                                    +
                                    input->bodies[j]->inertia[3] * raos_j[idx0+4]
                                    +
                                    input->bodies[j]->inertia[4] * raos_j[idx0+5]
                                ) * ang_2;
            hydro_force[5] = - (
                                    input->bodies[j]->inertia[2] * raos_j[idx0+3]
                                    +
                                    input->bodies[j]->inertia[4] * raos_j[idx0+4]
                                    +
                                    input->bodies[j]->inertia[5] * raos_j[idx0+5]
                                ) * ang_2;

            // Add moment due to translational forces
            cross(
                        &(raos_i[idx0+3]),
                        hydro_force,
                        mom_i
                );

            sv_add(
                        3,
                        &(qtf_values[idx0]),
                        mom_i,
                        &(qtf_values[idx0])
                    );

            // Add moment due to rotational force
            cross(
                        &(raos_i[idx0+3]),
                        &(hydro_force[3]),
                        mom_i
                );

            sv_add(
                        3,
                        &(qtf_values[idx0+3]),
                        mom_i,
                        &(qtf_values[idx0+3])
                    );
            
        }
    }

    // Calculate second order potential contribution using
    // the approximate formula by Pinkster
    cusfloat    ai                              = input->wave_amplitude;
    cusfloat    aj                              = input->wave_amplitude;
    cusfloat    aij                             = 0.0;
    cusfloat    bij                             = 0.0;
    cusfloat    cij                             = 0.0;
    cusfloat    dk                              = 0.0;
    cusfloat    dkwij                           = 0.0;
    cusfloat    fij                             = 0.0;
    cuscomplex  froude_krylov[input->dofs_np];  clear_vector( input->dofs_np, froude_krylov );
    cusfloat    g                               = input->grav_acc;
    cusfloat    ki                              = 0.0;
    cusfloat    kj                              = 0.0;
    cusfloat    h                               = input->water_depth;
    cusfloat    wij                             = ang_freq_i - ang_freq_j;

    if ( std::abs( wij ) > 1e-12 )
    {
        for ( int i=0; i<input->heads_np; i++ )
        {
            for ( int j=0; j<mesh_gp->meshes_np; j++ )
            {
                 // Define chunk index
                idx0 = i * ( input->dofs_np * input->bodies_np ) + j * input->bodies_np;

                // Calculate ith, jth wave numbers
                ki = w2k( ang_freq_i, input->water_depth, input->grav_acc );
                kj = w2k( ang_freq_j, input->water_depth, input->grav_acc );

                // Calculate wave numbers difference and the equivalent angular frequency
                dk  = ki - kj;
                dkwij = k2w( dk, input->water_depth, input->grav_acc );

                // Calculate equivalent Froude-Krylov force for frequency difference
                calculate_froude_krylov(
                                            input,
                                            mpi_config,
                                            mesh_gp,
                                            dkwij,
                                            froude_krylov
                                        );
                
                // Calculate correction coefficients
                cij = 2 * ki * kj * wij * ( 1 + std::tanh( ki * h ) * std::tanh( kj * h ) ) / ang_freq_i / ang_freq_j;
                bij = (
                            pow2s( ki ) / ang_freq_i / pow2s( std::cosh( ki * h ) )
                            -
                            pow2s( kj ) / ang_freq_j / pow2s( std::cosh( kj * h ) )
                        );
                aij = (
                            ( bij + cij )
                            /
                            ( pow2s( wij ) - dk * g * std::tanh( dk * h ) )
                        ) * 0.5 * g;
                fij = ai * aj * aij * wij / g;

                // Calcuate second order force contribution using Pinkster approximation
                for ( int k=0; k<input->dofs_np; k++ )
                {
                    qtf_values[idx0+k] = fij * froude_krylov[k];
                }
            }
        }
    }
}