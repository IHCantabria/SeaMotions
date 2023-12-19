
// Include local modules
#include "qtf.hpp"

#include "froude_krylov.hpp"
#include "../../math/gauss.hpp"
#include "../../math/euler_transforms.hpp"
#include "../../math/math_interface.hpp"
#include "../../math/math_tools.hpp"
#include "../../math/topology.hpp"
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
                                            cuscomplex*     qtf_values,
                                            cuscomplex*     qtf_wl,
                                            cuscomplex*     qtf_bern,
                                            cuscomplex*     qtf_acc,
                                            cuscomplex*     qtf_mom,
                                            MLGCmpx*        pot_gp,
                                            MLGCmpx*        vel_gp
                                    )
{
    // Define aux variables to be used along the function
    GaussPoints gp( input->gauss_order );
    int         idx0    = 0;
    int         idx1    = 0;
    cuscomplex  int_mod = 0.0;
    PanelGeom*  panel_k = nullptr;
    int         ngp     = input->gauss_order;
    int         ngpf_1d = input->gauss_np_factor_1d( );
    int         ngpf    = input->gauss_np_factor_2d( );

    // Clear QTF input vector to ensure that previous data will not be storaged
    // erroneously
    clear_vector( 
                    input->heads_np * mesh_gp->meshes_np * input->dofs_np,
                    qtf_values
                );
    clear_vector( 
                    input->heads_np * mesh_gp->meshes_np * input->dofs_np,
                    qtf_wl
                );
    clear_vector( 
                    input->heads_np * mesh_gp->meshes_np * input->dofs_np,
                    qtf_bern
                );
    clear_vector( 
                    input->heads_np * mesh_gp->meshes_np * input->dofs_np,
                    qtf_acc
                );
    clear_vector( 
                    input->heads_np * mesh_gp->meshes_np * input->dofs_np,
                    qtf_mom
                );

    // Calculate second order force due to relative wave
    // elevation at the WL
    cuscomplex wl_x( 0.0, 0.0 );
    for ( int i=0; i<input->heads_np; i++ )
    {
        for ( int j=0; j<mesh_gp->meshes_np; j++ )
        {
            idx0 = i * ( input->dofs_np * input->bodies_np ) + j * input->dofs_np;
            for ( int k=pot_gp->field_points_cnp[j]/ngpf_1d; k<pot_gp->field_points_cnp[j+1]/ngpf_1d; k++ )
            {
                panel_k             = mesh_gp->panels_wl[k];
                int_mod             = cuscomplex( 0.0, 0.0 );

                for ( int gpi=0; gpi<input->gauss_order; gpi++ )
                {
                    idx1                = i *  pot_gp->field_points_np + k * ngp + gpi;
                    int_mod             += mdrift_rel_we_i[idx1] * std::conj( mdrift_rel_we_j[idx1] ) * gp.weights[gpi] * panel_k->len_wl / 2.0;
                    wl_x                -= 0.25 * input->grav_acc * input->water_density * int_mod * panel_k->normal_vec[0];

                }


                for ( int r=0; r<input->dofs_np; r++ )
                {
                    qtf_values[idx0+r] += - 0.25 * input->grav_acc * input->water_density * int_mod * panel_k->normal_vec[r];
                    qtf_wl[idx0+r]     += - 0.25 * input->grav_acc * input->water_density * int_mod * panel_k->normal_vec[r];
                }

            }
        }
    }

    // Calculate second order force due to the bernouilly contribution
    for ( int i=0; i<input->heads_np; i++ )
    {
        for ( int j=0; j<mesh_gp->meshes_np; j++ )
        {
            idx0 = i * ( input->dofs_np * input->bodies_np ) + j * input->dofs_np;
            for ( int k= vel_gp->field_points_cnp[j]/ngpf; k<vel_gp->field_points_cnp[j+1]/ngpf; k++ )
            {
                int_mod             = cuscomplex( 0.0, 0.0 );
                panel_k             = mesh_gp->panels[k];

                for ( int gpi=0; gpi<ngp; gpi++ )
                {
                    for ( int gpj=0; gpj<ngp; gpj++ )
                    {
                        idx1        = i * vel_gp->field_points_np + k * pow2s( ngp ) + gpi * ngp + gpj;
                        int_mod     += (
                                            vel_x_i[idx1] * std::conj( vel_x_j[idx1] )
                                            +
                                            vel_y_i[idx1] * std::conj( vel_y_j[idx1] )
                                            +
                                            vel_z_i[idx1] * std::conj( vel_z_j[idx1] )
                                        ) * gp.weights[gpi] * gp.weights[gpj] * jacobi_det_2d( 
                                                                                                    panel_k->num_nodes,
                                                                                                    panel_k->xl,
                                                                                                    panel_k->yl,
                                                                                                    gp.roots[gpi],
                                                                                                    gp.roots[gpj]
                                                                                                );
                    }

                }
                
                for ( int r=0; r<input->dofs_np; r++ )
                {
                    qtf_values[idx0+r] += 0.25 * input->water_density * int_mod * panel_k->normal_vec[r];
                    qtf_bern[idx0+r]   += 0.25 * input->water_density * int_mod * panel_k->normal_vec[r];
                }
            }
        }
    }

    // Calculate second order force due to acceleration term
    cusfloat    cog_to_fp[3];   clear_vector( 3, cog_to_fp );
    cuscomplex  cog_to_fp_c[3]; clear_vector( 3, cog_to_fp );
    cuscomplex  point_disp[3];  clear_vector( 3, point_disp );
    cuscomplex  rao_trans[3];   clear_vector( 3, rao_trans );
    cuscomplex  rao_rot[3];     clear_vector( 3, rao_rot );
    cuscomplex  vel_x_acc;
    cuscomplex  vel_y_acc;
    cuscomplex  vel_z_acc;

    for ( int i=0; i<input->heads_np; i++ )
    {
        for ( int j=0; j<mesh_gp->meshes_np; j++ )
        {
            // Get index to locate RAO data
            idx0 = i * ( input->dofs_np * input->bodies_np ) + j * input->dofs_np;

            // Get RAO values for the current body
            for ( int r=0; r<3; r++ )
            {
                rao_trans[r]    = raos_i[idx0+r] * input->wave_amplitude;
                rao_rot[r]      = raos_i[idx0+3+r] * input->wave_amplitude;
            }

            for ( int k= vel_gp->field_points_cnp[j]/ngpf; k<vel_gp->field_points_cnp[j+1]/ngpf; k++ )
            {
                // Get handle to the current panel for calculations
                panel_k = mesh_gp->panels[k];
                int_mod = cuscomplex( 0.0, 0.0 );

                for ( int gpi=0; gpi<input->gauss_order; gpi++ )
                {
                    for ( int gpj=0; gpj<input->gauss_order; gpj++ )
                    {
                        // Define field points index
                        idx1 = i * vel_gp->field_points_np + k * pow2s( ngp ) + gpi * ngp + gpj;

                        // Define vector from cog to field point
                        sv_sub( 3, &(vel_gp->field_points[3*idx1]), panel_k->body_cog, cog_to_fp );
                        for ( int r=0; r<3; r++ )
                        {
                            cog_to_fp_c[r]  = cuscomplex( cog_to_fp[r], 0.0 );
                        }

                        // Calculate first order displacement of the panel centre
                        clear_vector( 3, point_disp );

                        cross(
                                    rao_rot,
                                    cog_to_fp_c,
                                    point_disp
                            );
                        sv_add(
                                    3,
                                    point_disp,
                                    rao_trans,
                                    point_disp
                                );

                        // Get velocity pressure term
                        vel_x_acc   = input->water_density * cuscomplex( 0.0, -ang_freq_j ) * vel_x_j[idx1];
                        vel_y_acc   = input->water_density * cuscomplex( 0.0, -ang_freq_j ) * vel_y_j[idx1];
                        vel_z_acc   = input->water_density * cuscomplex( 0.0, -ang_freq_j ) * vel_z_j[idx1];

                        // Calculate point displacement
                        
                        int_mod             += 0.25 * (
                                                        point_disp[0] * std::conj( vel_x_acc )
                                                        +
                                                        point_disp[1] * std::conj( vel_y_acc )
                                                        +
                                                        point_disp[2] * std::conj( vel_z_acc )
                                                        +
                                                        std::conj( point_disp[0] ) * vel_x_acc
                                                        +
                                                        std::conj( point_disp[1] ) * vel_y_acc
                                                        +
                                                        std::conj( point_disp[2] ) * vel_z_acc
                                                    ) * gp.weights[gpi] * gp.weights[gpj] * jacobi_det_2d( 
                                                                                                            panel_k->num_nodes,
                                                                                                            panel_k->xl,
                                                                                                            panel_k->yl,
                                                                                                            gp.roots[gpi],
                                                                                                            gp.roots[gpj]
                                                                                                        );
                    }
                }

                cuscomplex int_mod_2        = 0.25 * (
                                                        point_disp[0] * std::conj( vel_x_acc )
                                                        +
                                                        point_disp[1] * std::conj( vel_y_acc )
                                                        +
                                                        point_disp[2] * std::conj( vel_z_acc )
                                                        +
                                                        std::conj( point_disp[0] ) * vel_x_acc
                                                        +
                                                        std::conj( point_disp[1] ) * vel_y_acc
                                                        +
                                                        std::conj( point_disp[2] ) * vel_z_acc
                                                    ) * panel_k->area;

                for ( int r=0; r<input->dofs_np; r++ )
                {
                    qtf_values[idx0+r] += int_mod * panel_k->normal_vec[r];
                    qtf_acc[idx0+r]    += int_mod * panel_k->normal_vec[r];
                }
            }
        }
    }

    // Calculate second order force due to momentum
    cusfloat    ang_2                           = pow2s( ang_freq_j );
    cuscomplex  conj_vec[3];                    clear_vector( 3, conj_vec );
    cuscomplex  hydro_force[input->dofs_np];    clear_vector( input->dofs_np, hydro_force );
    cuscomplex  mom_i[3];                       clear_vector( 3, mom_i );

    for ( int i=0; i<input->heads_np; i++ )
    {
        for ( int j=0; j<mesh_gp->meshes_np; j++ )
        {
            // Define chunk index
            idx0 = i * ( input->dofs_np * input->bodies_np ) + j * input->dofs_np;

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
            cuscomplex      scale_f( 0.25, 0.0 );

            clear_vector(   3,                  mom_i                                                       );
            conj_vector(    3,                  &(raos_i[idx0+3]),      conj_vec                            );
            cross(          conj_vec,           hydro_force,            mom_i                               );
            svs_mult(       3,                  mom_i,                  scale_f,    mom_i                   );
            sv_add(         3,                  &(qtf_values[idx0]),    mom_i,      &(qtf_values[idx0])     );
            sv_add(         3,                  &(qtf_mom[idx0]),       mom_i,      &(qtf_mom[idx0])        );

            clear_vector(   3,                  mom_i                                                       );
            conj_vector(    3,                  hydro_force,            conj_vec                            );
            cross(          &(raos_i[idx0+3]),  conj_vec,               mom_i                               );
            svs_mult(       3,                  mom_i,                  scale_f,    mom_i                   );
            sv_add(         3,                  &(qtf_values[idx0]),    mom_i,      &(qtf_values[idx0])     );
            sv_add(         3,                  &(qtf_mom[idx0]),       mom_i,      &(qtf_mom[idx0])        );

            // Add moment due to rotational force
            clear_vector(   3,                  mom_i                                                       );
            conj_vector(    3,                  &(raos_i[idx0+3]),      conj_vec                            );
            cross(          conj_vec,           &(hydro_force[3]),      mom_i                               );
            svs_mult(       3,                  mom_i,                  scale_f,    mom_i                   );
            sv_add(         3,                  &(qtf_values[idx0+3]),  mom_i,      &(qtf_values[idx0+3])   );
            sv_add(         3,                  &(qtf_mom[idx0+3]),     mom_i,      &(qtf_mom[idx0+3])      );

            clear_vector(   3,                  mom_i                                                       );
            conj_vector(    3,                  &(hydro_force[3]),      conj_vec                            );
            cross(          &(raos_i[idx0+3]),  conj_vec,               mom_i                               );
            svs_mult(       3,                  mom_i,                  scale_f,    mom_i                   );
            sv_add(         3,                  &(qtf_values[idx0+3]),  mom_i,      &(qtf_values[idx0+3])   );
            sv_add(         3,                  &(qtf_mom[idx0+3]),     mom_i,      &(qtf_mom[idx0+3])      );
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

    // if ( std::abs( wij ) > 1e-12 )
    // {
    //     for ( int i=0; i<input->heads_np; i++ )
    //     {
    //         for ( int j=0; j<mesh_gp->meshes_np; j++ )
    //         {
    //              // Define chunk index
    //             idx0 = i * ( input->dofs_np * input->bodies_np ) + j * input->dofs_np;

    //             // Calculate ith, jth wave numbers
    //             ki = w2k( ang_freq_i, input->water_depth, input->grav_acc );
    //             kj = w2k( ang_freq_j, input->water_depth, input->grav_acc );

    //             // Calculate wave numbers difference and the equivalent angular frequency
    //             dk  = ki - kj;
    //             dkwij = k2w( dk, input->water_depth, input->grav_acc );

    //             // Calculate equivalent Froude-Krylov force for frequency difference
    //             calculate_froude_krylov(
    //                                         input,
    //                                         mpi_config,
    //                                         mesh_gp,
    //                                         dkwij,
    //                                         froude_krylov
    //                                     );
                
    //             // Calculate correction coefficients
    //             cij = 2 * ki * kj * wij * ( 1 + std::tanh( ki * h ) * std::tanh( kj * h ) ) / ang_freq_i / ang_freq_j;
    //             bij = (
    //                         pow2s( ki ) / ang_freq_i / pow2s( std::cosh( ki * h ) )
    //                         -
    //                         pow2s( kj ) / ang_freq_j / pow2s( std::cosh( kj * h ) )
    //                     );
    //             aij = (
    //                         ( bij + cij )
    //                         /
    //                         ( pow2s( wij ) - dk * g * std::tanh( dk * h ) )
    //                     ) * 0.5 * g;
    //             fij = ai * aj * aij * wij / g;

    //             // Calcuate second order force contribution using Pinkster approximation
    //             for ( int k=0; k<input->dofs_np; k++ )
    //             {
    //                 qtf_values[idx0+k] = fij * froude_krylov[k];
    //             }
    //         }
    //     }
    // }
}