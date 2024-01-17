
// Include local modules
#include "qtf_indirect_method.hpp"

#include "../../containers/simulation_data.hpp"
#include "../../containers/matlin_group.hpp"
#include "froude_krylov.hpp"
#include "../../math/integration.hpp"
#include "../../math/math_tools.hpp"
#include "../../waves/waves_common.hpp"
#include "../../waves/wave_dispersion_so.hpp"
#include "tools.hpp"


void    calculate_qtf_indirect_body_term(
                                            Input*          input,
                                            MeshGroup*      mesh_gp,
                                            cusfloat        ang_freq_i,
                                            cusfloat        ang_freq_j,
                                            int             qtf_type,
                                            cuscomplex*     raos_i,
                                            cuscomplex*     raos_j,
                                            cuscomplex*     fluid_body_pot_i,
                                            cuscomplex*     fluid_body_pot_j,
                                            cuscomplex*     fluid_body_vel_raddif_x_i,
                                            cuscomplex*     fluid_body_vel_raddif_y_i,
                                            cuscomplex*     fluid_body_vel_raddif_z_i,
                                            cuscomplex*     fluid_body_vel_raddif_x_j,
                                            cuscomplex*     fluid_body_vel_raddif_y_j,
                                            cuscomplex*     fluid_body_vel_raddif_z_j,
                                            cuscomplex*     fluid_body_vel_total_x_i,
                                            cuscomplex*     fluid_body_vel_total_y_i,
                                            cuscomplex*     fluid_body_vel_total_z_i,
                                            cuscomplex*     fluid_body_vel_total_x_j,
                                            cuscomplex*     fluid_body_vel_total_y_j,
                                            cuscomplex*     fluid_body_vel_total_z_j,
                                            SimulationData* sim_data,
                                            MLGCmpx*        body_gp,
                                            MLGCmpx*        wl_gp,
                                            cuscomplex*     qtf_body_force
                                        )
{
    // Clear output results vector
    clear_vector( 
                    pow2s( input->heads_np ) * input->bodies_np * input->dofs_np,
                    qtf_body_force
                );

    // Define second order wave dispersion object
    WaveDispersionSO*   wdso    = new WaveDispersionSO( 
                                                        input->wave_amplitude,
                                                        input->wave_amplitude,
                                                        ang_freq_i,
                                                        ang_freq_j,
                                                        input->heads[0],
                                                        input->heads[0],
                                                        input->water_depth,
                                                        input->grav_acc
                                                    );
    
    // Define aux variables to be used along the function
    cuscomplex  ca_i[3];                    clear_vector( 3, ca_i );
    cuscomplex  ca_j[3];                    clear_vector( 3, ca_j );
    cuscomplex  cb_i[3];                    clear_vector( 3, cb_i );
    cuscomplex  cb_j[3];                    clear_vector( 3, cb_j );
    GaussPoints gp( input->gauss_order );
    cuscomplex  field_point_pos_i[3];       clear_vector( 3, field_point_pos_i );
    cuscomplex  field_point_pos_j[3];       clear_vector( 3, field_point_pos_j );
    cuscomplex  field_point_vel_i[3];       clear_vector( 3, field_point_vel_i );
    cuscomplex  field_point_vel_j[3];       clear_vector( 3, field_point_vel_j );
    cuscomplex  fluid_body_vel_raddif[3];   clear_vector( 3, fluid_body_vel_raddif );
    cuscomplex  fluid_body_vel_total_i[3];  clear_vector( 3, fluid_body_vel_total_i );
    cuscomplex  fluid_body_vel_total_j[3];  clear_vector( 3, fluid_body_vel_total_j );
    int         idx0                        = 0;
    int         idx1_b                      = 0;
    int         idx1_i                      = 0;
    int         idx1_j                      = 0;
    int         idx2                        = 0;
    cuscomplex  int_mod                     = 0.0;
    cuscomplex  jac_rot_i[9];               clear_vector( 9, jac_rot_i );
    cuscomplex  jac_rot_j[9];               clear_vector( 9, jac_rot_j );
    PanelGeom*  panel_k                     = nullptr;
    int         ngp                         = input->gauss_order;
    int         ngpf_1d                     = input->gauss_np_factor_1d( );
    int         ngpf                        = input->gauss_np_factor_2d( );
    cuscomplex  normal_rot_i[3];            clear_vector( 3, normal_rot_i );
    cuscomplex  normal_rot_j[3];            clear_vector( 3, normal_rot_j );
    cuscomplex  normal_vec_c[3];            clear_vector( 3, normal_vec_c );
    cuscomplex  psi_ds                      = cuscomplex( 0.0, 0.0 );
    cuscomplex  rao_rot_i[3];               clear_vector( 3, rao_rot_i );
    cuscomplex  rao_rot_j[3];               clear_vector( 3, rao_rot_j );
    cuscomplex  rao_trans_i[3];             clear_vector( 3, rao_trans_i );
    cuscomplex  rao_trans_j[3];             clear_vector( 3, rao_trans_j );
    cusfloat    rho_w                       = input->water_density;
    cuscomplex  val_mod                     = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_inc[3];
    cusfloat    w_ds                        = 0.0;

    // Get difference or summation 
    if ( qtf_type == 0 )
    {
        w_ds = ang_freq_i - ang_freq_j;
    }
    else if ( qtf_type == 1 )
    {
        w_ds = ang_freq_i + ang_freq_j;
    }

    // Calculate first body term
    for ( int ih1=0; ih1<input->heads_np; ih1++ )
    {
        for ( int ih2=0; ih2<input->heads_np; ih2++ )
        {
            
            for ( int j=0; j<mesh_gp->meshes_np; j++ )
            {
                
                idx0 = (
                            ih1 * ( input->dofs_np * input->bodies_np * input->heads_np )
                            +
                            ih2 * ( input->dofs_np * input->bodies_np )
                            + 
                            j * input->dofs_np
                        );
                        
                

                for ( int k=body_gp->field_points_cnp[j]/ngpf; k<body_gp->field_points_cnp[j+1]/ngpf; k++ )
                {
                    int_mod             = cuscomplex( 0.0, 0.0 );
                    panel_k             = mesh_gp->panels[k];

                    // Translate normal vector to complex form
                    for ( int r=0; r<3; r++ )
                    {
                        normal_vec_c[r] = cuscomplex( panel_k->normal_vec[r], 0.0 );
                    }

                    // Loop over gauss points to perform the integration
                    for ( int gpi=0; gpi<ngp; gpi++ )
                    {
                        for ( int gpj=0; gpj<ngp; gpj++ )
                        {
                            // Get panel indexes
                            idx1_b  = k * pow2s( ngp ) + gpi * ngp + gpj;
                            idx1_i  = ih1 * body_gp->field_points_np + idx1_b;
                            idx1_j  = ih2 * body_gp->field_points_np + idx1_b;

                            // Calculate wave incident term
                            wave_inc[0] =   wave_potential_so_space_dx( 
                                                                            body_gp->field_points[3*idx1_b],
                                                                            body_gp->field_points[3*idx1_b+1],
                                                                            body_gp->field_points[3*idx1_b+2],
                                                                            wdso,
                                                                            qtf_type
                                                                        );

                            wave_inc[1] =   wave_potential_so_space_dy( 
                                                                            body_gp->field_points[3*idx1_b],
                                                                            body_gp->field_points[3*idx1_b+1],
                                                                            body_gp->field_points[3*idx1_b+2],
                                                                            wdso,
                                                                            qtf_type
                                                                        );

                            wave_inc[2] =   wave_potential_so_space_dz( 
                                                                            body_gp->field_points[3*idx1_b],
                                                                            body_gp->field_points[3*idx1_b+1],
                                                                            body_gp->field_points[3*idx1_b+2],
                                                                            wdso,
                                                                            qtf_type
                                                                        );

                            val_mod     =   (
                                                wave_inc[0] * mesh_gp->panels[k]->normal_vec[0]
                                                +
                                                wave_inc[1] * mesh_gp->panels[k]->normal_vec[1]
                                                +
                                                wave_inc[2] * mesh_gp->panels[k]->normal_vec[2]
                                            );

                            // Get RAO values for the current body
                            for ( int r=0; r<3; r++ )
                            {
                                rao_rot_i[r]    = raos_i[idx1_i+3+r] * input->wave_amplitude;
                                rao_rot_j[r]    = raos_j[idx1_j+3+r] * input->wave_amplitude;
                                rao_trans_i[r]  = raos_i[idx1_i+r]   * input->wave_amplitude;
                                rao_trans_j[r]  = raos_j[idx1_j+r]   * input->wave_amplitude;
                            }

                            // Calculate panel field points velocity
                            calculate_field_point_vel_rot(
                                                            rao_trans_i,
                                                            rao_rot_i,
                                                            &(body_gp->field_points[3*idx1_b]),
                                                            panel_k->body_cog,
                                                            ang_freq_i,
                                                            field_point_vel_i
                                                        );

                            calculate_field_point_vel_rot(
                                                            rao_trans_j,
                                                            rao_rot_j,
                                                            &(body_gp->field_points[3*idx1_b]),
                                                            panel_k->body_cog,
                                                            ang_freq_j,
                                                            field_point_vel_j
                                                        );

                            // Calculate the rotated normal vector
                            cross(
                                        rao_rot_i,
                                        normal_vec_c,
                                        normal_rot_i
                                    );

                            cross(
                                        rao_rot_j,
                                        normal_vec_c,
                                        normal_rot_j
                                    );

                            // Get total fluid velocity components
                            fluid_body_vel_total_i[0]   = fluid_body_vel_total_x_i[idx1_i];
                            fluid_body_vel_total_i[1]   = fluid_body_vel_total_y_i[idx1_i];
                            fluid_body_vel_total_i[2]   = fluid_body_vel_total_z_i[idx1_i];

                            fluid_body_vel_total_j[0]   = fluid_body_vel_total_x_j[idx1_j];
                            fluid_body_vel_total_j[1]   = fluid_body_vel_total_y_j[idx1_j];
                            fluid_body_vel_total_j[2]   = fluid_body_vel_total_z_j[idx1_j];

                            // Calculate velocity term
                            if ( qtf_type == 0 )
                            {
                                conj_vector( 3, normal_rot_i, normal_rot_i );
                                conj_vector( 3, field_point_vel_j, field_point_vel_j );
                                conj_vector( 3, fluid_body_vel_total_j, fluid_body_vel_total_j );
                            }

                            // Calcualte difference in between total rigid body point 
                            // velocity and field velocity
                            sv_sub( 
                                        3,
                                        field_point_vel_i,
                                        fluid_body_vel_total_i,
                                        ca_i
                                    );

                            sv_sub(
                                        3,
                                        field_point_vel_j,
                                        fluid_body_vel_total_j,
                                        ca_j
                                    );

                            // Calculate field points position
                            val_mod     -=  0.5 * (
                                                sv_dot( 3, ca_i, field_point_vel_i )
                                                +
                                                sv_dot( 3, ca_j, field_point_vel_j )
                                            );

                            // Apply integration weights and multiply by the jacobian determinant
                            int_mod     =  val_mod * gp.weights[gpi] * gp.weights[gpj] * jacobi_det_2d( 
                                                                                                            panel_k->num_nodes,
                                                                                                            panel_k->xl,
                                                                                                            panel_k->yl,
                                                                                                            gp.roots[gpi],
                                                                                                            gp.roots[gpj]
                                                                                                        );

                            // Loop over dofs to get the force for each 6 DOFs component
                            for ( int r=0; r<input->dofs_np; r++ )
                            {
                                // Get index to access the radiation potential data
                                idx2 =  r * body_gp->field_points_np + idx1_b;

                                // Get radiation potential
                                if ( qtf_type == 0 )
                                {
                                    psi_ds = fluid_body_pot_i[idx2] - fluid_body_pot_j[idx2];
                                }
                                else
                                {
                                    psi_ds = fluid_body_pot_i[idx2] + fluid_body_pot_j[idx2];
                                }

                                // Get total force
                                qtf_body_force[idx0+r] += cuscomplex( 0.0, w_ds * rho_w ) * psi_ds * int_mod;
                            }
                        }

                    }
                }
            }
        }
    }

    // Calculate second body term
    for ( int ih1=0; ih1<input->heads_np; ih1++ )
    {
        for ( int ih2=0; ih2<input->heads_np; ih2++ )
        {
            
            for ( int j=0; j<mesh_gp->meshes_np; j++ )
            {
                
                idx0 = (
                            ih1 * ( input->dofs_np * input->bodies_np * input->heads_np )
                            +
                            ih2 * ( input->dofs_np * input->bodies_np )
                            + 
                            j * input->dofs_np
                        );
                        
                

                for ( int k=body_gp->field_points_cnp[j]/ngpf; k<body_gp->field_points_cnp[j+1]/ngpf; k++ )
                {
                    int_mod             = cuscomplex( 0.0, 0.0 );
                    panel_k             = mesh_gp->panels[k];

                    // Translate normal vector to complex form
                    for ( int r=0; r<3; r++ )
                    {
                        normal_vec_c[r] = cuscomplex( panel_k->normal_vec[r], 0.0 );
                    }

                    // Loop over gauss points to perform the integration
                    for ( int gpi=0; gpi<ngp; gpi++ )
                    {
                        for ( int gpj=0; gpj<ngp; gpj++ )
                        {
                            // Get panel indexes
                            idx1_b  = k * pow2s( ngp ) + gpi * ngp + gpj;
                            idx1_i  = ih1 * body_gp->field_points_np + idx1_b;
                            idx1_j  = ih2 * body_gp->field_points_np + idx1_b;

                            // Get RAO values for the current body
                            for ( int r=0; r<3; r++ )
                            {
                                rao_rot_i[r]    = raos_i[idx1_i+3+r] * input->wave_amplitude;
                                rao_rot_j[r]    = raos_j[idx1_j+3+r] * input->wave_amplitude;
                                rao_trans_i[r]  = raos_i[idx1_i+r]   * input->wave_amplitude;
                                rao_trans_j[r]  = raos_j[idx1_j+r]   * input->wave_amplitude;
                            }

                            // Get total fluid velocity components
                            fluid_body_vel_total_i[0]   = fluid_body_vel_total_x_i[idx1_i];
                            fluid_body_vel_total_i[1]   = fluid_body_vel_total_y_i[idx1_i];
                            fluid_body_vel_total_i[2]   = fluid_body_vel_total_z_i[idx1_i];

                            fluid_body_vel_total_j[0]   = fluid_body_vel_total_x_j[idx1_j];
                            fluid_body_vel_total_j[1]   = fluid_body_vel_total_y_j[idx1_j];
                            fluid_body_vel_total_j[2]   = fluid_body_vel_total_z_j[idx1_j];

                            // Calculate panel field points displacement
                            calculate_field_point_rot(
                                                            rao_trans_i,
                                                            rao_rot_i,
                                                            &(body_gp->field_points[3*idx1_b]),
                                                            panel_k->body_cog,
                                                            field_point_pos_i
                                                        );

                            calculate_field_point_rot(
                                                            rao_trans_j,
                                                            rao_rot_j,
                                                            &(body_gp->field_points[3*idx1_b]),
                                                            panel_k->body_cog,
                                                            field_point_pos_j
                                                        );

                            // Loop over dofs to get the force for each 6 DOFs component
                            for ( int r=0; r<input->dofs_np; r++ )
                            {
                                // Get index to access the radiation potential data
                                idx2 =  r * body_gp->field_points_np + idx1_b;

                                // Get raddiation fluid velocity components
                                if ( qtf_type == 0 )
                                {
                                    fluid_body_vel_raddif[0]  = fluid_body_vel_raddif_x_i[idx1_i] - fluid_body_vel_raddif_x_j[idx1_j];
                                    fluid_body_vel_raddif[1]  = fluid_body_vel_raddif_y_i[idx1_i] - fluid_body_vel_raddif_y_j[idx1_j];
                                    fluid_body_vel_raddif[2]  = fluid_body_vel_raddif_z_i[idx1_i] - fluid_body_vel_raddif_z_j[idx1_j];
                                }
                                else
                                {
                                    fluid_body_vel_raddif[0]  = fluid_body_vel_raddif_x_i[idx1_i] + fluid_body_vel_raddif_x_j[idx1_j];
                                    fluid_body_vel_raddif[1]  = fluid_body_vel_raddif_y_i[idx1_i] + fluid_body_vel_raddif_y_j[idx1_j];
                                    fluid_body_vel_raddif[2]  = fluid_body_vel_raddif_z_i[idx1_i] + fluid_body_vel_raddif_z_j[idx1_j];
                                }

                                // Check if the terms needs to be conjugated
                                if ( qtf_type == 0 )
                                {
                                    conj_vector( 3, field_point_pos_j, field_point_pos_j );
                                    conj_vector( 3, fluid_body_vel_total_j, fluid_body_vel_total_j );
                                }

                                // Calculate field points position
                                val_mod     =  (
                                                    sv_dot( 3, fluid_body_vel_total_j, fluid_body_vel_raddif ) * sv_dot( 3, field_point_pos_i, normal_vec_c )
                                                    +
                                                    sv_dot( 3, fluid_body_vel_total_i, fluid_body_vel_raddif ) * sv_dot( 3, field_point_pos_j, normal_vec_c )
                                                    -
                                                    sv_dot( 3, fluid_body_vel_raddif, field_point_pos_i ) * sv_dot( 3, fluid_body_vel_total_j, normal_vec_c )
                                                    -
                                                    sv_dot( 3, fluid_body_vel_raddif, field_point_pos_j ) * sv_dot( 3, fluid_body_vel_total_i, normal_vec_c )
                                                );

                                // Apply integration weights and multiply by the jacobian determinant
                                int_mod     =  val_mod * gp.weights[gpi] * gp.weights[gpj] * jacobi_det_2d( 
                                                                                                                panel_k->num_nodes,
                                                                                                                panel_k->xl,
                                                                                                                panel_k->yl,
                                                                                                                gp.roots[gpi],
                                                                                                                gp.roots[gpj]
                                                                                                            );

                                // Get total force
                                qtf_body_force[idx0+r] += cuscomplex( 0.0, w_ds * rho_w / 2.0 ) * int_mod;
                            }
                        }
                    }
                }
            }
        }
    }

    // Delete local heap memory
    delete wdso;
}


void    calculate_secord_force_indirect(
                                            Input*      input,
                                            MeshGroup*  mesh_gp,
                                            cusfloat    ang_freq_i,
                                            cusfloat    ang_freq_j,
                                            int         qtf_type,
                                            cuscomplex* froude_krylov,
                                            cuscomplex* body_force,
                                            cuscomplex* fs_near_field,
                                            cuscomplex* fs_far_field,
                                            cuscomplex* secord_force_total
                                        )
{
    // Clear input data
    int data_np = pow2s( input->heads_np ) * input->bodies_np * input->dofs_np;

    clear_vector( data_np, froude_krylov        );
    clear_vector( data_np, body_force           );
    clear_vector( data_np, fs_near_field        );
    clear_vector( data_np, fs_far_field         );
    clear_vector( data_np, secord_force_total   );

    // Calculate second order Froude-Krylov force
    calculate_froude_krylov_so(
                                    input,
                                    mesh_gp,
                                    ang_freq_i,
                                    ang_freq_j,
                                    qtf_type,
                                    froude_krylov
                                );

    // Calculate diffraction force due to body term

    // Calculate diffraction force due to the near field free surface term

    // Calculate diffraction force due to the far field free surface term

    // Sum-up contributions
    sv_add( data_np, secord_force_total, froude_krylov, secord_force_total );
    sv_add( data_np, secord_force_total, body_force,    secord_force_total );
    sv_add( data_np, secord_force_total, fs_near_field, secord_force_total );
    sv_add( data_np, secord_force_total, fs_far_field,  secord_force_total );
    
}