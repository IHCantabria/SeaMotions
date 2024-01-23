
// Include local modules
#include "qtf_indirect_method.hpp"

#include "../../containers/simulation_data.hpp"
#include "../../containers/matlin_group.hpp"
#include "froude_krylov.hpp"
#include "../../green/kochin.hpp"
#include "../../math/integration.hpp"
#include "../../math/math_tools.hpp"
#include "../../math/special_math.hpp"
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
                                            cuscomplex*     fluid_body_pot_raddif_i,
                                            cuscomplex*     fluid_body_pot_raddif_j,
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
                                            cuscomplex*     fluid_wl_pot_raddif_i,
                                            cuscomplex*     fluid_wl_pot_raddif_j,
                                            cuscomplex*     fluid_wl_vel_x_total_i,
                                            cuscomplex*     fluid_wl_vel_y_total_i,
                                            cuscomplex*     fluid_wl_vel_z_total_i,
                                            cuscomplex*     fluid_wl_vel_x_total_j,
                                            cuscomplex*     fluid_wl_vel_y_total_j,
                                            cuscomplex*     fluid_wl_vel_z_total_j,
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

    /************************************************************************/
    /********************* Calculate first body term ************************/
    /************************************************************************/
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
                                    psi_ds = fluid_body_pot_raddif_i[idx2] - fluid_body_pot_raddif_j[idx2];
                                }
                                else
                                {
                                    psi_ds = fluid_body_pot_raddif_i[idx2] + fluid_body_pot_raddif_j[idx2];
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

    /************************************************************************/
    /******************** Calculate second body term ************************/
    /************************************************************************/
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

    /************************************************************************/
    /******************* Calculate thrid body term WL ***********************/
    /************************************************************************/
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

                for ( int k=wl_gp->field_points_cnp[j]/ngpf_1d; k<wl_gp->field_points_cnp[j+1]/ngpf_1d; k++ )
                {
                    panel_k = mesh_gp->panels_wl[k];
                    int_mod = cuscomplex( 0.0, 0.0 );

                    // Translate normal vector to complex form
                    for ( int r=0; r<3; r++ )
                    {
                        normal_vec_c[r] = cuscomplex( panel_k->normal_vec[r], 0.0 );
                    }

                    for ( int gpi=0; gpi<input->gauss_order; gpi++ )
                    {
                        // Get indexes for heading dependent data
                        idx1_b = k * ngp + gpi;
                        idx1_i = ih1 *  wl_gp->field_points_np + idx1_b;
                        idx1_j = ih2 *  wl_gp->field_points_np + idx1_b;

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
                                                        &(wl_gp->field_points[3*idx1_b]),
                                                        panel_k->body_cog,
                                                        field_point_pos_i
                                                    );

                        calculate_field_point_rot(
                                                        rao_trans_j,
                                                        rao_rot_j,
                                                        &(wl_gp->field_points[3*idx1_b]),
                                                        panel_k->body_cog,
                                                        field_point_pos_j
                                                    );

                        // Check for conjugated quantities
                        if ( qtf_type == 0 )
                        {
                            conj_vector( 3, field_point_pos_j, field_point_pos_j );
                            conj_vector( 3, fluid_body_vel_total_j, fluid_body_vel_total_j );
                        }

                        // Calculate cross product in between field point displacement and the 
                        // total velocity vector
                        cross( field_point_pos_i, fluid_body_vel_total_j, ca_i );
                        cross( field_point_pos_j, fluid_body_vel_total_i, ca_j );

                        // Loop over 6 DOFs in order to calculate the contribution for each of them
                        for ( int r=0; r<input->dofs_np; r++ )
                        {
                            // Get index for radiation values
                            idx2   = r * wl_gp->field_points_nb + idx1_b;

                            // Get radiation potential
                            if ( qtf_type == 0 )
                            {
                                psi_ds = fluid_wl_pot_raddif_i[idx2] - fluid_wl_pot_raddif_j[idx2];
                            }
                            else
                            {
                                psi_ds = fluid_wl_pot_raddif_i[idx2] + fluid_wl_pot_raddif_j[idx2];
                            }

                            // Get integrand value
                            sv_add(  3, ca_i, ca_j, cb_i );
                            svs_mult( 3, cb_i, psi_ds, cb_i );
                            val_mod = sv_dot( 3, cb_i, normal_vec_c );

                            // Get gauss integral chunk value
                            int_mod = val_mod * gp.weights[gpi] * panel_k->len_wl / 2.0;

                            // Add to QTF body force object and scale the gauss integral chunk value
                            qtf_body_force[idx0+r] += cuscomplex( 0.0, - w_ds * rho_w / 2.0 ) * int_mod;
                        }

                    }

                }
            }
        }
    }

    // Delete local heap memory
    delete wdso;
}


void    calculate_qtf_indirect_fs_near_term(
                                                Input*          input,
                                                MeshGroup*      mesh_gp,
                                                int             freq_pos_i,
                                                int             freq_pos_j,
                                                int             qtf_type,
                                                SimulationData* sim_data,
                                                MLGCmpx*        body_gp,
                                                cuscomplex*     qtf_fs_force
                                            )
{
    // Get required fields
    cusfloat    ang_freq_i          = input->angfreqs[freq_pos_i];
    cusfloat    ang_freq_j          = input->angfreqs[freq_pos_j];

    cuscomplex* pot_fk_freq_i       = &( sim_data->qtf_fs_pot_fk_freq[freq_pos_i*sim_data->qtf_fs_heads_np]         );
    cuscomplex* pot_fk_freq_j       = &( sim_data->qtf_fs_pot_fk_freq[freq_pos_j*sim_data->qtf_fs_heads_np]         );

    cuscomplex* pot_raddif_freq_i   = &( sim_data->qtf_fs_pot_raddif_freq[freq_pos_i*sim_data->qtf_fs_raddif_np]    );
    cuscomplex* pot_raddif_freq_j   = &( sim_data->qtf_fs_pot_raddif_freq[freq_pos_j*sim_data->qtf_fs_raddif_np]    );

    cuscomplex* pot_total_freq_i    = &( sim_data->qtf_fs_pot_total_freq[freq_pos_i*sim_data->qtf_fs_heads_np]      );
    cuscomplex* pot_total_freq_j    = &( sim_data->qtf_fs_pot_total_freq[freq_pos_j*sim_data->qtf_fs_heads_np]      );
    
    cuscomplex* vel_x_fk_freq_i     = &( sim_data->qtf_fs_vel_x_fk_freq[freq_pos_i*sim_data->qtf_fs_raddif_np]      );
    cuscomplex* vel_y_fk_freq_i     = &( sim_data->qtf_fs_vel_y_fk_freq[freq_pos_i*sim_data->qtf_fs_raddif_np]      );
    cuscomplex* vel_z_fk_freq_i     = &( sim_data->qtf_fs_vel_z_fk_freq[freq_pos_i*sim_data->qtf_fs_raddif_np]      );
    cuscomplex* vel_x_fk_freq_j     = &( sim_data->qtf_fs_vel_x_fk_freq[freq_pos_j*sim_data->qtf_fs_raddif_np]      );
    cuscomplex* vel_y_fk_freq_j     = &( sim_data->qtf_fs_vel_y_fk_freq[freq_pos_j*sim_data->qtf_fs_raddif_np]      );
    cuscomplex* vel_z_fk_freq_j     = &( sim_data->qtf_fs_vel_z_fk_freq[freq_pos_j*sim_data->qtf_fs_raddif_np]      );

    cuscomplex* vel_x_total_freq_i  = &( sim_data->qtf_fs_vel_x_total_freq[freq_pos_i*sim_data->qtf_fs_heads_np]    );
    cuscomplex* vel_y_total_freq_i  = &( sim_data->qtf_fs_vel_y_total_freq[freq_pos_i*sim_data->qtf_fs_heads_np]    );
    cuscomplex* vel_z_total_freq_i  = &( sim_data->qtf_fs_vel_z_total_freq[freq_pos_i*sim_data->qtf_fs_heads_np]    );
    cuscomplex* vel_x_total_freq_j  = &( sim_data->qtf_fs_vel_x_total_freq[freq_pos_j*sim_data->qtf_fs_heads_np]    );
    cuscomplex* vel_y_total_freq_j  = &( sim_data->qtf_fs_vel_y_total_freq[freq_pos_j*sim_data->qtf_fs_heads_np]    );
    cuscomplex* vel_z_total_freq_j  = &( sim_data->qtf_fs_vel_z_total_freq[freq_pos_j*sim_data->qtf_fs_heads_np]    );

    // Declare local variables
    cuscomplex  f0                  = cuscomplex( 0.0, 0.0 );
    cuscomplex  f1                  = cuscomplex( 0.0, 0.0 );
    cuscomplex  f2                  = cuscomplex( 0.0, 0.0 );
    GaussPoints gp( input->gauss_order );
    cusfloat    grav_acc            = input->grav_acc;
    int         idx0                = 0;
    int         idx1_b              = 0;
    int         idx1_i              = 0;
    int         idx1_j              = 0;
    int         idx2                = 0;
    cuscomplex  int_mod             = cuscomplex( 0.0, 0.0 );
    PanelGeom*  panel_k             = nullptr;
    cuscomplex  pot_fk_i            = cuscomplex( 0.0, 0.0 );
    cuscomplex  pot_fk_j            = cuscomplex( 0.0, 0.0 );
    cuscomplex  pot_pert_i          = cuscomplex( 0.0, 0.0 );
    cuscomplex  pot_pert_j          = cuscomplex( 0.0, 0.0 );
    cuscomplex  pot_total_i         = cuscomplex( 0.0, 0.0 );
    cuscomplex  pot_total_j         = cuscomplex( 0.0, 0.0 );
    cuscomplex  psi_ds              = cuscomplex( 0.0, 0.0 );
    int         ngp                 = input->gauss_order;
    int         ngpf                = input->gauss_np_factor_2d( );
    cusfloat    rho_w               = input->water_density;
    cusfloat    sf                  = 0.0;
    cuscomplex  val_mod_0           = cuscomplex( 0.0, 0.0 );
    cuscomplex  val_mod_1           = cuscomplex( 0.0, 0.0 );
    cuscomplex  vel_fk_i[3];        clear_vector( 3, vel_fk_i );
    cuscomplex  vel_fk_j[3];        clear_vector( 3, vel_fk_j );
    cuscomplex  vel_pert_i[3];      clear_vector( 3, vel_pert_i );
    cuscomplex  vel_pert_j[3];      clear_vector( 3, vel_pert_j );
    cuscomplex  vel_total_i[3];     clear_vector( 3, vel_total_i );
    cuscomplex  vel_total_j[3];     clear_vector( 3, vel_total_j );
    cusfloat    w2_i                = pow2s( ang_freq_i );
    cusfloat    w2_j                = pow2s( ang_freq_j );
    cusfloat    w_ds                = 0.0;

    // Calculate constants depending if there is QTF diff or QTF sum
    // requested
    if ( qtf_type == 0 )
    {
        sf      = 1.0;
        w_ds    = ang_freq_i - ang_freq_j;
    }
    else
    {
        sf      = -1.0;
        w_ds    = ang_freq_i + ang_freq_j;
    }

    // Calculate scaling factors
    f0  = cuscomplex( 0.0, w_ds );
    f1  = - cuscomplex( 0.0, ang_freq_i / 2.0 / grav_acc );
    f2  = sf * cuscomplex( 0.0, ang_freq_j / 2.0 / grav_acc );

    // Clear output results vector
    clear_vector( 
                    pow2s( input->heads_np ) * input->bodies_np * input->dofs_np,
                    qtf_fs_force
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

    cusfloat    k2_i    = pow2s( wdso->k0 );
    cusfloat    k2_j    = pow2s( wdso->k1 );

    // Calculate second order forces for the different headings combinations
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

                    // Loop over gauss points to perform the integration
                    for ( int gpi=0; gpi<ngp; gpi++ )
                    {
                        for ( int gpj=0; gpj<ngp; gpj++ )
                        {
                            // Get panel indexes
                            idx1_b  = k * pow2s( ngp ) + gpi * ngp + gpj;
                            idx1_i  = ih1 * body_gp->field_points_np + idx1_b;
                            idx1_j  = ih2 * body_gp->field_points_np + idx1_b;

                            // Get total fluid velocity components
                            pot_fk_i        = pot_fk_freq_i[idx1_i];
                            pot_fk_j        = pot_fk_freq_j[idx1_j];

                            pot_total_i     = pot_total_freq_i[idx1_i];
                            pot_total_j     = pot_total_freq_j[idx1_j];

                            pot_pert_i      = pot_total_i - pot_fk_i;
                            pot_pert_j      = pot_total_j - pot_fk_j;

                            vel_fk_i[0]     = vel_x_fk_freq_i[idx1_i];
                            vel_fk_i[1]     = vel_y_fk_freq_i[idx1_i];
                            vel_fk_i[2]     = vel_z_fk_freq_i[idx1_i];

                            vel_fk_j[0]     = vel_x_fk_freq_j[idx1_j];
                            vel_fk_j[1]     = vel_y_fk_freq_j[idx1_j];
                            vel_fk_j[2]     = vel_z_fk_freq_j[idx1_j];

                            vel_total_i[0]  = vel_x_total_freq_i[idx1_i];
                            vel_total_i[1]  = vel_y_total_freq_i[idx1_i];
                            vel_total_i[2]  = vel_z_total_freq_i[idx1_i];

                            vel_total_j[0]  = vel_x_total_freq_j[idx1_j];
                            vel_total_j[1]  = vel_y_total_freq_j[idx1_j];
                            vel_total_j[2]  = vel_z_total_freq_j[idx1_j];

                            vel_pert_i[0]   = vel_total_i[0] - vel_fk_i[0];
                            vel_pert_i[1]   = vel_total_i[1] - vel_fk_i[1];
                            vel_pert_i[2]   = vel_total_i[2] - vel_fk_i[2];

                            vel_pert_j[0]   = vel_total_j[0] - vel_fk_j[0];
                            vel_pert_j[1]   = vel_total_j[1] - vel_fk_j[1];
                            vel_pert_j[2]   = vel_total_j[2] - vel_fk_j[2];

                            // Check if the terms needs to be conjugated
                            if ( qtf_type == 0 )
                            {
                                pot_fk_j    = std::conj( pot_fk_j );
                                pot_pert_j  = std::conj( pot_pert_j );
                                pot_total_j = std::conj( pot_total_j );

                                conj_vector( 3, vel_fk_j, vel_fk_j );
                                conj_vector( 3, vel_pert_j, vel_pert_j );
                                conj_vector( 3, vel_total_j, vel_total_j );
                            }

                            // Calculate field points position
                            val_mod_0   =       (
                                                    f0 * ( sv_dot( 3, vel_total_i, vel_pert_j ) + sv_dot( 3, vel_pert_i, vel_fk_j ) )
                                                    +
                                                    f1 * ( pot_total_i * ( - w2_j * vel_pert_j[2] ) + pot_pert_i * ( - w2_j * vel_fk_j[2] ) )
                                                    +
                                                    f2 * ( pot_total_j * ( - w2_i * vel_pert_i[2] ) + pot_pert_j * ( - w2_i * vel_fk_i[2] ) )
                                                );

                            val_mod_1   =  0.5 * (
                                                    -
                                                    cuscomplex( 0.0, ang_freq_i ) * pot_total_i * k2_j * pot_pert_j
                                                    -
                                                    cuscomplex( 0.0, ang_freq_i ) * pot_pert_i * k2_j * pot_fk_j
                                                    +
                                                    sf * cuscomplex( 0.0, ang_freq_j ) * pot_total_j * k2_i * pot_pert_i
                                                    +
                                                    sf * cuscomplex( 0.0, ang_freq_j ) * pot_pert_j * k2_i * pot_fk_i
                                                );

                            // Loop over dofs to get the force for each 6 DOFs component
                            for ( int r=0; r<input->dofs_np; r++ )
                            {
                                // Get index to access the radiation potential data
                                idx2 =  r * body_gp->field_points_np + idx1_b;

                                // Get raddiation fluid velocity components
                                if ( qtf_type == 0 )
                                {
                                    psi_ds = pot_raddif_freq_i[idx2] - pot_raddif_freq_j[idx2];
                                }
                                else
                                {
                                    psi_ds = pot_raddif_freq_i[idx2] + pot_raddif_freq_j[idx2];
                                }

                                // Apply integration weights and multiply by the jacobian determinant
                                int_mod     =  ( val_mod_0 + val_mod_1 ) * gp.weights[gpi] * gp.weights[gpj] * jacobi_det_2d( 
                                                                                                                panel_k->num_nodes,
                                                                                                                panel_k->xl,
                                                                                                                panel_k->yl,
                                                                                                                gp.roots[gpi],
                                                                                                                gp.roots[gpj]
                                                                                                            );

                                // Get total force
                                qtf_fs_force[idx0+r]    += cuscomplex( 0.0, w_ds * rho_w / grav_acc ) * psi_ds * int_mod;
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


void    calculate_qtf_indirect_fs_far_term(
                                                Input*          input,
                                                MeshGroup*      mesh_gp,
                                                int             freq_pos_i,
                                                int             freq_pos_j,
                                                int             qtf_type,
                                                SimulationData* sim_data,
                                                MLGCmpx*        body_gp,
                                                cuscomplex*     qtf_fs_force
                                            )
{
    // Define local variables
    int         idx0                = 0;

    // Get required fields
    cusfloat    ang_freq_i          = input->angfreqs[freq_pos_i];
    cusfloat    ang_freq_j          = input->angfreqs[freq_pos_j];

    cusfloat    grav_acc            = input->grav_acc;
    cusfloat    rho_w               = input->water_density;

    // Define second order wave dispersion object
    WaveDispersionSO*   wdso        = new WaveDispersionSO( 
                                                            input->wave_amplitude,
                                                            input->wave_amplitude,
                                                            ang_freq_i,
                                                            ang_freq_j,
                                                            input->heads[0],
                                                            input->heads[0],
                                                            input->water_depth,
                                                            input->grav_acc
                                                        );

    // Define common scaling factors
    cusfloat    k_ds    = 0.0;
    cusfloat    sf0     = 0.0;
    cusfloat    sf1     = 0.0;
    cusfloat    w_ds    = 0.0;

    if ( qtf_type == 0 )
    {
        k_ds    = wdso->k_diff_mod;
        sf0     = 1.0;
        sf1     = - 1.0;
        w_ds    = wdso->w_diff;
    }
    else
    {
        k_ds    = wdso->k_sum_mod;
        sf0     = - 1.0;
        sf1     = 1.0;
        w_ds    = wdso->w_sum;
    }

    cuscomplex  sf2     = cuscomplex( 0.0, 1.0 ) * w_ds * rho_w / grav_acc;
    cuscomplex  sf3     = (
                                cuscomplex( 0.0, 1.0 )
                                *
                                input->wave_amplitude
                                *
                                grav_acc
                                *
                                8.0 
                                * 
                                PI
                                *
                                std::sqrt( k_ds )
                                *
                                wave_vertical_profile_mod_fo( 
                                                                k_ds,
                                                                input->water_depth,
                                                                0.0
                                                            )
                            );

    cuscomplex  kappa_1 = sf0 * cuscomplex( 0.0, -1.0 ) * w_ds * wdso->k0 * wdso->k1;
    cuscomplex  kappa_2 =  (
                                (
                                    cuscomplex( 0.0, 1.0 )
                                    *
                                    w_ds
                                    *
                                    pow2s( ang_freq_i ) / grav_acc
                                    *
                                    pow2s( ang_freq_j ) / grav_acc
                                )
                                +
                                sf0 * (
                                            cuscomplex( 0.0, 1.0 ) * ang_freq_i * ang_freq_j / 2.0
                                            *
                                            (
                                                pow2s( wdso->k0 ) / ( ang_freq_i * pow2s( std::cosh( wdso->k0 * input->water_depth ) ) )
                                                +
                                                sf1 * pow2s( wdso->k1 ) / ( ang_freq_j * pow2s( std::cosh( wdso->k1 * input->water_depth ) ) )
                                            )
                                        )
                            );

    // Calculate second order forces for the different headings combinations
    for ( int ih1=0; ih1<input->heads_np; ih1++ )
    {
        for ( int ih2=0; ih2<input->heads_np; ih2++ )
        {
            
            for ( int j=0; j<mesh_gp->meshes_np; j++ )
            {
                // Get index for the storage of the body force
                idx0 = (
                            ih1 * ( input->dofs_np * input->bodies_np * input->heads_np )
                            +
                            ih2 * ( input->dofs_np * input->bodies_np )
                            + 
                            j * input->dofs_np
                        );

                //
            }
        }
    }

    // Delete local heap memory
    delete wdso; 
}


cuscomplex  calculate_r1_integral(
                                    cusfloat            R,
                                    WaveDispersionSO*   wdso,
                                    int                 l_order,
                                    int                 qtf_type
                                )
{
    // Calculate scaling factors
    cuscomplex sf0( 1.0, 0.0 );
    if ( qtf_type == 1 )
    {
        sf0 = std::exp( cuscomplex( 0.0, PI / 2.0 ) );
    }

    cusfloat    ep_l    = ep_n( l_order );
    cuscomplex  sfi     = std::pow( cuscomplex( 0.0, 1.0 ), l_order );

    // Calculate alpha and beta parameters
    cusfloat    alpha   = 0.0;
    cusfloat    beta    = wdso->k0;

    if ( qtf_type == 0 )
    {
        alpha   = wdso->k_diff_mod - wdso->k1;
    }
    else
    {
        alpha   = wdso->k_sum_mod + wdso->k1;
    }

    // Calculate 0 to inf analytical integral
    cuscomplex  int_value_ana   = besseljn_expi_int( 
                                                        alpha,
                                                        beta,
                                                        static_cast<cusfloat>( l_order )
                                                    );

    // Calculate 0 to R numerical integral
    cusfloat    cos_value_num   = romberg_quadrature(
                                                        besseljn_cos_kernel,
                                                        0.0,
                                                        R,
                                                        1e-6
                                                    );
    
    cusfloat    sin_value_num   = romberg_quadrature(
                                                        besseljn_sin_kernel,
                                                        0.0,
                                                        R,
                                                        1e-6
                                                    );

    cuscomplex  int_value_num   = cuscomplex( cos_value_num, sin_value_num );

    // Calculate R to infinite integral
    cuscomplex  int_value       = sf0 * ep_l * sfi * ( int_value_ana - int_value_num );

    return int_value;
}


cuscomplex  calculate_r2_integral(
                                    cusfloat            R,
                                    WaveDispersionSO*   wdso,
                                    int                 l_order,
                                    int                 qtf_type
                                )
{
    // Define scaling factors
    cuscomplex  sf0     = std::exp( cuscomplex( 0.0, PI / 2.0 ) );
    cusfloat    ep_l    = ep_n( l_order );
    cuscomplex  sfi     = std::pow( cuscomplex( 0.0, 1.0 ), l_order );

    if ( qtf_type == 0 )
    {
        sfi = std::conj( sfi );
    }

    // Calculate alpha and beta
    cusfloat    alpha   = 0.0;
    cusfloat    beta    = wdso->k1;
    if ( qtf_type == 0 )
    {
        alpha = wdso->k_diff_mod + wdso->k0;
    }
    else
    {
        alpha = wdso->k_sum_mod + wdso->k0;
    }

    // Calculate 0 to inf analytical integral
    cuscomplex  int_value_ana   = besseljn_expi_int( 
                                                        alpha,
                                                        beta,
                                                        static_cast<cusfloat>( l_order )
                                                    );

    // Calculate 0 to R numerical integral
    cusfloat    cos_value_num   = romberg_quadrature(
                                                        besseljn_cos_kernel,
                                                        0.0,
                                                        R,
                                                        1e-6
                                                    );
    
    cusfloat    sin_value_num   = romberg_quadrature(
                                                        besseljn_sin_kernel,
                                                        0.0,
                                                        R,
                                                        1e-6
                                                    );

    cuscomplex  int_value_num   = cuscomplex( cos_value_num, sin_value_num );

    // Calculate R to infinite integral
    cuscomplex  int_value       = sf0 * ep_l * sfi * ( int_value_ana - int_value_num );

    return int_value;
}


cuscomplex  calculate_theta_integral(
                                        Input*      input,
                                        cusfloat    beta,
                                        int         l_order,
                                        int         qtf_type,
                                        cuscomplex* kochin_cos_pert_j,
                                        cuscomplex* kochin_sin_pert_j,
                                        cuscomplex* kochin_cos_rad_i,
                                        cuscomplex* kochin_cos_rad_j,
                                        cuscomplex* kochin_sin_rad_i,
                                        cuscomplex* kochin_sin_rad_j,
                                        cuscomplex* body_force
                                    )
{
    // Clear incoming vector in order to avoid taking into account spurious data
    clear_vector( 
                    input->heads_np * input->bodies_np * input->dofs_np,
                    body_force
                );

    // Define local variables
    cuscomplex  cpm     = cuscomplex( 0.0, 0.0 );
    cuscomplex  spm     = cuscomplex( 0.0, 0.0 );
    cuscomplex  crn     = cuscomplex( 0.0, 0.0 );
    cuscomplex  srn     = cuscomplex( 0.0, 0.0 );
    int         idx0    = 0;
    int         idx1    = 0;
    cusfloat    t0_int  = 0.0;
    cusfloat    t1_int  = 0.0;
    cusfloat    t2_int  = 0.0;
    cusfloat    t3_int  = 0.0;

    // Loop over values to calculate the resulting forces
    for ( int ih=0; ih<input->heads_np; ih++ )
    {
        for ( int ib=0; ib<input->bodies_np; ib++ )
        {
            for ( int id=0; id<input->dofs_np; id++ )
            {
                // Define index to stoage the data
                idx0    = ( 
                                ih * ( input->bodies_np * input->dofs_np )
                                +
                                ib * input->dofs_np
                                +
                                id
                            );

                // Loop over kochin coefficients to calculate the force
                for ( int m=0; m<input->kochin_np; m++ )
                {
                    // Get current index for the mth perturbation coefficient
                    idx1 = (
                                ih * ( input->bodies_np * input->kochin_np )
                                +
                                ib * input->kochin_np
                                +
                                m
                            );
                    
                    // Get perturbation series coefficients
                    if ( qtf_type == 0 )
                    {
                        cpm = std::conj( kochin_cos_pert_j[idx1] );
                        spm = std::conj( kochin_sin_pert_j[idx1] );
                    }
                    else
                    {
                        cpm = kochin_cos_pert_j[idx1];
                        spm = kochin_sin_pert_j[idx1];
                    }

                    
                    for ( int n=0; n<input->kochin_np; n++ )
                    {
                        // Get current index for the nth perturbation coefficient
                        idx1 = (
                                    id * ( input->bodies_np * input->kochin_np )
                                    +
                                    ib * input->kochin_np
                                    +
                                    m
                                );

                        // Get radiation coefficient
                        if ( qtf_type == 0 )
                        {
                            crn = ( kochin_cos_rad_i[idx1] - kochin_cos_rad_j[idx1] );
                            srn = ( kochin_sin_rad_i[idx1] - kochin_sin_rad_j[idx1] );
                        }
                        else
                        {
                            crn = ( kochin_cos_rad_i[idx1] + kochin_cos_rad_j[idx1] );
                            srn = ( kochin_sin_rad_i[idx1] + kochin_sin_rad_j[idx1] );
                        }

                        // Calculate angular integrals
                        t0_int  = (
                                        calculate_kochin_cosexp_t0( 2*PI, beta, l_order, m, n )
                                        -
                                        calculate_kochin_cosexp_t0( 0.0, beta, l_order, m, n )
                                    );

                        t1_int  = (
                                        calculate_kochin_cosexp_t1( 2*PI, beta, l_order, m, n )
                                        -
                                        calculate_kochin_cosexp_t1( 0.0, beta, l_order, m, n )
                                    );

                        t2_int  = (
                                        calculate_kochin_cosexp_t2( 2*PI, beta, l_order, m, n )
                                        -
                                        calculate_kochin_cosexp_t2( 0.0, beta, l_order, m, n )
                                    );

                        t3_int  = (
                                        calculate_kochin_cosexp_t3( 2*PI, beta, l_order, m, n )
                                        -
                                        calculate_kochin_cosexp_t3( 0.0, beta, l_order, m, n )
                                    );

                        // Calculate body force
                        body_force[idx0] += (
                                                cpm * crn * t0_int
                                                +
                                                cpm * srn * t1_int
                                                +
                                                spm * crn * t2_int
                                                +
                                                spm * srn * t3_int
                                            );
                    }
                }
            }
        }
    }
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