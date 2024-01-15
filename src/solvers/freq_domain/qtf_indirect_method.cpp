
// Include local modules
#include "qtf_indirect_method.hpp"

#include "../../containers/simulation_data.hpp"
#include "froude_krylov.hpp"
#include "../../math/integration.hpp"


void    calculate_qtf_indirect_body_term(
                                            Input*          input,
                                            MeshGroup*      mesh_gp,
                                            cusfloat        ang_freq_i,
                                            cusfloat        ang_freq_j,
                                            SimulationData* sim_data,
                                            cuscomplex*     qtf_body_force
                                        )
{
    // Define aux variables to be used along the function
    GaussPoints gp( input->gauss_order );
    int         idx0        = 0;
    int         idx1_i      = 0;
    int         idx1_j      = 0;

    // Calculate first body term
    // for ( int ih1=0; ih1<input->heads_np; ih1++ )
    // {
    //     for ( int ih2=0; ih2<input->heads_np; ih2++ )
    //     {
    //         for ( int j=0; j<mesh_gp->meshes_np; j++ )
    //         {
                
    //             idx0 = (
    //                         ih1 * ( input->dofs_np * input->bodies_np * input->heads_np )
    //                         +
    //                         ih2 * ( input->dofs_np * input->bodies_np )
    //                         + 
    //                         j * input->dofs_np
    //                     );
                        
                

    //             for ( int k=vel_gp->field_points_cnp[j]/ngpf; k<vel_gp->field_points_cnp[j+1]/ngpf; k++ )
    //             {
    //                 int_mod             = cuscomplex( 0.0, 0.0 );
    //                 panel_k             = mesh_gp->panels[k];

    //                 for ( int gpi=0; gpi<ngp; gpi++ )
    //                 {
    //                     for ( int gpj=0; gpj<ngp; gpj++ )
    //                     {
    //                             idx1_i    = ih1 * vel_gp->field_points_np + k * pow2s( ngp ) + gpi * ngp + gpj;
    //                             idx1_j    = ih2 * vel_gp->field_points_np + k * pow2s( ngp ) + gpi * ngp + gpj;
                            

    //                         if ( qtf_type == 0 )
    //                         {
    //                             val_mod     =   (
    //                                                 vel_x_i[idx1_i] * std::conj( vel_x_j[idx1_j] )
    //                                                 +
    //                                                 vel_y_i[idx1_i] * std::conj( vel_y_j[idx1_j] )
    //                                                 +
    //                                                 vel_z_i[idx1_i] * std::conj( vel_z_j[idx1_j] )
    //                                             );
    //                         }
    //                         else
    //                         {
    //                             val_mod     =   (
    //                                                 vel_x_i[idx1_i] * vel_x_j[idx1_j]
    //                                                 +
    //                                                 vel_y_i[idx1_i] * vel_y_j[idx1_j]
    //                                                 +
    //                                                 vel_z_i[idx1_i] * vel_z_j[idx1_j]
    //                                             );
    //                         }
    //                         int_mod     +=  val_mod * gp.weights[gpi] * gp.weights[gpj] * jacobi_det_2d( 
    //                                                                                                         panel_k->num_nodes,
    //                                                                                                         panel_k->xl,
    //                                                                                                         panel_k->yl,
    //                                                                                                         gp.roots[gpi],
    //                                                                                                         gp.roots[gpj]
    //                                                                                                     );
    //                     }

    //                 }
                    
    //                 for ( int r=0; r<input->dofs_np; r++ )
    //                 {
    //                     qtf_values[idx0+r] += 0.25 * input->water_density * int_mod * panel_k->normal_vec[r];
    //                     qtf_bern[idx0+r]   += 0.25 * input->water_density * int_mod * panel_k->normal_vec[r];
    //                 }
    //             }
    //         }
    //     }
    // }
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