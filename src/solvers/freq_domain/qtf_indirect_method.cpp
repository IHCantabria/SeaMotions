
// Include local modules
#include "qtf_indirect_method.hpp"

#include "froude_krylov.hpp"


void    calculate_secord_force_indirect(
                                            Input*      input,
                                            MeshGroup*  mesh_gp,
                                            cusfloat    ang_freq_i,
                                            cusfloat    ang_freq_j,
                                            bool        is_diff,
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
                                    is_diff,
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