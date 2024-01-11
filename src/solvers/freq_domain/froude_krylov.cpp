
// Include local modules
#include "../../math/integration.hpp"
#include "froude_krylov.hpp"
#include "../../waves/waves_common.hpp"


void    calculate_froude_krylov(
                                    Input*          input,
                                    MpiConfig*      mpi_config,
                                    MeshGroup*      mesh_gp,
                                    cusfloat        ang_freq,
                                    cuscomplex*     froude_krylov
                                )
{
    // Calculate wave number
    cusfloat k      = w2k( ang_freq, input->water_depth, input->grav_acc );
    cusfloat rho_w  = input->water_density;

    // Generate pressure vector to storage the pressures over the panels
    int max_panels  = 0;
    for ( int i=0; i<mesh_gp->meshes_np; i++ )
    {
        // Get ith mesh panels
        if ( mesh_gp->panels_np[i] > max_panels )
        {
            max_panels = mesh_gp->panels_np[i];
        }
    }
    cuscomplex* pressure = generate_empty_vector<cuscomplex>( max_panels * mesh_gp->meshes_np );

    // Loop around headings to get the Froude-Krylov force
    // for each of them
    int         elem_end_pos    = 0;
    int         elem_start_pos  = 0;
    int         index           = 0;
    int         index_1         = 0;
    cuscomplex  press_i         = cuscomplex( 0.0, 0.0 );

    for ( int ih=0; ih<input->heads_np; ih++ )
    {
        // Definition of the target function to integrate incident
        // wave potential for the desised heading
        auto    target_fcn  =   [ 
                                    input,
                                    ang_freq, 
                                    k,
                                    ih
                                ]
                                (
                                    cusfloat ,
                                    cusfloat ,
                                    cusfloat x,
                                    cusfloat y,
                                    cusfloat z
                                )
                                {
                                    return wave_potential_airy_space(
                                                                        1.0,
                                                                        ang_freq,
                                                                        k,
                                                                        input->water_depth,
                                                                        input->grav_acc,
                                                                        x,
                                                                        y,
                                                                        z,
                                                                        input->heads[ih]
                                                                    );
                                };
        
        // Loop around bodies to calculate the pressure of the incident wave
        // over each panel of the mesh
        for ( int ib=0; ib<mesh_gp->meshes_np; ib++ )
        {
            // Calculate MPI data chunks
            elem_end_pos    = 0;
            elem_start_pos  = 0;
            mpi_config->get_1d_bounds( 
                                            mesh_gp->panels_np[ib], 
                                            elem_start_pos, 
                                            elem_end_pos 
                                        );

            // Loop over panels integrating incident wave pressure
            for ( int ie=elem_start_pos; ie<elem_end_pos; ie++ )
            {
                // Integrate pressure over panel
                index           = max_panels * ib + ie;
                if ( mesh_gp->panels[mesh_gp->panels_cnp[ib]+ie]->type == 0 )
                {
                    pressure[index] = adaptive_quadrature_panel(
                                                                    mesh_gp->panels[mesh_gp->panels_cnp[ib]+ie],
                                                                    target_fcn,
                                                                    input->press_abs_err,
                                                                    input->press_rel_err,
                                                                    input->is_block_adaption,
                                                                    true,
                                                                    input->gauss_order
                                                                );
                }
            }
        }

        // Loop over bodies to calculate the projection of the Froude-Krylov force
        // over each DOF
        for ( int ib=0; ib<mesh_gp->meshes_np; ib++ )
        {
            // Calculate MPI data chunks
            elem_end_pos    = 0;
            elem_start_pos  = 0;
            mpi_config->get_1d_bounds( 
                                            mesh_gp->panels_np[ib], 
                                            elem_start_pos, 
                                            elem_end_pos 
                                        );

            for ( int id=0; id<input->dofs_np; id++ )
            {
                // Generate index for the current Froude-Krylov force and
                // to clean the memory space
                index                   = (
                                                ih * ( mesh_gp->meshes_np * input->dofs_np )
                                                +
                                                ib * input->dofs_np
                                                +
                                                id
                                            );
                
                froude_krylov[index]    = 0.0;
                // Loop over panels to sum all the pressures
                for ( int ie=elem_start_pos; ie<elem_end_pos; ie++ )
                {
                    // Integrate pressure over panel
                    index_1                 = max_panels * ib + ie;
                    press_i                 = pressure[index_1] * mesh_gp->panels[mesh_gp->panels_cnp[ib]+ie]->normal_vec[id];
                    froude_krylov[index]    += cuscomplex( 0.0, -ang_freq * rho_w ) * press_i;
                }
            }
        }
    }

    // Deallocate local heap memory
    mkl_free( pressure );
}
