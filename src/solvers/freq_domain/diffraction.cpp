
// Include local modules
#include "diffraction.hpp"
#include "../../math/integration.hpp"


void    calculate_diffraction_forces_lin(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cuscomplex*     panel_pot,
                                                cusfloat        w,
                                                cuscomplex*     wave_diffrac
                                        )
{
    // Define local variables
    int         dofs_np         = input->dofs_np;
    int         elem_end_pos    = 0;
    int         elem_start_pos  = 0;
    int         index           = 0;
    int         index_1         = 0;
    cuscomplex  press_i         = cuscomplex( 0.0, 0.0 );
    cusfloat    rho_w           = input->water_density;
    
    // Loop over headings to define the wave diffraction forces
    for ( int ih=0; ih<input->heads_np; ih++ )
    {
        // Get hydromechanic coefficients for all the degrees of freedom
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

            for ( int id=0; id<dofs_np; id++ )
            {
                index = (
                            ih * ( dofs_np * mesh_gp->meshes_np )
                            +
                            ib * dofs_np
                            + 
                            id
                        );
                wave_diffrac[index]   = 0.0;
                for ( int ie=elem_start_pos; ie<elem_end_pos; ie++ )
                {
                    // Integrate pressure over panel
                    index_1             = (
                                                input->dofs_np * mesh_gp->panels_tnp 
                                                +
                                                ih * mesh_gp->panels_tnp
                                                +
                                                mesh_gp->panels_cnp[ib]
                                                +
                                                ie
                                            );
                    press_i             = (
                                                panel_pot[index_1]
                                                *
                                                mesh_gp->panels[mesh_gp->panels_cnp[ib]+ie]->normal_vec[id]
                                                *
                                                mesh_gp->panels[mesh_gp->panels_cnp[ib]+ie]->area
                                            );
                    wave_diffrac[index] += cuscomplex( 0.0, -w * rho_w ) * press_i;
                }
            }
        }
    }
}


void    calculate_diffraction_forces_nlin(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                HMFInterface*   hmf_interf,
                                                cusfloat        w,
                                                cuscomplex*     wave_diffrac
                                        )
{
    // Define local variables
    int         dofs_np = input->dofs_np;
    cusfloat    rho_w   = input->water_density;

    // Generate lambda function for the integration
    auto target_fcn = [hmf_interf]
                    (
                        cusfloat xi, 
                        cusfloat eta, 
                        cusfloat x,
                        cusfloat y,
                        cusfloat z
                    )
                    {
                        return (*hmf_interf)( xi, eta, x, y, z );
                    };

    // Allocate space for pressure vector
    int max_panels  = 0;
    for ( int i=0; i<mesh_gp->meshes_np; i++ )
    {
        // Get ith mesh panels
        if ( mesh_gp->panels_np[i] > max_panels )
        {
            max_panels = mesh_gp->panels_np[i];
        }
    }
    cuscomplex* pressure    = generate_empty_vector<cuscomplex>( mesh_gp->meshes_np * max_panels );

    // Loop over first dimension of degrees of freedrom
    int         elem_end_pos    = 0;
    int         elem_start_pos  = 0;
    int         index           = 0;
    int         index_1         = 0;
    cuscomplex  press_i         = 0.0;
    // cuscomplex  pressure = 0.0;
    for ( int ih=0; ih<input->heads_np; ih++ )
    {
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

            // Set start index for the sources evaluation
            hmf_interf->set_start_index_i( 
                                                mesh_gp->source_nodes_tnp * ( dofs_np + ih ),
                                                0,
                                                mesh_gp->source_nodes_tnp
                                            );

            // Loop over panels to integrate the wave radiation
            // pressure along the floating object external shape
            for ( int ie=elem_start_pos; ie<elem_end_pos; ie++ )
            {
                // Set new panel
                hmf_interf->set_panel( mesh_gp->panels[mesh_gp->panels_cnp[ib]+ie] );

                // Integrate pressure over panel
                index           = max_panels * ib + ie;
                if ( mesh_gp->panels[mesh_gp->panels_cnp[ib]+ie]->type == DIFFRAC_PANEL_CODE )
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

        // Get hydromechanic coefficients for all the degrees of freedom
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

            for ( int id=0; id<dofs_np; id++ )
            {
                index = (
                            ih * ( dofs_np * mesh_gp->meshes_np )
                            +
                            ib * dofs_np
                            + 
                            id
                        );

                wave_diffrac[index]   = 0.0;
                for ( int ie=elem_start_pos; ie<elem_end_pos; ie++ )
                {
                    // Integrate pressure over panel
                    index_1             = max_panels * ib + ie;
                    press_i             = pressure[index_1] * mesh_gp->panels[mesh_gp->panels_cnp[ib]+ie]->normal_vec[id];
                    wave_diffrac[index] += cuscomplex( 0.0, -w * rho_w ) * press_i;
                }
            }
        }
    }
}