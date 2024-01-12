
// Include local modules
#include "../../math/math_interface.hpp"
#include "../../math/integration.hpp"
#include "../../waves/wave_dispersion_base_fo.hpp"


template<typename T>
void    calculate_influence_field_mat(
                                                Input*          input,
                                                MeshGroup*      mesh_gp,
                                                T*              field_funct,
                                                MLGCmpx*        field_gp
                                    )
{
    // Define potential funcions objects interface
    auto    wave_fcn    =   [field_funct]
                            (cusfloat xi, cusfloat eta, cusfloat x, cusfloat y, cusfloat z) -> cuscomplex
                            {
                                return (*field_funct)( xi, eta, x, y, z );
                            };

    // Generate potential matrix
    int         count       = 0;
    bool        is_diffrac  = false;
    cuscomplex  field_wave_term( 0.0, 0.0 );
    for ( int i=0; i<field_gp->field_points_np; i++ )
    {
        // Change field point
        field_funct->set_field_point(
                                        &(field_gp->field_points[3*i])
                                    );

        for ( int j=field_gp->start_col; j<field_gp->end_col; j++ )
        {
            // Change source point
            field_funct->set_source_i(
                                        mesh_gp->source_nodes[j],
                                        1.0
                                    );

            // Check if its necessary to discard field row
            is_diffrac  =   true;
            if ( field_gp->is_sysmat_field )
            {
                is_diffrac  = ( mesh_gp->panels[i]->type == DIFFRAC_PANEL_CODE );
            }
            
            // Compute steady and wave terms over the panel
            if ( 
                    is_diffrac
                    &&
                    mesh_gp->panels[j]->type == DIFFRAC_PANEL_CODE
                )
            {
                field_wave_term     = adaptive_quadrature_panel(
                                                                    mesh_gp->panels[j],
                                                                    wave_fcn,
                                                                    input->pot_abs_err,
                                                                    input->pot_rel_err,
                                                                    input->is_block_adaption,
                                                                    false,
                                                                    input->gauss_order
                                                                );

                field_gp->sysmat[count] = field_gp->sysmat_steady[count] + field_wave_term / 4.0 / PI;

            }
            else
            {
                field_gp->sysmat[count] = cuscomplex( 0.0, 0.0 );
            }
            
            count++;
        }
    }
}


template<typename T, typename V>
void    calculate_fields_lin(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                T*              field_func,
                                                V*              fk_fcn,
                                                cusfloat        ang_freq,
                                                cuscomplex*     intensities,
                                                cuscomplex*     raos,
                                                MLGCmpx*        field_gp,
                                                cuscomplex*     field_fk_p0,
                                                cuscomplex*     field_raddif_p0,
                                                cuscomplex*     field_total
                                    )
{
    // Declare auxiliar variables to use in the function
    int index       =   0;

    /***************************************************************/
    /************* Radiation-Diffraction Potential *****************/
    /***************************************************************/

    // Calculate potential influence coeffcients matrix
    calculate_influence_field_mat(
                                    input,
                                    mesh_gp,
                                    field_func,
                                    field_gp
                                );

    // Calculate radiation-diffraction panels potential
    calculate_fields_raddif_lin(
                                    input,
                                    intensities,
                                    field_gp
                                );

    /***************************************************************/
    /******** Sum panel potentials from all processes **************/
    /***************************************************************/
    MPI_Reduce(
                    field_gp->field_values,
                    field_raddif_p0,
                    field_gp->fields_np * field_gp->sysmat_nrows,
                    mpi_cuscomplex,
                    MPI_SUM,
                    mpi_config->proc_root,
                    MPI_COMM_WORLD
                );

    /***************************************************************/
    /************** Calculate Inicident Wave Field *****************/
    /***************************************************************/

    if ( mpi_config->is_root( ) )
    {
        // Calculate incident wave potential
        cusfloat    k   =   w2k( 
                                    ang_freq,
                                    input->water_depth,
                                    input->grav_acc
                                );
        
        for ( int i=0; i<input->heads_np; i++ )
        {
            for ( int j=0; j<field_gp->field_points_np; j++ )
            {
                index                   =   field_gp->field_points_np * i + j;
                field_fk_p0[index]      =   fk_fcn(
                                                        input->wave_amplitude,
                                                        ang_freq,
                                                        k,
                                                        input->water_depth,
                                                        input->grav_acc,
                                                        field_gp->field_points[3*j],
                                                        field_gp->field_points[3*j+1],
                                                        field_gp->field_points[3*j+2],
                                                        input->heads[i]
                                                    );
            }
        }
    }

    /***************************************************************/
    /****************** Compose total potential ********************/
    /***************************************************************/

    if ( mpi_config->is_root( ) )
    {
        calculate_total_field(
                                    input,
                                    ang_freq,
                                    field_gp,
                                    raos,
                                    field_fk_p0,
                                    field_raddif_p0,
                                    field_total
                                );
    }

}


void    calculate_fields_raddif_lin(
                                                Input*          input,
                                                cuscomplex*     intensities,
                                                MLGCmpx*        field_gp
                                    )
{
    // Loop over RHS to compute all the panel potentials the panels potentials
    cuscomplex  alpha( 1.0, 0.0 );
    cuscomplex  beta( 0.0, 0.0 );
    int         icnx = 1;
    int         icny = 1;
    for ( int i=0; i<( input->dofs_np + input->heads_np ); i++ )
    {
        cblas_gemv<cuscomplex>( 
                                    CblasRowMajor,
                                    CblasNoTrans,
                                    field_gp->sysmat_nrows,
                                    field_gp->sysmat_ncols,
                                    &alpha,
                                    field_gp->sysmat,
                                    field_gp->sysmat_ncols,
                                    &(intensities[i*field_gp->sysmat_ncols+field_gp->start_col]),
                                    icnx,
                                    &beta,
                                    &(field_gp->field_values[i*field_gp->sysmat_nrows]),
                                    icny
                                );
    }
}


void    calculate_pertubation_field(
                                                Input*          input,
                                                cuscomplex*     raos,
                                                cusfloat        ang_freq,
                                                MLGCmpx*        field_gp,
                                                cuscomplex*     field_raddif,
                                                cuscomplex*     field_total
                                        )
{
    // Declare local variables
    int index   = 0;
    int index_2 = 0;

    // Add diffraction forces
    for ( int i=0; i<input->heads_np; i++ )
    {
        for ( int j=0; j<field_gp->field_points_np; j++ )
        {
            index                   = field_gp->field_points_np * i + j;
            index_2                 = ( input->dofs_np + i ) * field_gp->field_points_np + j;
            field_total[index]      += field_raddif[index_2];
        }
    }

    // Calculate radiation field
    calculate_raddiation_field(
                                    input,
                                    raos,
                                    ang_freq,
                                    field_gp,
                                    field_raddif,
                                    field_total
                                );
}


void    calculate_raddiation_field(
                                                Input*          input,
                                                cuscomplex*     raos,
                                                cusfloat        ang_freq,
                                                MLGCmpx*        field_gp,
                                                cuscomplex*     field_raddif,
                                                cuscomplex*     field_total
                                        )
{
    // Declare local variables
    int index   = 0;
    int index_2 = 0;
    int index_3 = 0;

    // Add radiation forces
    for ( int i=0; i<input->heads_np; i++ )
    {
        for ( int j=0; j<field_gp->field_points_nb; j++ )
        {
            for ( int k=0; k<input->dofs_np; k++ )
            {
                for ( int r=field_gp->field_points_cnp[j]; r<field_gp->field_points_cnp[j+1]; r++ )
                {
                    index               = i * field_gp->field_points_np + r;
                    index_2             = k * field_gp->field_points_np + r;
                    index_3             = i * ( input->dofs_np * field_gp->field_points_nb ) + j * input->dofs_np + k;
                    field_total[index]  += cuscomplex( 0.0, -1 ) * ang_freq * raos[index_3] * field_raddif[index_2];
                }
            }
        }
    }
}


void    calculate_total_field(
                                                Input*          input,
                                                cusfloat        ang_freq,
                                                MLGCmpx*        field_gp,
                                                cuscomplex*     raos,
                                                cuscomplex*     field_fk,
                                                cuscomplex*     field_raddif,
                                                cuscomplex*     field_total
                                    )
{
    // Declare local variables
    int index   = 0;

    // Add Froude-Krylov potential
    for ( int i=0; i<input->heads_np; i++ )
    {
        for ( int j=0; j<field_gp->field_points_np; j++ )
        {
            index               =   field_gp->field_points_np * i + j;
            field_total[index]  =   field_fk[index];
        }
    }

    // Calculate pertubation field    
    calculate_pertubation_field(
                                    input,
                                    raos,
                                    ang_freq,
                                    field_gp,
                                    field_raddif,
                                    field_total
                                );
}