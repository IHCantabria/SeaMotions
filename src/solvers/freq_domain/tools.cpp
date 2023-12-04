
// Include local modules
#include "tools.hpp"

#include "../../containers/matlin_group.hpp"
#include "../../inout/input.hpp"
#include "../../math/gauss.hpp"
#include "../../math/math_interface.hpp"


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


void    define_gauss_points_diffrac_panels(
                                                Input*      input,
                                                MeshGroup*  mesh_gp,
                                                MLGCmpx*    mat_gp
                                            )
{
    int         _count_pot_np       = 0;
    cusfloat    field_point_i[3]    = { 0.0, 0.0, 0.0 };
    GaussPoints gp( input->gauss_order );
    for ( int i=0; i<mesh_gp->panels_tnp; i++ )
    {
        if ( mesh_gp->panels[i]->type == DIFFRAC_PANEL_CODE )
        {
            for ( int j=0; j<input->gauss_order; j++ )
            {
                for ( int k=0; k<input->gauss_order; k++ )
                {
                    mesh_gp->panels[i]->local_to_global(
                                                            gp.roots[j],
                                                            gp.roots[k],
                                                            field_point_i
                                                        );
                    copy_vector( 3, field_point_i, &(mat_gp->field_points[3*_count_pot_np]) );
                    _count_pot_np++;
                }
            }
        }
    }

    for ( int i=0; i<mesh_gp->meshes_np+1; i++ )
    {
        mat_gp->field_points_cnp[i] = pow2s( input->gauss_order ) * mesh_gp->panels_raddif_cnp[i];
    }
}


void    define_gauss_points_wl(
                                                Input*      input,
                                                MeshGroup*  mesh_gp,
                                                MLGCmpx*    mat_gp
                                )
{
    // Get all WL line centers to be more accesible through a vector
    int         count_lines = 0;
    cusfloat    tv[3];      clear_vector( 3, tv );
    GaussPoints gp( input->gauss_order );

    for ( int i=0; i<mesh_gp->panels_wl_tnp; i++ )
    {
        // Create direction vector
        tv[0]   = ( mesh_gp->panels_wl[i]->x_wl[1] - mesh_gp->panels_wl[i]->x_wl[0] ) / mesh_gp->panels_wl[i]->len_wl;
        tv[1]   = ( mesh_gp->panels_wl[i]->y_wl[1] - mesh_gp->panels_wl[i]->y_wl[0] ) / mesh_gp->panels_wl[i]->len_wl;

        for ( int j=0; j<input->gauss_order; j++ )
        {
            mat_gp->field_points[3*count_lines]     = tv[0] * ( gp.roots[j] + 1.0 ) / 2.0 * mesh_gp->panels_wl[i]->len_wl + mesh_gp->panels_wl[i]->x_wl[0];
            mat_gp->field_points[3*count_lines+1]   = tv[1] * ( gp.roots[j] + 1.0 ) / 2.0 * mesh_gp->panels_wl[i]->len_wl + mesh_gp->panels_wl[i]->y_wl[0];
            mat_gp->field_points[3*count_lines+2]   = 0.0;

            count_lines++;
        }
    }

    for ( int i=0; i<mesh_gp->meshes_np+1; i++ )
    {
        mat_gp->field_points_cnp[i] = input->gauss_np_factor_1d( ) * mesh_gp->panels_wl_cnp[i];
    }

    // Get the radius from the WL line center to the body COG
    int ngpf  = input->gauss_np_factor_1d( );
    int index = 0;
    for ( int i=0; i<mesh_gp->meshes_np; i++ )
    {
        for ( int j=mat_gp->field_points_cnp[i]; j<mat_gp->field_points_cnp[i+1]; j++ )
        {
            sv_sub(
                        3,
                        &(mat_gp->field_points[3*j]),
                        input->bodies[i]->cog,
                        &(mat_gp->cog_to_field_points[3*j])
                    );
        }
    }
}