
// Include general usage libraries
#include <fstream>

// Include local modules
#include "tools.hpp"

#include "../../containers/matlin_group.hpp"
#include "../../inout/input.hpp"
#include "../../math/gauss.hpp"


void        calculate_field_point_rot(
                                                cuscomplex*     raos_trans,
                                                cuscomplex*     raos_rot,
                                                cusfloat*       field_point,
                                                cusfloat*       cog,
                                                cuscomplex*     point_disp
                                    )
{
    // Define vector from cog to field point
    cuscomplex cog_to_fp_c[3];      clear_vector( 3, cog_to_fp_c );
    for ( int r=0; r<3; r++ )
    {
        cog_to_fp_c[r]  = cuscomplex( field_point[r] - cog[r], 0.0 );
    }

    // Calculate first order displacement of the panel centre
    clear_vector( 3, point_disp );

    cross(
                raos_rot,
                cog_to_fp_c,
                point_disp
        );
    sv_add(
                3,
                point_disp,
                raos_trans,
                point_disp
            );
}


void        calculate_field_point_vel_rot(
                                                cuscomplex*     raos_trans,
                                                cuscomplex*     raos_rot,
                                                cusfloat*       field_point,
                                                cusfloat*       cog,
                                                cusfloat        ang_freq,
                                                cuscomplex*     point_disp
                                        )
{
    // Calculate point position
    calculate_field_point_rot(
                                    raos_trans,
                                    raos_rot,
                                    field_point,
                                    cog,
                                    point_disp
                                );

    // Scale with angular frequency in order to
    // obtain the point velocity
    svs_mult( 
                3,
                point_disp,
                cuscomplex( 0.0, -ang_freq ),
                point_disp
            );
}


std::string compose_dof_path( 
                                                std::string     base_path,
                                                int             dofs_num,
                                                int             ang_freq_num
                            )
{
    std::stringstream ss0;
    ss0 << base_path << "dof_" << dofs_num << "/";
    ss0 << "dof_" << dofs_num << "_ang_freq_" << ang_freq_num << ".dat";

    return ss0.str( );
}


void        define_gauss_points_diffrac_panels(
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


void        define_gauss_points_wl(
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


void        storage_radiation_potential( 
                                                std::string dof_base_path,
                                                int         field_points_np,
                                                cuscomplex* rad_potential
                                    )
{
    std::ofstream outfile( dof_base_path );
    CHECK_FILE_UNIT_STATUS( outfile, dof_base_path );

    outfile << field_points_np << " " << 1 << std::endl;
    for ( int i=0; i<field_points_np; i++ )
    {
        outfile << rad_potential[i].real( ) << " ";
        outfile << rad_potential[i].imag( ) << std::endl;
    }

    outfile.close( );
}