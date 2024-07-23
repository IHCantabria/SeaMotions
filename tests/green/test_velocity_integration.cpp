
// Include local modules
#include "../../src/containers/panel_geom.hpp"
#include "../../src/interfaces/gwfdx_interface.hpp"
#include "../../src/interfaces/gwfdy_interface.hpp"
#include "../../src/interfaces/gwfdz_interface.hpp"
#include "../../src/green/source.hpp"
#include "../../src/math/gauss.hpp"
#include "../../src/math/math_tools.hpp"
#include "../../src/math/integration.hpp"
#include "../../src/mesh/mesh.hpp"
#include "../../src/mesh/mesh_group.hpp"
#include "../../src/tools.hpp"


struct _Input
{
public:
    // Define class attributes
    cusfloat    ang_freq            = 0.0;
    int         bodies_np           = 0;
    int         gauss_order         = 0;
    cusfloat    gfdn_abs_err        = 0.0;
    cusfloat    gfdn_rel_err        = 0.0;
    cusfloat    grav_acc            = 9.81;
    bool        is_block_adaption   = false;
    cusfloat    period              = 0.0;
    cusfloat    water_depth         = 0.0;

};

// Declare local module functions
void    velocity_integral(
                                            _Input*     input,
                                            MeshGroup*  mesh_gp,
                                            PanelGeom*  panel,
                                            cuscomplex* vel_x,
                                            cuscomplex* vel_y,
                                            cuscomplex* vel_z
                        );


void    velocity_integral_convergence(
                                            _Input*     input,
                                            MeshGroup*  mesh_gp,
                                            PanelGeom*  panel
                                    );


void    velocity_steady(
                                            _Input*     input,
                                            PanelGeom*  panel,
                                            PanelGeom*  panel_mirror,
                                            cusfloat*   field_point,
                                            cuscomplex* vel_x,
                                            cuscomplex* vel_y,
                                            cuscomplex* vel_z
                        );


void    launch_test_1( 
                        _Input*     input,
                        MeshGroup*  mesh_gp
                    )
{
    // Define integration panel
    PanelGeom*  panel = new PanelGeom;

    panel->x[0] = -0.5;
    panel->x[1] = 0.5;
    panel->x[2] = 0.5;
    panel->x[3] = -0.5;

    panel->y[0] = -0.5;
    panel->y[1] = -0.5;
    panel->y[2] = 0.5;
    panel->y[3] = 0.5;

    panel->z[0] = -9.5;
    panel->z[1] = -9.5;
    panel->z[2] = -9.5;
    panel->z[3] = -9.5;

    panel->num_nodes = 4;

    cusfloat cog[3] = { 0.0, 0.0, 0.0 };
    panel->calculate_properties( cog );

    // Launch convergence study
    velocity_integral_convergence( 
                                        input,
                                        mesh_gp,
                                        panel
                                    );

    // Delete heap memory allocated objects
    delete panel;
}


void    launch_test_2( 
                        _Input*     input,
                        MeshGroup*  mesh_gp
                    )
{
    // Define integration panel
    PanelGeom*  panel = new PanelGeom;

    panel->x[0] = -5.0;
    panel->x[1] = -3.0;
    panel->x[2] = -3.0;
    panel->x[3] = -5.0;

    panel->y[0] = -5.0;
    panel->y[1] = -5.0;
    panel->y[2] = -3.0;
    panel->y[3] = -3.0;

    panel->z[0] = -10.0;
    panel->z[1] = -10.0;
    panel->z[2] = -10.0;
    panel->z[3] = -10.0;

    panel->num_nodes = 4;

    cusfloat cog[3] = { 0.0, 0.0, 0.0 };
    panel->calculate_properties( cog );

    // Launch convergence study
    velocity_integral_convergence( 
                                        input,
                                        mesh_gp,
                                        panel
                                    );

    // Delete heap memory allocated objects
    delete panel;
}


int     main(int argc, char* argv[])
{
    // Read command line arguments
    if (!check_num_cmd_args(argc, 1))
    {
        return 1;
    }

    std::string mesh_fipath(argv[1]);

    // Define problem inputs´
    _Input* input               = new _Input;

    input->period               = 7.0;
    input->ang_freq             = 2.0 * PI / input->period;
    input->bodies_np            = 1;
    input->gauss_order          = 1;
    input->gfdn_abs_err         = 1e-3;
    input->gfdn_rel_err         = 1e-3;
    input->grav_acc             = 9.81;
    input->is_block_adaption    = true;
    input->water_depth          = 1000.0;

    // Generate mesh of panels to calculate influence velocities
    cusfloat    cog[3]  = { 0.0, 0.0, 0.0 };
    Mesh*   msh         = new   Mesh( 
                                        mesh_fipath,
                                        std::string( "box" ),
                                        cog,
                                        DIFFRAC_PANEL_CODE                                
                                    );
    msh->define_source_nodes(
                                0,
                                cog
                            );

    Mesh** all_meshes   = new Mesh*[input->bodies_np];
    for ( int i=0; i<input->bodies_np; i++ )
    {
        all_meshes[i] = msh;
    }

    MeshGroup*  mesh_gp =   new   MeshGroup( 
                                                all_meshes,
                                                input->bodies_np,
                                                false
                                            );

    mesh_gp->define_mirror_panels( );

    // Launch test 1
    launch_test_1( 
                        input,
                        mesh_gp
                    );

    // Launch test 2
    launch_test_2( 
                        input,
                        mesh_gp
                    );
    
    // Delete heap memory
    delete input;
    delete mesh_gp;
    delete msh;

    return 0;
}


void    velocity_integral(
                                        _Input*     input,
                                        MeshGroup*  mesh_gp,
                                        PanelGeom*  panel,
                                        cuscomplex* vel_x,
                                        cuscomplex* vel_y,
                                        cuscomplex* vel_z,
                                        cuscomplex* vel_mod
                        )
{
    GWFDxInterface* green_wave_dx   = new   GWFDxInterface( 
                                                                mesh_gp->source_nodes[0],
                                                                mesh_gp->source_nodes[0],
                                                                input->ang_freq,
                                                                input->water_depth,
                                                                input->grav_acc
                                                            );
    
    GWFDyInterface* green_wave_dy   = new   GWFDyInterface( 
                                                                mesh_gp->source_nodes[0],
                                                                mesh_gp->source_nodes[0],
                                                                input->ang_freq,
                                                                input->water_depth,
                                                                input->grav_acc
                                                            );

    GWFDzInterface* green_wave_dz   = new   GWFDzInterface( 
                                                                mesh_gp->source_nodes[0],
                                                                mesh_gp->source_nodes[0],
                                                                input->ang_freq,
                                                                input->water_depth,
                                                                input->grav_acc
                                                            );

    auto            wave_dx_fcn     =   [green_wave_dx]
                                        (cusfloat xi, cusfloat eta, cusfloat x, cusfloat y, cusfloat z) -> cuscomplex
                                        {
                                            return (*green_wave_dx)( xi, eta, x, y, z );
                                        };
    
    auto            wave_dy_fcn     =   [green_wave_dy]
                                        (cusfloat xi, cusfloat eta, cusfloat x, cusfloat y, cusfloat z) -> cuscomplex
                                        {
                                            return (*green_wave_dy)( xi, eta, x, y, z );
                                        };

    auto            wave_dz_fcn     =   [green_wave_dz]
                                        (cusfloat xi, cusfloat eta, cusfloat x, cusfloat y, cusfloat z) -> cuscomplex
                                        {
                                            return (*green_wave_dz)( xi, eta, x, y, z );
                                        };
    
    // Clear input variables to avoid using previous data
    (*vel_x)    = cuscomplex( 0.0, 0.0 );
    (*vel_y)    = cuscomplex( 0.0, 0.0 );
    (*vel_z)    = cuscomplex( 0.0, 0.0 );
    (*vel_mod)  = cuscomplex( 0.0, 0.0 );

    // Define auxiliar variables to get ith velocity field
    GaussPoints gp( input->gauss_order );
    cusfloat    field_point_i[3] = { 0.0, 0.0, 0.0 };
    cuscomplex  vel_x_i( 0.0, 0.0 );    cuscomplex  vel_x_j( 0.0, 0.0 );
    cuscomplex  vel_y_i( 0.0, 0.0 );    cuscomplex  vel_y_j( 0.0, 0.0 );
    cuscomplex  vel_z_i( 0.0, 0.0 );    cuscomplex  vel_z_j( 0.0, 0.0 );

    for ( int gpi=0; gpi<input->gauss_order; gpi++ )
    {
        for ( int gpj=0; gpj<input->gauss_order; gpj++ )
        {
            // Calculate new field point
            panel->local_to_global(
                                        gp.roots[gpi],
                                        gp.roots[gpj],
                                        field_point_i
                                    );
            
            // Change field point
            green_wave_dx->set_field_point(
                                                    field_point_i
                                                );
            green_wave_dy->set_field_point(
                                                    field_point_i
                                                );
            green_wave_dz->set_field_point(
                                                    field_point_i
                                                );

            // Loop over the panels to integrate the velocity field
            vel_x_j = cuscomplex( 0.0, 0.0 );
            vel_y_j = cuscomplex( 0.0, 0.0 );
            vel_z_j = cuscomplex( 0.0, 0.0 );
            for ( int i=0; i<mesh_gp->panels_tnp; i++ )
            {
                vel_x_i = cuscomplex( 0.0, 0.0 );
                vel_y_i = cuscomplex( 0.0, 0.0 );
                vel_z_i = cuscomplex( 0.0, 0.0 );

                // Integrate steady velocity part
                velocity_steady(
                                    input,
                                    mesh_gp->panels[i],
                                    mesh_gp->panels_mirror[i],
                                    panel->center,
                                    &vel_x_i,
                                    &vel_y_i,
                                    &vel_z_i
                                );

                // Integrate wave velocity part
                vel_x_i     += adaptive_quadrature_panel(
                                                            mesh_gp->panels[i],
                                                            wave_dx_fcn,
                                                            input->gfdn_abs_err,
                                                            input->gfdn_rel_err,
                                                            input->is_block_adaption,
                                                            false,
                                                            input->gauss_order
                                                        );
                vel_y_i     += adaptive_quadrature_panel(
                                                            mesh_gp->panels[i],
                                                            wave_dy_fcn,
                                                            input->gfdn_abs_err,
                                                            input->gfdn_rel_err,
                                                            input->is_block_adaption,
                                                            false,
                                                            input->gauss_order
                                                        );
                vel_z_i     += adaptive_quadrature_panel(
                                                            mesh_gp->panels[i],
                                                            wave_dz_fcn,
                                                            input->gfdn_abs_err,
                                                            input->gfdn_rel_err,
                                                            input->is_block_adaption,
                                                            false,
                                                            input->gauss_order
                                                        );
                
                // Add velocities from each panel to the total contribution
                vel_x_j     +=  vel_x_i;
                vel_y_j     +=  vel_y_i;
                vel_z_j     +=  vel_z_i;
            }
            std::cout << "GPI: " << gpi << " - GPJ: " << gpj << " - vel_x: " << vel_x_j << " - vel_y: " << vel_y_j << " - vel_z: " << vel_z_j << std::endl;
            (*vel_x)    +=  vel_x_j * gp.weights[gpi] * gp.weights[gpj];
            (*vel_y)    +=  vel_y_j * gp.weights[gpi] * gp.weights[gpj];
            (*vel_z)    +=  vel_z_j * gp.weights[gpi] * gp.weights[gpj];
            (*vel_mod)  +=  (
                                vel_x_j * std::conj( vel_x_j )
                                +
                                vel_y_j * std::conj( vel_y_j )
                                +
                                vel_z_j * std::conj( vel_z_j )
                            ) * gp.weights[gpi] * gp.weights[gpj];
        }
    }

    // Delete heap memory objects
    delete green_wave_dx;
    delete green_wave_dy;
    delete green_wave_dz;
}


void    velocity_integral_convergence(
                                            _Input*     input,
                                            MeshGroup*  mesh_gp,
                                            PanelGeom*  panel
                                    )
{
    // Define variables to storage integral velocity over the panel
    const   int gp_np           = 4;
            int gp[gp_np]       = { 1, 2, 3, 4 };
    cuscomplex  vel_x[gp_np];   clear_vector( gp_np, vel_x );
    cuscomplex  vel_y[gp_np];   clear_vector( gp_np, vel_y );
    cuscomplex  vel_z[gp_np];   clear_vector( gp_np, vel_z );
    cuscomplex  vel_mod[gp_np]; clear_vector( gp_np, vel_mod );

    // Loop over Gauss Orders
    for ( int i=0; i<gp_np; i++ )
    {
        // Change Gauss integration order
        input->gauss_order = gp[i];

        // Calculate velocity integral
        velocity_integral(
                            input,
                            mesh_gp,
                            panel,
                            &(vel_x[i]),
                            &(vel_y[i]),
                            &(vel_z[i]),
                            &(vel_mod[i])
                        );
    }

    // Print out results
    for ( int i=0; i<gp_np; i++ )
    {
        std::cout << "GP: " << i;
        std::cout << " - Vel.X:     " << vel_x[i];
        std::cout << " - Vel.Y:     " << vel_y[i];
        std::cout << " - Vel.Z:     " << vel_z[i];
        std::cout << " - Press.Dyn: " << vel_mod[i];
        std::cout << std::endl;
    }
}


void    velocity_steady(
                                            _Input*     input,
                                            PanelGeom*  panel,
                                            PanelGeom*  panel_mirror,
                                            cusfloat*   field_point,
                                            cuscomplex* vel_x,
                                            cuscomplex* vel_y,
                                            cuscomplex* vel_z
                        )
{
    // Reset velocity values
    cusfloat    field_point_i[3];   clear_vector( 3, field_point_i );
    cusfloat    vel_0[3];           clear_vector( 3, vel_0 );
    cusfloat    vel_1[3];           clear_vector( 3, vel_1 );
    cusfloat    vel_2[3];           clear_vector( 3, vel_2 );
    cusfloat    vel_3[3];           clear_vector( 3, vel_3 );
    cusfloat    vel_4[3];           clear_vector( 3, vel_4 );
    cusfloat    vel_5[3];           clear_vector( 3, vel_5 );
    
    // Calcualte velocity corresponding to the r0 source
    calculate_source_velocity_newman(
                                        panel,
                                        field_point, 
                                        0,
                                        0, 
                                        vel_0
                                    );

    // Calculate velocity corresponding to the r1 source
    field_point_i[0]    =   field_point[0];
    field_point_i[1]    =   field_point[1];
    field_point_i[2]    =   field_point[2] + 2 * input->water_depth;
    calculate_source_velocity_newman(
                                        panel_mirror,
                                        field_point_i, 
                                        0,
                                        0, 
                                        vel_1
                                    );
    
    // Calculate velocity corresponding to the r2 source
    field_point_i[0]    =   field_point[0];
    field_point_i[1]    =   field_point[1];
    field_point_i[2]    =   field_point[2];
    calculate_source_velocity_newman(
                                        panel_mirror,
                                        field_point_i, 
                                        0,
                                        0, 
                                        vel_2
                                    );

    // Calculate velocity corresponding to the r3 source
    field_point_i[0]    =   field_point[0];
    field_point_i[1]    =   field_point[1];
    field_point_i[2]    =   field_point[2] + 2.0 * input->water_depth;
    calculate_source_velocity_newman(
                                        panel,
                                        field_point_i, 
                                        0,
                                        0, 
                                        vel_3
                                    );

    // Calculate velocity corresponding to the r4 source
    field_point_i[0]    =   field_point[0];
    field_point_i[1]    =   field_point[1];
    field_point_i[2]    =   -field_point[2] + 2.0 * input->water_depth;
    calculate_source_velocity_newman(
                                        panel_mirror,
                                        field_point_i, 
                                        0,
                                        0, 
                                        vel_4
                                    );

    // Calculate velocity corresponding to the r5 source
    field_point_i[0]    =   field_point[0];
    field_point_i[1]    =   field_point[1];
    field_point_i[2]    =   field_point[2] + 4.0 * input->water_depth;
    calculate_source_velocity_newman(
                                        panel_mirror,
                                        field_point_i, 
                                        0,
                                        0, 
                                        vel_5
                                    );
    
    // Compose total velocity vector
    (*vel_x)    =  ( vel_0[0] + vel_1[0] + vel_2[0] + vel_3[0] + vel_4[0] + vel_5[0] ) / 4.0 / PI;
    (*vel_y)    =  ( vel_0[1] + vel_1[1] + vel_2[1] + vel_3[1] + vel_4[1] + vel_5[1] ) / 4.0 / PI;
    (*vel_z)    =  ( vel_0[2] + vel_1[2] + vel_2[2] + vel_3[2] + vel_4[2] + vel_5[2] ) / 4.0 / PI;
}