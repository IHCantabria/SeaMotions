
// Include general usage libraries
#include <iostream>

// Include local modules
#include "../../src/containers/panel_geom.hpp"
#include "../../src/green/source.hpp"
#include "../../src/interfaces/grfdx_interface.hpp"
#include "../../src/interfaces/grfdy_interface.hpp"
#include "../../src/interfaces/grfdz_interface.hpp"
#include "../../src/math/math_tools.hpp"
#include "../../src/math/integration.hpp"


cusfloat ABS_EPS = 1e-3;
cusfloat REL_EPS = 1e-3;


void test_0( void )
{
    // Create panel for checking 
    PanelGeom* panel = new PanelGeom;

    panel->x[0] = -1.0;
    panel->x[1] =  1.0;
    panel->x[2] =  1.0;
    panel->x[3] = -1.0;

    panel->y[0] = -1.0;
    panel->y[1] = -1.0;
    panel->y[2] =  1.0;
    panel->y[3] =  1.0;

    panel->z[0] =  0.0;
    panel->z[1] =  0.0;
    panel->z[2] =  0.0;
    panel->z[3] =  0.0;
    
    panel->num_nodes = 4;

    cusfloat cog[3] = { 0.0, 0.0, 0.0 };
    panel->calculate_properties( cog );

    // Define field points
    const int fp_np = 9;
    cusfloat field_points[3*fp_np];

    field_points[0] = 0.0;
    field_points[1] = 0.0;
    field_points[2] = 0.0;

    field_points[3] = 0.25;
    field_points[4] = 0.0;
    field_points[5] = 0.0;

    field_points[6] = 0.45;
    field_points[7] = 0.0;
    field_points[8] = 0.0;

    field_points[9]  = 0.6;
    field_points[10] = 0.0;
    field_points[11] = 0.0;

    field_points[12] = 1.0;
    field_points[13] = 0.0;
    field_points[14] = 0.0;

    field_points[15] = 1.5;
    field_points[16] = 0.0;
    field_points[17] = 0.0;

    field_points[18] = 3.0;
    field_points[19] = 0.0;
    field_points[20] = 0.0;

    field_points[21] = 6.0;
    field_points[22] = 0.0;
    field_points[23] = 0.0;

    field_points[24] = 20.0;
    field_points[25] = 0.0;
    field_points[26] = 0.0;
    
    // Allocate memory for the velocity vectors computed
    // used the three different formulations
    const int   vel_np                  = 3 * fp_np;
    cusfloat    hess_vel[3*fp_np];      clear_vector( vel_np, hess_vel );
    cusfloat    nw_vel[3*fp_np];        clear_vector( vel_np, nw_vel );
    cusfloat    quad_vel[3*fp_np];      clear_vector( vel_np, quad_vel );
    cuscomplex  vel_aux( 0.0, 0.0 );

    std::cout << "***************************************" << std::endl;
    std::cout << "***************  TEST_0 ***************" << std::endl;
    std::cout << "***************************************" << std::endl;
    for ( int i=0; i<fp_np; i++ )
    {
        // Calculate velocity vector using Hess & Smith formulation
        calculate_source_velocity_hess(
                                            panel,
                                            &(field_points[3*i]), 
                                            0,
                                            0, 
                                            &(hess_vel[3*i])
                                        );

        // Calculate velocity vector using Newman formulation
        calculate_source_velocity_newman(
                                            panel,
                                            &(field_points[3*i]), 
                                            0,
                                            0, 
                                            &(nw_vel[3*i])
                                        );

        // Define lambda functions for the derivatives of the 
        // Rankine sources
        auto rank_dx_fcn    =   [field_points, i]
                                (cusfloat xi, cusfloat eta, cusfloat x, cusfloat y, cusfloat z)
                                {
                                    cusfloat dx = field_points[3*i] - x;
                                    cusfloat dy = field_points[3*i+1] - y;
                                    cusfloat dz = field_points[3*i+2] - z;
                                    cusfloat R  = std::sqrt( pow2s( dx ) + pow2s( dy ) + pow2s( dz ) );

                                    return dx / pow3s( R );
                                };

        auto rank_dy_fcn    =   [field_points, i]
                                (cusfloat xi, cusfloat eta, cusfloat x, cusfloat y, cusfloat z)
                                {
                                    cusfloat dx = field_points[3*i] - x;
                                    cusfloat dy = field_points[3*i+1] - y;
                                    cusfloat dz = field_points[3*i+2] - z;
                                    cusfloat R  = std::sqrt( pow2s( dx ) + pow2s( dy ) + pow2s( dz ) );

                                    return dy / pow3s( R );
                                };
        
        auto rank_dz_fcn    =   [field_points, i]
                                (cusfloat xi, cusfloat eta, cusfloat x, cusfloat y, cusfloat z)
                                {
                                    cusfloat dx = field_points[3*i] - x;
                                    cusfloat dy = field_points[3*i+1] - y;
                                    cusfloat dz = field_points[3*i+2] - z;
                                    cusfloat R  = std::sqrt( pow2s( dx ) + pow2s( dy ) + pow2s( dz ) );

                                    return dz / pow3s( R );
                                };

        // vel_aux =   adaptive_quadrature_panel(
        //                                         panel,
        //                                         rank_dx_fcn,
        //                                         ABS_EPS,
        //                                         REL_EPS,
        //                                         false,
        //                                         true,
        //                                         1
        //                                     );
        // quad_vel[3*i]   = vel_aux.real( );
        
        // vel_aux =   adaptive_quadrature_panel(
        //                                         panel,
        //                                         rank_dy_fcn,
        //                                         ABS_EPS,
        //                                         REL_EPS,
        //                                         false,
        //                                         true,
        //                                         1
        //                                     );
        // quad_vel[3*i+1] = vel_aux.real( );
        
        // vel_aux =   adaptive_quadrature_panel(
        //                                         panel,
        //                                         rank_dz_fcn,
        //                                         ABS_EPS,
        //                                         REL_EPS,
        //                                         false,
        //                                         true,
        //                                         1
        //                                     );
        // quad_vel[3*i+2] = vel_aux.real( );

        // Print out results
        std::cout << std::endl;
        std::cout << "-> FIELD POINT[" << i << "]: "; print_vector( 3, &(field_points[3*i]), 0, 3 );
        std::cout << " --> Hess:    "; print_vector( 3, &(hess_vel[3*i]), 0, 6 );
        std::cout << " --> Newman:  "; print_vector( 3, &(nw_vel[3*i]), 0, 6 );
        std::cout << " --> Quad:    "; print_vector( 3, &(quad_vel[3*i]), 0, 6 );
        std::cout << " --> Is inside: " << panel->is_inside( &(field_points[3*i]) ) << std::endl;
        std::cout << std::endl;
    }
}


void test_1( void )
{
    // Create panel for checking 
    PanelGeom* panel = new PanelGeom;

    panel->x[0] = -1.0;
    panel->x[1] =  1.0;
    panel->x[2] =  1.0;
    panel->x[3] = -1.0;

    panel->y[0] =  1.0;
    panel->y[1] =  1.0;
    panel->y[2] =  1.0;
    panel->y[3] =  1.0;

    panel->z[0] = -1.0;
    panel->z[1] = -1.0;
    panel->z[2] =  1.0;
    panel->z[3] =  1.0;
    
    panel->num_nodes = 4;

    cusfloat cog[3] = { 0.0, 0.0, 0.0 };
    panel->calculate_properties( cog );

    // Define field points
    const int fp_np = 9;
    cusfloat field_points[3*fp_np];

    field_points[0] = 0.0;
    field_points[1] = 1.0;
    field_points[2] = 0.0;

    field_points[3] = 0.25;
    field_points[4] = 1.0;
    field_points[5] = 0.0;

    field_points[6] = 0.45;
    field_points[7] = 1.0;
    field_points[8] = 0.0;

    field_points[9]  = 0.6;
    field_points[10] = 1.0;
    field_points[11] = 0.0;

    field_points[12] = 1.0;
    field_points[13] = 1.0;
    field_points[14] = 0.0;

    field_points[15] = 1.5;
    field_points[16] = 1.0;
    field_points[17] = 0.0;

    field_points[18] = 3.0;
    field_points[19] = 1.0;
    field_points[20] = 0.0;

    field_points[21] = 6.0;
    field_points[22] = 1.0;
    field_points[23] = 0.0;

    field_points[24] = 20.0;
    field_points[25] = 1.0;
    field_points[26] = 0.0;
    
    // Allocate memory for the velocity vectors computed
    // used the three different formulations
    const int   vel_np                  = 3 * fp_np;
    cusfloat    hess_vel[3*fp_np];      clear_vector( vel_np, hess_vel );
    cusfloat    nw_vel[3*fp_np];        clear_vector( vel_np, nw_vel );
    cusfloat    quad_vel[3*fp_np];      clear_vector( vel_np, quad_vel );
    cuscomplex  vel_aux( 0.0, 0.0 );

    std::cout << "***************************************" << std::endl;
    std::cout << "***************  TEST_1 ***************" << std::endl;
    std::cout << "***************************************" << std::endl;
    for ( int i=0; i<fp_np; i++ )
    {
        // Calculate velocity vector using Hess & Smith formulation
        calculate_source_velocity_hess(
                                            panel,
                                            &(field_points[3*i]), 
                                            0,
                                            0, 
                                            &(hess_vel[3*i])
                                        );

        // Calculate velocity vector using Newman formulation
        calculate_source_velocity_newman(
                                            panel,
                                            &(field_points[3*i]), 
                                            0,
                                            0, 
                                            &(nw_vel[3*i])
                                        );

        // Define lambda functions for the derivatives of the 
        // Rankine sources
        auto rank_dx_fcn    =   [field_points, i]
                                (cusfloat xi, cusfloat eta, cusfloat x, cusfloat y, cusfloat z)
                                {
                                    cusfloat dx = field_points[3*i] - x;
                                    cusfloat dy = field_points[3*i+1] - y;
                                    cusfloat dz = field_points[3*i+2] - z;
                                    cusfloat R  = std::sqrt( pow2s( dx ) + pow2s( dy ) + pow2s( dz ) );

                                    return dx / pow3s( R );
                                };

        auto rank_dy_fcn    =   [field_points, i]
                                (cusfloat xi, cusfloat eta, cusfloat x, cusfloat y, cusfloat z)
                                {
                                    cusfloat dx = field_points[3*i] - x;
                                    cusfloat dy = field_points[3*i+1] - y;
                                    cusfloat dz = field_points[3*i+2] - z;
                                    cusfloat R  = std::sqrt( pow2s( dx ) + pow2s( dy ) + pow2s( dz ) );

                                    return dy / pow3s( R );
                                };
        
        auto rank_dz_fcn    =   [field_points, i]
                                (cusfloat xi, cusfloat eta, cusfloat x, cusfloat y, cusfloat z)
                                {
                                    cusfloat dx = field_points[3*i] - x;
                                    cusfloat dy = field_points[3*i+1] - y;
                                    cusfloat dz = field_points[3*i+2] - z;
                                    cusfloat R  = std::sqrt( pow2s( dx ) + pow2s( dy ) + pow2s( dz ) );

                                    return dz / pow3s( R );
                                };

        // vel_aux =   adaptive_quadrature_panel(
        //                                         panel,
        //                                         rank_dx_fcn,
        //                                         ABS_EPS,
        //                                         REL_EPS,
        //                                         false,
        //                                         true,
        //                                         1
        //                                     );
        // quad_vel[3*i]   = vel_aux.real( );
        
        // vel_aux =   adaptive_quadrature_panel(
        //                                         panel,
        //                                         rank_dy_fcn,
        //                                         ABS_EPS,
        //                                         REL_EPS,
        //                                         false,
        //                                         true,
        //                                         1
        //                                     );
        // quad_vel[3*i+1] = vel_aux.real( );
        
        // vel_aux =   adaptive_quadrature_panel(
        //                                         panel,
        //                                         rank_dz_fcn,
        //                                         ABS_EPS,
        //                                         REL_EPS,
        //                                         false,
        //                                         true,
        //                                         1
        //                                     );
        // quad_vel[3*i+2] = vel_aux.real( );

        // Print out results
        std::cout << std::endl;
        std::cout << "-> FIELD POINT[" << i << "]: "; print_vector( 3, &(field_points[3*i]), 0, 3 );
        std::cout << " --> Hess:    "; print_vector( 3, &(hess_vel[3*i]), 0, 6 );
        std::cout << " --> Newman:  "; print_vector( 3, &(nw_vel[3*i]), 0, 6 );
        std::cout << " --> Quad:    "; print_vector( 3, &(quad_vel[3*i]), 0, 6 );
        std::cout << " --> Is inside: " << panel->is_inside( &(field_points[3*i]) ) << std::endl;
        std::cout << std::endl;
    }
}


void test_2( void )
{
    // Create panel for checking 
    PanelGeom* panel = new PanelGeom;

    panel->x[0] = -1.0;
    panel->x[1] =  1.0;
    panel->x[2] =  1.0;
    panel->x[3] = -1.0;

    panel->y[0] =  -1.0;
    panel->y[1] =  -1.0;
    panel->y[2] =  1.0;
    panel->y[3] =  1.0;

    panel->z[0] =  0.0;
    panel->z[1] =  0.0;
    panel->z[2] =  0.0;
    panel->z[3] =  0.0;
    
    panel->num_nodes = 4;

    cusfloat cog[3] = { 0.0, 0.0, 0.0 };
    panel->calculate_properties( cog );

    // Define field points
    const int fp_np = 7;
    cusfloat field_points[3*fp_np];

    field_points[0] = 0.0;
    field_points[1] = 0.0;
    field_points[2] = 0.0;

    field_points[3] = 1e-6;
    field_points[4] = 1e-6;
    field_points[5] = 0.0;

    field_points[6] = 1e-5;
    field_points[7] = 1e-5;
    field_points[8] = 0.0;

    field_points[9]  = 1e-4;
    field_points[10] = 1e-4;
    field_points[11] = 0.0;

    field_points[12] = 1e-3;
    field_points[13] = 1e-3;
    field_points[14] = 0.0;

    field_points[15] = 1e-2;
    field_points[16] = 1e-2;
    field_points[17] = 0.0;

    field_points[18] = 1e-1;
    field_points[19] = 1e-1;
    field_points[20] = 0.0;
    
    // Allocate memory for the velocity vectors computed
    // used the three different formulations
    const int   vel_np                  = 3 * fp_np;
    cusfloat    hess_vel[3*fp_np];      clear_vector( vel_np, hess_vel );
    cusfloat    nw_vel[3*fp_np];        clear_vector( vel_np, nw_vel );
    cusfloat    quad_vel[3*fp_np];      clear_vector( vel_np, quad_vel );
    cuscomplex  vel_aux( 0.0, 0.0 );

    std::cout << "***************************************" << std::endl;
    std::cout << "***************  TEST_1 ***************" << std::endl;
    std::cout << "***************************************" << std::endl;
    for ( int i=0; i<fp_np; i++ )
    {
        // Calculate velocity vector using Hess & Smith formulation
        calculate_source_velocity_hess(
                                            panel,
                                            &(field_points[3*i]), 
                                            0,
                                            0, 
                                            &(hess_vel[3*i])
                                        );

        // Calculate velocity vector using Newman formulation
        calculate_source_velocity_newman(
                                            panel,
                                            &(field_points[3*i]), 
                                            0,
                                            0, 
                                            &(nw_vel[3*i])
                                        );

        // Define lambda functions for the derivatives of the 
        // Rankine sources
        auto rank_dx_fcn    =   [field_points, i]
                                (cusfloat xi, cusfloat eta, cusfloat x, cusfloat y, cusfloat z)
                                {
                                    cusfloat dx = field_points[3*i] - x;
                                    cusfloat dy = field_points[3*i+1] - y;
                                    cusfloat dz = field_points[3*i+2] - z;
                                    cusfloat R  = std::sqrt( pow2s( dx ) + pow2s( dy ) + pow2s( dz ) );

                                    return dx / pow3s( R );
                                };

        auto rank_dy_fcn    =   [field_points, i]
                                (cusfloat xi, cusfloat eta, cusfloat x, cusfloat y, cusfloat z)
                                {
                                    cusfloat dx = field_points[3*i] - x;
                                    cusfloat dy = field_points[3*i+1] - y;
                                    cusfloat dz = field_points[3*i+2] - z;
                                    cusfloat R  = std::sqrt( pow2s( dx ) + pow2s( dy ) + pow2s( dz ) );

                                    return dy / pow3s( R );
                                };
        
        auto rank_dz_fcn    =   [field_points, i]
                                (cusfloat xi, cusfloat eta, cusfloat x, cusfloat y, cusfloat z)
                                {
                                    cusfloat dx = field_points[3*i] - x;
                                    cusfloat dy = field_points[3*i+1] - y;
                                    cusfloat dz = field_points[3*i+2] - z;
                                    cusfloat R  = std::sqrt( pow2s( dx ) + pow2s( dy ) + pow2s( dz ) );

                                    return dz / pow3s( R );
                                };

        // vel_aux =   adaptive_quadrature_panel(
        //                                         panel,
        //                                         rank_dx_fcn,
        //                                         ABS_EPS,
        //                                         REL_EPS,
        //                                         false,
        //                                         true,
        //                                         1
        //                                     );
        // quad_vel[3*i]   = vel_aux.real( );
        
        // vel_aux =   adaptive_quadrature_panel(
        //                                         panel,
        //                                         rank_dy_fcn,
        //                                         ABS_EPS,
        //                                         REL_EPS,
        //                                         false,
        //                                         true,
        //                                         1
        //                                     );
        // quad_vel[3*i+1] = vel_aux.real( );
        
        // vel_aux =   adaptive_quadrature_panel(
        //                                         panel,
        //                                         rank_dz_fcn,
        //                                         ABS_EPS,
        //                                         REL_EPS,
        //                                         false,
        //                                         true,
        //                                         1
        //                                     );
        // quad_vel[3*i+2] = vel_aux.real( );

        // Print out results
        std::cout << std::endl;
        std::cout << "-> FIELD POINT[" << i << "]: "; print_vector( 3, &(field_points[3*i]), 0, 3 );
        std::cout << " --> Hess:    "; print_vector( 3, &(hess_vel[3*i]), 0, 6 );
        std::cout << " --> Newman:  "; print_vector( 3, &(nw_vel[3*i]), 0, 6 );
        std::cout << " --> Quad:    "; print_vector( 3, &(quad_vel[3*i]), 0, 6 );
        std::cout << " --> Is inside: " << panel->is_inside( &(field_points[3*i]) ) << std::endl;
        std::cout << std::endl;
    }
}


int main( void )
{   
    // test_0( );
    // test_1( );
    test_2( );

    return 0;
}