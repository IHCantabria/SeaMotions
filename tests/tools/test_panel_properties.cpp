
/*
 * Copyright (c) 2025 Sergio Fernández Ruano / IHCantabria
 *
 * This file is part of SeaMotions Software.
 *
 * SeaMotions is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SeaMotions is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */


// Include local modules
#include "../../src/config.hpp"
#include "../../src/math/math_tools.hpp"


struct PanelProps
{
    cusfloat area;
    cusfloat centroid[3];
    cusfloat moments_fo[3];
    cusfloat moments_so[3];
    int      num_nodes;
    cusfloat x[4];
    cusfloat y[4];
    cusfloat z[4];
};


void define_quadrilateral_0( PanelProps& panel )
{
    // Define node positions
    panel.num_nodes = 4;
    panel.x[0] = -3.6610856;    panel.y[0] = -2.4552272;    panel.z[0] = 0.0;
    panel.x[1] = 4.0588555;     panel.y[1] = -3.2909639;    panel.z[1] = 0.0;
    panel.x[2] = 4.6821168;     panel.y[2] = 3.1399595;     panel.z[2] = 0.0;
    panel.x[3] = -5.8991602;    panel.y[3] = 4.7831029;     panel.z[3] = 0.0;

    // Define expected geometric properties
    panel.area          = 61.5402647;
    panel.centroid[0]   = -0.273091179;
    panel.centroid[1]   = 0.725623761;
    panel.moments_fo[0] = -16.8061035;
    panel.moments_fo[1] = 44.6550783;
    panel.moments_so[0] = 278.648016; 
    panel.moments_so[1] = 447.558618;
    panel.moments_so[2] = -100.858527;
    
}


void define_quadrilateral_1( PanelProps& panel )
{
    // Define node positions
    panel.num_nodes = 4;
    panel.x[0] = -5.5770783;    panel.y[0] = -1.6548556;    panel.z[0] = 0.0;
    panel.x[1] = 2.8936093;     panel.y[1] = -3.6662898;    panel.z[1] = 0.0;
    panel.x[2] = 1.7745720;     panel.y[2] = 0.5265589;     panel.z[2] = 0.0;
    panel.x[3] = -2.3049565;    panel.y[3] = 2.9204489;     panel.z[3] = 0.0;

    // Define expected geometric properties
    panel.area          = 29.8818132;
    panel.centroid[0]   = -1.07128445;
    panel.centroid[1]   = -0.624712501;
    panel.moments_fo[0] = -32.0119218;
    panel.moments_fo[1] = -18.6675423;
    panel.moments_so[0] = 70.9665805; 
    panel.moments_so[1] = 145.09503;
    panel.moments_so[2] = -2.87223728;

}


void define_quadrilateral_2( PanelProps& panel )
{
    // Define node positions
    panel.num_nodes = 4;
    panel.x[0] = -7.3160090;    panel.y[0] = -0.1111885;    panel.z[0] = 0.0;
    panel.x[1] = 0.5597474;     panel.y[1] = -2.4342504;    panel.z[1] = 0.0;
    panel.x[2] = 6.6648751;     panel.y[2] = 0.1012899;     panel.z[2] = 0.0;
    panel.x[3] = -0.6584451;    panel.y[3] = 2.8634706;     panel.z[3] = 0.0;

    // Define expected geometric properties
    panel.area          = 37.1628297;
    panel.centroid[0]   = -0.249943854;
    panel.centroid[1]   = 0.13977487;
    panel.moments_fo[0] = -9.28862089;
    panel.moments_fo[1] = 5.19442969;
    panel.moments_so[0] = 44.3584639; 
    panel.moments_so[1] = 307.378911;
    panel.moments_so[2] = -6.42268923;

}


void define_triangle_0( PanelProps& panel )
{
    // Define node positions
    panel.num_nodes = 3;
    panel.x[0] = -0.6666667;    panel.y[0] = -0.6666667;    panel.z[0] = 0.0;
    panel.x[1] = 1.333333;      panel.y[1] = -0.6666667;    panel.z[1] = 0.0;
    panel.x[2] = -0.6666667;    panel.y[2] = 1.333333;      panel.z[2] = 0.0;

    // Define expected geometric properties
    panel.area          = 2.0;
    panel.centroid[0]   = 0.0;
    panel.centroid[1]   = 0.0;
    panel.moments_fo[0] = 0.0;
    panel.moments_fo[1] = 0.0;
    panel.moments_so[0] = 0.444444444; 
    panel.moments_so[1] = 0.444444444;
    panel.moments_so[2] = -0.222222222;

}


void define_triangle_1( PanelProps& panel )
{
    // Define node positions
    panel.num_nodes = 3;
    panel.x[0] = -0.3187483;    panel.y[0] = -0.82024656;   panel.z[0] = 0.0;
    panel.x[1] = 3.1109834;     panel.y[1] = -0.5014973;    panel.z[1] = 0.0;
    panel.x[2] = -2.7922351;    panel.y[2] = 1.3217429;     panel.z[2] = 0.0;

    // Define expected geometric properties
    panel.area          = 4.06743283;
    panel.centroid[0]   = 0.0;
    panel.centroid[1]   = 0.0;
    panel.moments_fo[0] = 0.0;
    panel.moments_fo[1] = 0.0;
    panel.moments_so[0] = 0.905446754; 
    panel.moments_so[1] = 5.95756708;
    panel.moments_so[2] = -1.69114194;

}


void define_triangle_2( PanelProps& panel )
{
    // Define node positions
    panel.num_nodes = 3;
    panel.x[0] = -1.1594326;    panel.y[0] = -2.0753431;    panel.z[0] = 0.0;
    panel.x[1] = 2.2576768;     panel.y[1] = -2.1179682;    panel.z[1] = 0.0;
    panel.x[2] = 1.0484606;     panel.y[2] = 6.6170103;     panel.z[2] = 0.0;

    // Define expected geometric properties
    panel.area          = 14.898417;
    panel.centroid[0]   = 0.715568283;
    panel.centroid[1]   = 0.807899675;
    panel.moments_fo[0] = 10.6608346;
    panel.moments_fo[1] = 12.0364262;
    panel.moments_so[0] = 72.5701644; 
    panel.moments_so[1] = 15.0834054;
    panel.moments_so[2] = 12.1238025;

}


template<typename T, typename U>
bool test_panel_properties(
                                T* define_func_panel,
                                U* calculate_panel_geometric_properties,
                                const std::optional<cusfloat*>& ref_sys = std::nullopt
                            )
{
    // Create panel properties structure
    PanelProps panel_props;

    // Define panel using the input function
    define_func_panel( panel_props );

    // Calculate panel geometric properties
    cusfloat area           = 0.0;
    cusfloat centroid[3]    = { 0.0, 0.0, 0.0 };
    cusfloat moments_fo[3]  = { 0.0, 0.0, 0.0 };
    cusfloat moments_so[3]  = { 0.0, 0.0, 0.0 };

    calculate_panel_geometric_properties(
                                                panel_props.x,
                                                panel_props.y,
                                                panel_props.z,
                                                area,
                                                centroid,
                                                moments_fo,
                                                moments_so,
                                                ref_sys
                                            );

    // Check results
    const cusfloat abs_tol = 1e-4;
    
    bool test_area  = assert_scalar_equality( area, panel_props.area, abs_tol );
    bool test_cog   = assert_vector_equality( 3, centroid, panel_props.centroid, abs_tol );
    bool test_mfo   = assert_vector_equality( 3, moments_fo, panel_props.moments_fo, abs_tol );
    bool test_mso   = assert_vector_equality( 3, moments_so, panel_props.moments_so, abs_tol );
    bool test_pass  = test_area && test_cog && test_mfo && test_mso;

    return test_pass;
    
}


int main( void )
{
    // Declare test variables
    bool test_passed = false;

    // Launch Triangle 0 test
    test_passed = test_panel_properties( define_triangle_0, triangle_geom_properties );
    if ( !test_passed )
    {
        std::cerr << "test_panel_properties: Triangle 0 test failed." << std::endl;
        return 1;
    }

    // Launch Triangle 1 test
    test_passed = test_panel_properties( define_triangle_1, triangle_geom_properties );
    if ( !test_passed )
    {
        std::cerr << "test_panel_properties: Triangle 1 test failed." << std::endl;
        return 1;
    }

    // Launch Triangle 2 test
    cusfloat ref_sys[3] = { 0.0, 0.0, 0.0 };
    test_passed = test_panel_properties( define_triangle_2, triangle_geom_properties, ref_sys );
    if ( !test_passed )
    {
        std::cerr << "test_panel_properties: Triangle 2 test failed." << std::endl;
        return 1;
    }

    // Launch Quadrilateral 0 test
    test_passed = test_panel_properties( define_quadrilateral_0, quad_geom_properties );
    if ( !test_passed )
    {
        std::cerr << "test_panel_properties: Quadrilateral 0 test failed." << std::endl;
        return 1;
    }

    // Launch Quadrilateral 1 test
    test_passed = test_panel_properties( define_quadrilateral_1, quad_geom_properties );
    if ( !test_passed )
    {
        std::cerr << "test_panel_properties: Quadrilateral 1 test failed." << std::endl;
        return 1;
    }

    // Launch Quadrilateral 2 test
    test_passed = test_panel_properties( define_quadrilateral_2, quad_geom_properties );
    if ( !test_passed )
    {
        std::cerr << "test_panel_properties: Quadrilateral 2 test failed." << std::endl;
        return 1;
    }

    
    return 0;

}