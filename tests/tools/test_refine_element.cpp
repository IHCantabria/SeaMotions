
// Include general usage libraries
#include <fstream>
#include <string>

// Include local modules
#include "../../src/containers/panel_geom.hpp"
#include "../../src/math/math_tools.hpp"
#include "../../src/mesh/tools.hpp"

/***********************************************/
/*********** Declare Module objects ************/
/***********************************************/
void    ref_square( 
                                    void 
                    );

void    ref_triangle( 
                                    void 
                    );

void    write_element_i(
                                    std::ofstream&  outfile,
                                    PanelGeom*      panel
                       );   

void    write_element_refinement(   
                                    std::string     file_path,
                                    PanelGeom*      panel,
                                    PanelGeomList*  panel_list
                               );


/***********************************************/
/************ Define Module objects ************/
/***********************************************/
int main( void )
{
    // Launch square refinement test
    ref_square( );

    // Launch triangle refinement test
    ref_triangle( );

    return 0;
}


void ref_square( void )
{
    // Define parent geometry
    PanelGeom panel;

    panel.num_nodes = 4;

    panel.x[0] = -1.0;
    panel.x[1] =  1.0;
    panel.x[2] =  1.0;
    panel.x[3] = -1.0;

    panel.y[0] = -1.0;
    panel.y[1] = -1.0;
    panel.y[2] =  1.0;
    panel.y[3] =  1.0;

    panel.z[0] =  0.0;
    panel.z[1] =  0.0;
    panel.z[2] =  0.0;
    panel.z[3] =  0.0;

    // Calculate local properties
    panel.calculate_properties(  );

    // Refine element
    PanelGeomList* panel_list;
    refine_element( &panel, panel_list );

    // Calculate sumation area for the panel refinement
    cusfloat cum_area = 0.0;
    for ( int i=0; i<panel_list->panel_np; i++ )
    {
        cum_area += panel_list->panels[i]->area;
    }

    // Check if the summation area
    if ( !assert_scalar_equality( panel.area, cum_area, 1e-9 ) )
    {
        std::cerr << "Test test_refine_element/ref_square failed!" << std::endl;
        std::cerr << " - Area: " << panel.area << " - Refinement Area: " << cum_area << std::endl;
        throw std::runtime_error( "" );
    }

}


void ref_triangle( void )
{
    // Define parent geometry
    PanelGeom panel;

    panel.num_nodes = 3;

    panel.x[0] =  0.0;
    panel.x[1] =  1.0;
    panel.x[2] = -1.0;

    panel.y[0] =  0.0;
    panel.y[1] =  0.0;
    panel.y[2] =  1.0;

    panel.z[0] =  0.0;
    panel.z[1] =  0.0;
    panel.z[2] =  0.0;

    // Calculate local properties
    panel.calculate_properties(  );

    // Refine element
    PanelGeomList* panel_list;
    refine_element( &panel, panel_list );

    // Calculate sumation area for the panel refinement
    cusfloat cum_area = 0.0;
    for ( int i=0; i<panel_list->panel_np; i++ )
    {
        cum_area += panel_list->panels[i]->area;
    }

    // Check if the summation area
    if ( !assert_scalar_equality( panel.area, cum_area, 1e-9 ) )
    {
        std::cerr << "Test test_refine_element/ref_square failed!" << std::endl;
        std::cerr << " - Area: " << panel.area << " - Refinement Area: " << cum_area << std::endl;
        throw std::runtime_error( "" );
    }

}


void write_element_i(
                        std::ofstream&  outfile,
                        PanelGeom*      panel
                    )
{
    outfile << "Num.Nodes: " << panel->num_nodes << "\n";

    for ( int i=0; i<panel->num_nodes; i++ )
    {
        outfile << panel->x[i] << " " << panel->y[i] << " " << panel->z[i] << "\n";
    }
}

void write_element_refinement( 
                                std::string     file_path,
                                PanelGeom*      panel,
                                PanelGeomList*  panel_list
                            )
{
    // Open file unit
    std::ofstream outfile( file_path );

    // Define Parent element geometry
    outfile << "Parent_Element\n";
    write_element_i( outfile, panel );

    // Define child elements geometry
    outfile << "Child_Elements\n";
    outfile << "Num.Elements: " << panel_list->panel_np << "\n";
    for ( int i=0; i<panel_list->panel_np; i++ )
    {
        outfile << "Element_" << i << "\n";
        write_element_i( outfile, panel_list->panels[i] );
    }

    // Close file unit
    outfile.close( );
}