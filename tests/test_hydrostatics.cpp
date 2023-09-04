
// Include general usage libraries
#include <fstream>
#include <iostream>

// Include local modules
#include "../src/containers/panel_geom.hpp"
#include "../src/math/integration.hpp"
#include "../src/mesh/mesh.hpp"
#include "../src/hydrostatics.hpp"
#include "../src/tools.hpp"


// Define precision level
cusfloat ABS_EPS = 1E-2;
cusfloat REL_EPS = 1E-2;


/******************************************************/
/************** Define Auxiliary Macros ***************/
/******************************************************/
#define CHECK_STAB_PARAMETER( condition, signal_name )                                                  \
    if ( condition )                                                                                    \
    {                                                                                                   \
        std::cerr << std::endl;                                                                         \
        std::cerr << "ERROR - Signal: " << signal_name;                                                 \
        std::cerr << " does not match the reference values with the required precision" << std::endl;   \
        std::cerr << "TEST test_hydrostatics/" << test_name << " failed!" << std::endl;                 \
        throw std::exception( "" );                                                                     \
    }                                                                                                   \


/******************************************************/
/*************** Define Reference Data ****************/
/******************************************************/
struct RefData
{
private:
    // Define class methods
    void _load_data( std::string fipath )
    {
        // Define auxiliar variablest to read lines
        std::string aux_str;

        // Open file unit
        std::ifstream infile( fipath );

        // Read gravitational acceleration
        infile >> aux_str >> this->grav_acc;

        // Read water density
        infile >> aux_str >> this->rho_water;

        // Read structural mass
        infile >> aux_str >> this->mass;

        // Read cetre of gravity
        infile >> aux_str >> this->cog[0] >> this->cog[1] >> this->cog[2];

        // Read radius of inertia
        infile >> aux_str >> this->rad_inertia[0] >> this->rad_inertia[1] >> this->rad_inertia[2];

        // Read volume
        infile >> aux_str >> this->volume;

        // Read centre of buoyancy
        infile >> aux_str >> this->cob[0] >> this->cob[1] >> this->cob[2];

        // Read water line area
        infile >> aux_str >> this->wl_area;

        // Read water line area centre of gravity
        infile >> aux_str >> this->wl_area_cog[0] >> this->wl_area_cog[1] >> this->wl_area_cog[2];

        // Read water line area interia around X axis
        infile >> aux_str >> this->wl_area_ixx;

        // Read water line area interia around Y axis
        infile >> aux_str >> this->wl_area_iyy;

        // Read KB
        infile >> aux_str >> this->kb;

        // Read metacentric radius around X axis
        infile >> aux_str >> this->bmx;

        // Read metacentric radius around Y axis
        infile >> aux_str >> this->bmy;

        // Read metacentric height around X axis
        infile >> aux_str >> this->gmx;

        // Read metacentric height around Y axis
        infile >> aux_str >> this->gmy;

        // Close file unit
        infile.close( );
    }

public:
    // Define class attributes
    cusfloat bmx            = 0.0;
    cusfloat bmy            = 0.0;
    cusfloat cob[3]         = { 0.0, 0.0, 0.0 };
    cusfloat cog[3]         = { 0.0, 0.0, 0.0 };
    cusfloat gmx            = 0.0;
    cusfloat gmy            = 0.0;
    cusfloat grav_acc       = 0.0;
    cusfloat kb             = 0.0;
    cusfloat mass           = 0.0;
    cusfloat rad_inertia[3] = { 0.0, 0.0, 0.0 };
    cusfloat rho_water      = 0.0;
    cusfloat volume         = 0.0;
    cusfloat wl_area        = 0.0;
    cusfloat wl_area_cog[3] = { 0.0, 0.0, 0.0 };
    cusfloat wl_area_ixx    = 0.0;
    cusfloat wl_area_iyy    = 0.0;

    // Define class constructors and destructor
    RefData( std::string fipath )
    {
        this->_load_data( fipath );
    }

};



/******************************************************/
/*************** Define Test Functions ****************/
/******************************************************/
void launch_test( 
                    std::string test_name,
                    std::string props_fipath,
                    std::string mesh_fipath
                )
{
    // Read hydrostatic properties
    RefData ref_data( props_fipath );

    // Read mesh
    Mesh mesh( mesh_fipath );

    // Calculate hydrostatics
    Hydrostatics hydro( 
                            &mesh,
                            ref_data.rho_water,
                            ref_data.grav_acc,
                            ref_data.mass,
                            ref_data.cog,
                            ref_data.rad_inertia
                        );

    // Check displacement
    CHECK_STAB_PARAMETER(
                            !assert_scalar_equality( 
                                                        ref_data.mass, 
                                                        hydro.displacement, 
                                                        ABS_EPS,
                                                        REL_EPS
                                                    ),
                            "DISPLACEMENT"
                        );

    // Check COB
    CHECK_STAB_PARAMETER(
                            !assert_vector_equality( 
                                                        3,
                                                        ref_data.cob, 
                                                        hydro.cob, 
                                                        ABS_EPS,
                                                        REL_EPS
                                                    ),
                            "COB"
                        );

    // Check water line area
    CHECK_STAB_PARAMETER(
                            !assert_scalar_equality( 
                                                        ref_data.wl_area, 
                                                        hydro.wl_area, 
                                                        ABS_EPS,
                                                        REL_EPS
                                                    ),
                            "WATERLINE AREA"
                        );

    // Check water line area centre of gravity
    CHECK_STAB_PARAMETER(
                            !assert_vector_equality( 
                                                        3,
                                                        ref_data.wl_area_cog, 
                                                        hydro.wl_area_cog, 
                                                        ABS_EPS,
                                                        REL_EPS
                                                    ),
                            "WATERLINE AREA COG"
                        );

    // Check water line area inertia around X axis
    CHECK_STAB_PARAMETER(
                            !assert_scalar_equality( 
                                                        ref_data.wl_area_ixx, 
                                                        hydro.wl_area_ixx, 
                                                        ABS_EPS,
                                                        REL_EPS
                                                    ),
                            "WATERLINE AREA IXX"
                        );
    
    // Check water line area inertia around Y axis
    CHECK_STAB_PARAMETER(
                            !assert_scalar_equality( 
                                                        ref_data.wl_area_iyy, 
                                                        hydro.wl_area_iyy, 
                                                        ABS_EPS,
                                                        REL_EPS
                                                    ),
                            "WATERLINE AREA IYY"
                        );

    // Check KB
    CHECK_STAB_PARAMETER(
                            !assert_scalar_equality( 
                                                        ref_data.kb, 
                                                        hydro.kb, 
                                                        ABS_EPS,
                                                        REL_EPS
                                                    ),
                            "KB"
                        );

    // Check metracentric radius around X axis
    CHECK_STAB_PARAMETER(
                            !assert_scalar_equality( 
                                                        ref_data.bmx, 
                                                        hydro.bmx, 
                                                        ABS_EPS,
                                                        REL_EPS
                                                    ),
                            "BMX"
                        );
    
    // Check metracentric radius around Y axis
    CHECK_STAB_PARAMETER(
                            !assert_scalar_equality( 
                                                        ref_data.bmy, 
                                                        hydro.bmy, 
                                                        ABS_EPS,
                                                        REL_EPS
                                                    ),
                            "BMY"
                        );

    // Check metracentric height around X axis
    CHECK_STAB_PARAMETER(
                            !assert_scalar_equality( 
                                                        ref_data.gmx, 
                                                        hydro.gmx, 
                                                        ABS_EPS,
                                                        REL_EPS
                                                    ),
                            "GMX"
                        );

    // Check metracentric height around Y axis
    CHECK_STAB_PARAMETER(
                            !assert_scalar_equality( 
                                                        ref_data.gmy, 
                                                        hydro.gmy, 
                                                        ABS_EPS,
                                                        REL_EPS
                                                    ),
                            "GMY"
                        );

}


int main( int argc, char* argv[ ] )
{
    // Read command line arguments
    if ( !check_num_cmd_args( argc, 4 ) )
    {
        return 1;
    }

    std::string dike_props_fipath( argv[1] );
    std::string dike_mesh_fipath( argv[2] );
    std::string semisub_props_fipath( argv[3] );
    std::string semisub_mesh_fipath( argv[4] );

    // Launch test for Dike floating object
    launch_test( 
                    "DIKE",
                    dike_props_fipath,
                    dike_mesh_fipath
                );

    // Launch test for Dike floating object
    launch_test( 
                    "SEMI-SUBMERSIBLE",
                    semisub_props_fipath,
                    semisub_mesh_fipath
                );

    return 0;
}