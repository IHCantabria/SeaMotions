
// Include general usage libraries
#include <filesystem>
#include <iostream>
#include <sstream>

// Include local modules
#include "../containers/body_def.hpp"
#include "input.hpp"
#include "reader.hpp"
#include "../tools.hpp"
#include "../version.hpp"

// Manage namespaces
namespace fs = std::filesystem;


void read_bodies( 
                    std::string folder_path,
                    Input*      input
                )
{
    // Allocate space for the bodies description
    input->bodies = new BodyDef*[input->bodies_np];

    // Loop over body definitions to load them
    for ( int i=0; i<input->bodies_np; i++ )
    {
        // Create new body instance
        input->bodies[i] = new BodyDef;

        // Load body data
        read_body(
                    input,
                    folder_path,
                    input->bodies_finame[i],
                    input->bodies[i]
                );

    }

    // Set bodies as read
    input->is_bodies = true;
}


void read_body(
                    Input*      input,
                    std::string folder_path,
                    std::string target_file,
                    BodyDef*    body
                )
{
    // Generate auxiliar variables
    int         line_count      = 0;
    std::string read_signal     = "";
    std::string target_signal   = "";
    std::string _version        = "";

    // Generate case.input.dat file path
    fs::path folder_path_( folder_path );
    fs::path file_name( target_file );
    fs::path file_path = folder_path_ / file_name;

    // Open file unit
    std::ifstream infile;
    infile.open( file_path, std::ios_base::in );
    CHECK_FILE_UNIT_STATUS( infile, file_path );

    // Read file header line
    _skip_header( infile, line_count, 1 );

    // Read file version
    target_signal   = "version";
    read_signal     = _read_channel_value( infile, _version );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );
    CHECK_INPUT_FILE_VERSION( VERSION_LABEL, _version, file_path );

    //////////////////////////////////////////////
    /************* Mass Properties **************/
    //////////////////////////////////////////////

    // Skip header
    _skip_header( infile, line_count, 3 );

    // Read mass
    target_signal   = "BodyMass";
    read_signal     = _read_channel_value( infile, body->mass );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read COGX
    target_signal   = "COGX";
    read_signal     = _read_channel_value( infile, body->cog[0] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read COGY
    target_signal   = "COGY";
    read_signal     = _read_channel_value( infile, body->cog[1] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read COGZ
    target_signal   = "COGZ";
    read_signal     = _read_channel_value( infile, body->cog[2] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read switch to define how the intertia is defined, by the 
    // radius of gyrantion or by the inertial matrix components
    target_signal   = "IBR";
    read_signal     = _read_channel_value( infile, body->interia_by_rad );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read IXX
    target_signal   = "IXX";
    read_signal     = _read_channel_value( infile, body->inertia[0] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read IXY
    target_signal   = "IXY";
    read_signal     = _read_channel_value( infile, body->inertia[1] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read IXZ
    target_signal   = "IXZ";
    read_signal     = _read_channel_value( infile, body->inertia[2] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read IYY
    target_signal   = "IYY";
    read_signal     = _read_channel_value( infile, body->inertia[3] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read IYZ
    target_signal   = "IYZ";
    read_signal     = _read_channel_value( infile, body->inertia[4] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read IZZ
    target_signal   = "IZZ";
    read_signal     = _read_channel_value( infile, body->inertia[5] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read RXX
    target_signal   = "RXX";
    read_signal     = _read_channel_value( infile, body->rad_inertia[0] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read RYY
    target_signal   = "RYY";
    read_signal     = _read_channel_value( infile, body->rad_inertia[1] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read RZZ
    target_signal   = "RZZ";
    read_signal     = _read_channel_value( infile, body->rad_inertia[2] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    if ( body->interia_by_rad )
    {
        body->inertia[0] = body->mass * pow2s( body->rad_inertia[0] );
        body->inertia[3] = body->mass * pow2s( body->rad_inertia[1] );
        body->inertia[5] = body->mass * pow2s( body->rad_inertia[2] );
    }

    //////////////////////////////////////////////
    /************* Mesh Description *************/
    //////////////////////////////////////////////

    // Skip header
    _skip_header( infile, line_count, 3 );

    // Read mesh body name
    target_signal   = "BodyName";
    read_signal     = _read_channel_value( infile, body->mesh_body_name );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read mesh file name
    target_signal   = "MeshFile";
    read_signal     = _read_channel_value( infile, body->mesh_finame );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read usage of internal lid
    target_signal   = "LidType";
    read_signal     = _read_channel_value( infile, body->lid_type );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Load mesh for the current body
    fs::path mesh_foname_( std::string( "mesh" ) );
    fs::path mesh_finame_( body->mesh_finame );
    fs::path mesh_fipath_ = folder_path_ / mesh_foname_ / mesh_finame_;

    body->mesh      = new   Mesh( 
                                    mesh_fipath_.string( ),
                                    body->mesh_body_name,
                                    body->cog,
                                    DIFFRAC_PANEL_CODE                                
                                  );
    body->is_mesh   = true;

    // Load lid if it is required by the user
    if ( body->lid_type == 1 )
    {  
        // Define lid name
        std::stringstream ss;
        ss << body->mesh_body_name << "_lid";
        
        // Load lid mesh
        Mesh* lid_mesh = new Mesh(
                                    mesh_fipath_.string( ),
                                    ss.str( ),
                                    body->cog,
                                    LID_PANEL_CODE
                                );

        // Joint body and lid meshes
        std::vector<Mesh*> meshes;
        meshes.push_back( body->mesh );
        meshes.push_back( lid_mesh );

        body->mesh  = new Mesh(
                                    meshes,
                                    body->cog
                                );

        // Delete body and lid meshes
        for ( auto mi: meshes )
        {
            delete mi;
        }
    }

    // Load QTF free surface
    if ( input->out_qtf_so_model > 0 )
    {
        // Define QTF free surface name
        std::stringstream ss;
        ss << body->mesh_body_name << "_fs_qtf";

        // Load QTF free surface mesh
        body->mesh_fs_qtf       = new Mesh(
                                            mesh_fipath_.string( ),
                                            ss.str( ),
                                            body->cog,
                                            LID_PANEL_CODE
                                        );
        body->is_mesh_fs_qtf    = true;
    }

    // Close file unit
    infile.close();

}


void    read_case( 
                    std::string folder_path,
                    Input*      input
                )
{
    // Generate auxiliar variables
    int         line_count      = 0;
    std::string read_signal     = "";
    std::string target_file     = "case.input.dat";
    std::string target_signal   = "";
    std::string _version        = "";

    // Generate case.input.dat file path
    fs::path folder_path_( folder_path );
    fs::path file_name( "case.input.dat" );
    fs::path file_path = folder_path_ / file_name;

    // Open file unit
    std::ifstream infile;
    infile.open( file_path, std::ios_base::in );
    CHECK_FILE_UNIT_STATUS( infile, file_path );

    // Read file header line
    _skip_header( infile, line_count, 1 );

    // Read file version
    target_signal   = "version";
    read_signal     = _read_channel_value( infile, _version );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );
    CHECK_INPUT_FILE_VERSION( VERSION_LABEL, _version, file_path );

    //////////////////////////////////////////////
    /************** Solver Controls *************/
    //////////////////////////////////////////////

    // Skip header
    _skip_header( infile, line_count, 3 );

    // Read flag to block the quadrature adaptation algorithm
    target_signal   = "BlockAdapt";
    read_signal     = _read_channel_value( infile, input->is_block_adaption );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read flag to use the fast solver configuration
    target_signal   = "FastSolver";
    read_signal     = _read_channel_value( infile, input->is_fast_solver );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read gauss order for the numerical integration
    target_signal   = "GaussOrder";
    read_signal     = _read_channel_value( infile, input->gauss_order );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read absolute error for green function normal derivative integration over panel
    target_signal   = "GFDnAbsErr";
    read_signal     = _read_channel_value( infile, input->gfdn_abs_err );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read relative error for green function normal derivative integration over panel
    target_signal   = "GFDnRelErr";
    read_signal     = _read_channel_value( infile, input->gfdn_rel_err );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read number of Kotchin expansion coefficients
    target_signal   = "KochinNC";
    read_signal     = _read_channel_value( infile, input->kochin_np );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read method to solver the logarithmic singularity
    target_signal   = "LogSingAna";
    read_signal     = _read_channel_value( infile, input->is_log_sin_ana );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read polynomial order to interpolate the solution
    target_signal   = "PolyOrder";
    read_signal     = _read_channel_value( infile, input->poly_order );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read absolute error for potential integration over panel
    target_signal   = "PotAbsErr";
    read_signal     = _read_channel_value( infile, input->pot_abs_err );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read relative error for potential integration over panel
    target_signal   = "PotRelErr";
    read_signal     = _read_channel_value( infile, input->pot_rel_err );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read absolute error for pressure integration over panel
    target_signal   = "PressAbsErr";
    read_signal     = _read_channel_value( infile, input->press_abs_err );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read relative error for pressure integration over panel
    target_signal   = "PressRelErr";
    read_signal     = _read_channel_value( infile, input->press_rel_err );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read QTF second order potential model
    target_signal   = "QTFSOModel";
    read_signal     = _read_channel_value( infile, input->out_qtf_so_model );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read water line points detection precision
    target_signal   = "WLDetPrec";
    read_signal     = _read_channel_value( infile, input->wl_det_prec );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    //////////////////////////////////////////////
    /************** Body Definition *************/
    //////////////////////////////////////////////

    // Skip header
    _skip_header( infile, line_count, 3 );

    // Read bodies name definition
    target_signal = "BodyFN";
    _read_channel_list(
                            infile,
                            target_file,
                            target_signal,
                            line_count,
                            input->bodies_finame
                        );
    input->bodies_np = input->bodies_finame.size( );

    //////////////////////////////////////////////
    /************** Output Channels *************/
    //////////////////////////////////////////////

    // Skip header
    _skip_header( infile, line_count, 3 );

    // Read output flag for diffraction force
    target_signal   = "OutDiffrac";
    read_signal     = _read_channel_value( infile, input->out_diffrac );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for Froude-Krylov force
    target_signal   = "OutFK";
    read_signal     = _read_channel_value( infile, input->out_fk );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for hydromechanic coefficients
    target_signal   = "OutHydMech";
    read_signal     = _read_channel_value( infile, input->out_hydmech );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for hydrostatic stiffness matrix
    target_signal   = "OutHydStiff";
    read_signal     = _read_channel_value( infile, input->out_hydstiff );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for pressure over the panels
    target_signal   = "OutPress";
    read_signal     = _read_channel_value( infile, input->out_pressure );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for mean drift
    target_signal   = "OutMDrift";
    read_signal     = _read_channel_value( infile, input->out_mdrift );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for quadratic transfer functions
    target_signal   = "OutQTF";
    read_signal     = _read_channel_value( infile, input->out_qtf );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for quadratic transfer functions components
    target_signal   = "OutQTFComp";
    read_signal     = _read_channel_value( infile, input->out_qtf_comp );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for response amplitude operator
    target_signal   = "OutRAOs";
    read_signal     = _read_channel_value( infile, input->out_raos );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for sources intensity over the panels
    target_signal   = "OutSources";
    read_signal     = _read_channel_value( infile, input->out_sources );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for structural mass
    target_signal   = "OutStMass";
    read_signal     = _read_channel_value( infile, input->out_struct_mass );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for wave excitation forces
    target_signal   = "OutWex";
    read_signal     = _read_channel_value( infile, input->out_wex );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );
    
    //////////////////////////////////////////////
    /************** Site Conditions *************/
    //////////////////////////////////////////////

    // Skip header
    _skip_header( infile, line_count, 3 );

    // Read water density
    target_signal   = "RhoW";
    read_signal     = _read_channel_value( infile, input->water_density );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );
    
    // Read gravitational acceleration
    target_signal   = "GravAcc";
    read_signal     = _read_channel_value( infile, input->grav_acc );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read end simulation time
    target_signal   = "WaterDepth";
    read_signal     = _read_channel_value( infile, input->water_depth );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    //////////////////////////////////////////////
    /****************** Heading *****************/
    //////////////////////////////////////////////

    // Skip header
    _skip_header( infile, line_count, 3 );

    // Read headings units
    target_signal   = "HeadUnits";
    read_signal     = _read_channel_value( infile, input->heads_units );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read list contraction
    _read_list_contraction(
                                infile,
                                line_count,
                                target_file,
                                input->heads
                            );
    input->heads_np = input->heads.size( );

    //////////////////////////////////////////////
    /*************** Frequencies ****************/
    //////////////////////////////////////////////

    // Skip header
    _skip_header( infile, line_count, 3 );

    // Read headings units
    target_signal   = "FreqUnit";
    read_signal     = _read_channel_value( infile, input->freqs_unit );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read list contraction
    _read_list_contraction(
                                infile,
                                line_count,
                                target_file,
                                input->angfreqs
                            );
    input->angfreqs_np = input->angfreqs.size( );
 
    // Close file unit
    infile.close();
}


std::string _read_channel_name( std::ifstream& infile )
{
    // Generate local auxiliar variables
    std::string aux_str, line, channel_name;

    // Read line from file unit
    std::getline(infile, line);

    // Get signal value and name from the line
    // read
    std::istringstream iss(line);
    iss >> channel_name >> aux_str;

    return channel_name;
}


Input* read_input_files( std::string folder_path )
{
    // Instantiate an Input object
    Input* input = new Input();

    // Save case folder path
    input->case_fopath = folder_path;

    // Read case.input.dat file
    read_case( folder_path, input );

    // Read bodies
    read_bodies( folder_path, input );

    // Configure inputs
    input->configure( );

    return input;
}


void _skip_header( 
                    std::ifstream&  infile, 
                    int&            line_count, 
                    int             np 
                )
{
    std::string line;
    for ( int i=0; i<np; i++)
    {
        std::getline( infile,  line );
    }
    line_count += np;
}