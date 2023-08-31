
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
                    folder_path,
                    input->bodies_finame[i],
                    input->bodies[i]
                );

    }

    // Set bodies as read
    input->is_bodies = true;
}


void read_body(
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
    /************* Mesh Description *************/
    //////////////////////////////////////////////

    // Skip header
    _skip_header( infile, line_count, 3 );

    // Read mesh
    target_signal   = "MeshFile";
    read_signal     = _read_channel_value( infile, body->mesh_finame );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    fs::path mesh_foname_( std::string( "mesh" ) );
    fs::path mesh_finame_( body->mesh_finame );
    fs::path mesh_fipath_ = folder_path_ / mesh_foname_ / mesh_finame_;
    body->mesh      = new Mesh( mesh_fipath_.string( ) );
    body->is_mesh   = true;

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