
// Include general usage libraries
#include <filesystem>
#include <iostream>
#include <sstream>

// Include local modules
#include "input.hpp"
#include "reader.hpp"
#include "../tools.hpp"
#include "../version.hpp"

// Manage namespaces
namespace fs = std::filesystem;


// void read_body( Input* input )
// {
//     // Generate auxiliar variables
//     int     line_count      = 0;
//     string  read_signal     = "";
//     string  target_file     = "body_dynamics.input.dat";
//     string  target_signal   = "";
//     string  _version        = "";

//     // Compose file path of the target file
//     fs::path file_name( target_file );
//     fs::path file_path = input->plant_descr_path / file_name;

//     // Open file unit
//     ifstream infile;
//     infile.open( file_path, ios_base::in );
//     CHECK_FILE_UNIT_STATUS( infile, file_path );

//     // Read file header
//     _skip_header( infile, line_count, 1 );

//     // Read file version
//     target_signal   = "version";
//     read_signal     = _read_channel_value( infile, _version );
//     CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );
//     CHECK_INPUT_FILE_VERSION( VERSION_LABEL, _version, file_path );

//     //////////////////////////////////////////////
//     /************** Body Properties *************/
//     //////////////////////////////////////////////
//     _skip_header( infile, line_count, 3 );

//     // Read rigid body model: Rigid or flexible
//     target_signal   = "BodyModel";
//     read_signal     = _read_channel_value( infile, input->body_model );
//     CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

//     if ( input->body_model > 1)
//     {
//         cerr << "Unknown body model: " << input->body_model << endl;
//         cerr << ". Acepted values are-> 0: Rigid - 1: Flexible" << endl;
//         exit( 10 );
//     }
//     else if ( input->body_model == 0 )
//     {
//         input->body_dofs_np = 1;
//     }

//     //////////////////////////////////////////////
//     /********* Rigid Body properties ************/
//     //////////////////////////////////////////////
//     _skip_header( infile, line_count, 3 );

//     // Read rigid body COG X position
//     target_signal   = "RBCOGX";
//     read_signal     = _read_channel_value( infile, input->cog_x );
//     CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

//     // Read rigid body COG Y position
//     target_signal   = "RBCOGY";
//     read_signal     = _read_channel_value( infile, input->cog_y );
//     CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

//     // Read inertia matrix
//     target_signal   = "RigidInertiaMat";
//     read_signal     = _read_channel_name( infile );
//     CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );
//     _read_channel_matrix( infile, 3, 3, input->inertia );

//     // Read damping matrix
//     target_signal   = "RigidDampingMat";
//     read_signal     = _read_channel_name( infile );
//     CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );
//     _read_channel_matrix( infile, 3, 3, input->damping );

//     // Close file unit
//     infile.close();

// }


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