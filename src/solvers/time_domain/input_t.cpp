
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

// Include general usage libraries
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

// Include local modules
#include "input_t.hpp"
#include "../../inout/reader.hpp"
#include "../../math/math_tools.hpp"
#include "../../tools.hpp"
#include "../../version.hpp"

// Manage namespaces
namespace fs = std::filesystem;


InputT::~InputT( )
{
    if ( this->bodies != nullptr )
    {
        for ( int i=0; i<this->bodies_np; i++ )
        {
            delete this->bodies[i];
        }
        delete [] this->bodies;
        this->bodies = nullptr;
    }
}


void InputT::load( const std::string& folder_path )
{
    // Save case folder path
    this->case_fopath = folder_path;
    this->folder_path = folder_path;

    // Read case.input.dat file
    this->read_case( folder_path );

    // Read bodies
    this->read_bodies( folder_path );

    // Initialize class
    this->initialize( );
}


void InputT::read_bodies( const std::string& folder_path )
{
    // Allocate space for the bodies description
    this->bodies = new BodyDef*[this->bodies_np];

    // Loop over body definitions to load them
    for ( int i=0; i<this->bodies_np; i++ )
    {
        this->bodies[i] = new BodyDef;

        read_body(
                    folder_path,
                    this->bodies_finame[i],
                    this->bodies[i]
                );
    }
}


void InputT::read_body(
                            const std::string& folder_path,
                            const std::string& target_file,
                            BodyDef* body
                      )
{
    // Generate auxiliar variables
    int         line_count      = 0;
    std::string read_signal     = "";
    std::string target_signal   = "";
    std::string _version        = "";

    // Generate body file path
    fs::path folder_path_( folder_path );
    fs::path file_name( target_file );
    fs::path file_path = folder_path_ / file_name;

    // Open file unit
    std::ifstream infile;
    infile.open( file_path, std::ios_base::in );
    CHECK_FILE_UNIT_STATUS( infile, file_path );

    // Read file header line
    skip_header( infile, line_count, 1 );

    // Read file version
    target_signal   = "version";
    read_signal     = read_channel_value( infile, _version );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );
    CHECK_INPUT_FILE_VERSION( VERSION_LABEL, _version, file_path );

    //////////////////////////////////////////////
    /*********** General Properties *************/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read per-DOF kinematic restrictions (0 = free, 1 = fixed)
    target_signal   = "RSurge";
    read_signal     = read_channel_value( infile, body->dof_restrictions[0] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    target_signal   = "RSway";
    read_signal     = read_channel_value( infile, body->dof_restrictions[1] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    target_signal   = "RHeave";
    read_signal     = read_channel_value( infile, body->dof_restrictions[2] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    target_signal   = "RRoll";
    read_signal     = read_channel_value( infile, body->dof_restrictions[3] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    target_signal   = "RPitch";
    read_signal     = read_channel_value( infile, body->dof_restrictions[4] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    target_signal   = "RYaw";
    read_signal     = read_channel_value( infile, body->dof_restrictions[5] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    //////////////////////////////////////////////
    /************* Mass Properties **************/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read mass
    target_signal   = "BodyMass";
    read_signal     = read_channel_value( infile, body->mass );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read COGX
    target_signal   = "COGX";
    read_signal     = read_channel_value( infile, body->cog[0] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read COGY
    target_signal   = "COGY";
    read_signal     = read_channel_value( infile, body->cog[1] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read COGZ
    target_signal   = "COGZ";
    read_signal     = read_channel_value( infile, body->cog[2] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read switch to define how the inertia is defined
    target_signal   = "IBR";
    read_signal     = read_channel_value( infile, body->interia_by_rad );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read IXX
    target_signal   = "IXX";
    read_signal     = read_channel_value( infile, body->inertia[0] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read IXY
    target_signal   = "IXY";
    read_signal     = read_channel_value( infile, body->inertia[1] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read IXZ
    target_signal   = "IXZ";
    read_signal     = read_channel_value( infile, body->inertia[2] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read IYY
    target_signal   = "IYY";
    read_signal     = read_channel_value( infile, body->inertia[3] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read IYZ
    target_signal   = "IYZ";
    read_signal     = read_channel_value( infile, body->inertia[4] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read IZZ
    target_signal   = "IZZ";
    read_signal     = read_channel_value( infile, body->inertia[5] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read RXX
    target_signal   = "RXX";
    read_signal     = read_channel_value( infile, body->rad_inertia[0] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read RYY
    target_signal   = "RYY";
    read_signal     = read_channel_value( infile, body->rad_inertia[1] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read RZZ
    target_signal   = "RZZ";
    read_signal     = read_channel_value( infile, body->rad_inertia[2] );
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
    skip_header( infile, line_count, 3 );

    // Read mesh body name
    target_signal   = "BodyName";
    read_signal     = read_channel_value( infile, body->mesh_body_name );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read mesh file name
    target_signal   = "MeshFile";
    read_signal     = read_channel_value( infile, body->mesh_finame );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read usage of internal lid
    target_signal   = "LidType";
    read_signal     = read_channel_value( infile, body->lid_type );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    //////////////////////////////////////////////
    /********** Morison Coefficients ************/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read Morison switch
    target_signal   = "UseMC";
    read_signal     = read_channel_value( infile, body->use_morison_elements );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    body->morison_elements_names.clear( );
    body->morison_elements_np = 0;

    // Load mesh for the current body
    std::vector<Mesh*> _mesh_total;
    fs::path mesh_foname_{ MESH_FOLDER_NAME };
    fs::path mesh_finame_{ body->mesh_finame };
    fs::path mesh_fipath_ = folder_path_ / mesh_foname_ / mesh_finame_;

    body->mesh          = new Mesh(
                                        mesh_fipath_.string( ),
                                        body->mesh_body_name,
                                        body->cog,
                                        body->is_fix,
                                        PanelTypeE::DIFFRAC
                                    );
    body->is_mesh       = true;
    body->mesh_items_np = 1;
    _mesh_total.push_back( body->mesh );

    // Write mesh
    fs::path mesh_foname_1_{ RESULTS_MESH_FOLDER_NAME };
    fs::path results_foname_{ RESULTS_FOLDER_NAME };
    fs::path results_mesh_fipath_ = folder_path_ / results_foname_;
    fs::path plot_mesh_fipath_ = results_mesh_fipath_ / mesh_foname_1_;

    if ( !fs::exists( results_mesh_fipath_ ) )
    {
        fs::create_directory( results_mesh_fipath_ );
    }

    if ( !fs::exists( plot_mesh_fipath_ ) )
    {
        fs::create_directory( plot_mesh_fipath_ );
    }

    body->mesh->write( plot_mesh_fipath_.string( ) );

    // Load lid if required
    if ( body->lid_type == 1 )
    {
        std::stringstream ss;
        ss << body->mesh_body_name << "_int_lid";

        body->mesh_int_lid  = new Mesh(
                                            mesh_fipath_.string( ),
                                            ss.str( ),
                                            body->cog,
                                            body->is_fix,
                                            PanelTypeE::INT_LID
                                        );
        body->is_mesh_int_lid   = true;
        body->mesh_items_np++;
        _mesh_total.push_back( body->mesh_int_lid );
    }

    // QTF free surface is NOT loaded in the time domain solver

    // Compose total mesh
    body->mesh_total    = new Mesh(
                                        body->mesh->name,
                                        _mesh_total,
                                        body->cog,
                                        body->is_fix
                                    );
    body->is_mesh_total = true;

    body->mesh_total->write( plot_mesh_fipath_.string( ), "total" );

    // Close file unit
    infile.close( );
}


void InputT::read_case( const std::string& folder_path )
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
    skip_header( infile, line_count, 1 );

    // Read file version
    target_signal   = "version";
    read_signal     = read_channel_value( infile, _version );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );
    CHECK_INPUT_FILE_VERSION( VERSION_LABEL, _version, file_path );

    //////////////////////////////////////////////
    /************** Solver Controls *************/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read total simulation time
    target_signal   = "SimTime";
    read_signal     = read_channel_value( infile, this->sim_time );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read time step
    target_signal   = "Dt";
    read_signal     = read_channel_value( infile, this->dt );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read Duhamel (radiation memory) integral switch
    target_signal   = "IsMemEff";
    read_signal     = read_channel_value( infile, this->use_duhamel );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    //////////////////////////////////////////////
    /************** Body Definition *************/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read bodies name definition
    target_signal = "BodyFN";
    read_channel_list(
                            infile,
                            target_file,
                            target_signal,
                            line_count,
                            this->bodies_finame
                        );
    this->bodies_np = static_cast<int>( this->bodies_finame.size( ) );

    //////////////////////////////////////////////
    /********* Field Points Definition **********/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read use field points flag
    target_signal   = "UseFP";
    read_signal     = read_channel_value( infile, this->use_field_points );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read field points file name
    target_signal   = "FieldPN";
    read_signal     = read_channel_value( infile, this->field_points_finame );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    if ( !this->use_field_points )
    {
        this->field_points_finame.clear( );
    }

    //////////////////////////////////////////////
    /************** Output Channels *************/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read output flag for Froude-Krylov force
    target_signal   = "OutFK";
    read_signal     = read_channel_value( infile, this->out_fk );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for hydrodynamic force
    target_signal   = "OutFH";
    read_signal     = read_channel_value( infile, this->out_fh );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for total force
    target_signal   = "OutFT";
    read_signal     = read_channel_value( infile, this->out_ft );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for potential
    target_signal   = "OutPot";
    read_signal     = read_channel_value( infile, this->out_potential );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for pressure
    target_signal   = "OutPress";
    read_signal     = read_channel_value( infile, this->out_pressure );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for mesh
    target_signal   = "OutMesh";
    read_signal     = read_channel_value( infile, this->out_mesh );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for sources
    target_signal   = "OutSources";
    read_signal     = read_channel_value( infile, this->out_sources );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read output flag for structural mass
    target_signal   = "OutStMass";
    read_signal     = read_channel_value( infile, this->out_struct_mass );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    //////////////////////////////////////////////
    /************** Site Conditions *************/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read water density
    target_signal   = "RhoW";
    read_signal     = read_channel_value( infile, this->water_density );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read gravitational acceleration
    target_signal   = "GravAcc";
    read_signal     = read_channel_value( infile, this->grav_acc );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read water depth
    target_signal   = "WaterDepth";
    read_signal     = read_channel_value( infile, this->water_depth );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read wave amplitude
    target_signal   = "WaveA";
    read_signal     = read_channel_value( infile, this->wave_amp );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read wave period
    target_signal   = "WaveP";
    read_signal     = read_channel_value( infile, this->wave_period );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read wave heading
    target_signal   = "WaveH";
    read_signal     = read_channel_value( infile, this->wave_heading );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Close file unit
    infile.close( );
}


void InputT::initialize( )
{
    this->ang_freq = ( this->wave_period > static_cast<cusfloat>( 0.0 ) )
                        ? ( static_cast<cusfloat>( 2.0 ) * PI / this->wave_period )
                        : static_cast<cusfloat>( 0.0 );
}
