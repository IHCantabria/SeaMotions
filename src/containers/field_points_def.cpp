
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
#include <filesystem>
#include <fstream>
#include <string>

// Include local modules
#include "field_points_def.hpp"
#include "../inout/reader.hpp"
#include "../tools.hpp"
#include "../version.hpp"


FieldMeshDataConfig FieldPointsDef::get_config( void ) const
{
    FieldMeshDataConfig config;
    config.mesh_file_path  = this->mesh_finame;
    config.body_name       = this->mesh_body_name;
    config.out_components  = this->out_components;
    config.out_potential   = this->out_potential;
    config.out_pressure    = this->out_pressure;
    config.out_velocity    = this->out_velocity;

    return config;

}


void FieldPointsDef::_load_mesh( 
                                    std::string fopath
                                )
{
    // Alias namespaces
    namespace fs = std::filesystem;

    // Load mesh for the current body
    fs::path mesh_foname_( std::string( "mesh" ) );
    fs::path mesh_finame_( this->mesh_finame );
    fs::path mesh_fipath_ = fopath / mesh_foname_ / mesh_finame_;

    this->mesh      = new Mesh( 
                                    mesh_fipath_.string( ),
                                    this->mesh_body_name,
                                    PanelTypeE::FIELD_POINT
                                );
    this->is_mesh   = true;

}


FieldPointsDef::FieldPointsDef( 
                                    std::string fopath,
                                    std::string finame
                                )
{
    // Parse input file
    this->_parse_input_file( fopath, finame );

    // Load mesh
    this->_load_mesh( fopath );

}


FieldPointsDef::~FieldPointsDef( void )
{
    if ( this->is_mesh )
    {
        delete this->mesh;
        this->mesh = nullptr;
    }
}


void FieldPointsDef::_parse_input_file( 
                                            std::string folder_path,
                                            std::string target_file
                                        )
{ 
    // Generate auxiliar variables
    int         line_count      = 0;
    std::string read_signal     = "";
    std::string target_signal   = "";
    std::string _version        = "";

    // Alias namespaces
    namespace fs = std::filesystem;

    // Generate case.input.dat file path
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
    /************* Output Channels **************/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read Output Components
    target_signal   = "OutComps";
    read_signal     = read_channel_value( infile, this->out_components );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read Output Potential
    target_signal   = "OutPot";
    read_signal     = read_channel_value( infile, this->out_potential );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read Output Pressure
    target_signal   = "OutPress";
    read_signal     = read_channel_value( infile, this->out_pressure );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read Output Velocity
    target_signal   = "OutVel";
    read_signal     = read_channel_value( infile, this->out_velocity );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    //////////////////////////////////////////////
    /************* Mesh Description *************/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read mesh body name
    target_signal   = "BodyName";
    read_signal     = read_channel_value( infile, this->mesh_body_name );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read mesh file name
    target_signal   = "MeshFile";
    read_signal     = read_channel_value( infile, this->mesh_finame );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Close file unit
    infile.close();
}