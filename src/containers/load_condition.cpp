
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

// Include local modules
#include "../inout/reader.hpp"
#include "load_condition.hpp"
#include "../tools.hpp"
#include "../version.hpp"


        LoadCondition::LoadCondition( 
                                                std::string     fopath_in,
                                                std::string     finame_in
                                        )
{
    // Read load condition file
    this->_read_load_condition_file( 
                                        this->_fopath,
                                        this->_finame
                                    );
}

 void   LoadCondition::_read_load_condition_file( 
                                                        std::string    folder_path,
                                                        std::string    target_file
                                                )
{
    // Generate auxiliar variables
    int         line_count      = 0;
    std::string read_signal     = "";
    std::string target_signal   = "";
    std::string _version        = "";

    // Store file path and name
    this->_fopath               = folder_path;
    this->_finame               = target_file;

    // Generate case.input.dat file path
    std::filesystem::path folder_path_( folder_path );
    std::filesystem::path file_name( target_file );
    std::filesystem::path file_path = folder_path_ / file_name;

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

    // Read Is Fix
    target_signal   = "UseMass";
    read_signal     = read_channel_value( infile, this->use_mass );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    //////////////////////////////////////////////
    /************* Mass Properties **************/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read mass
    target_signal   = "BodyMass";
    read_signal     = read_channel_value( infile, this->mass );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read COGX
    target_signal   = "COGX";
    read_signal     = read_channel_value( infile, this->cog[0] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read COGY
    target_signal   = "COGY";
    read_signal     = read_channel_value( infile, this->cog[1] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read COGZ
    target_signal   = "COGZ";
    read_signal     = read_channel_value( infile, this->cog[2] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read switch to define how the intertia is defined, by the 
    // radius of gyrantion or by the inertial matrix components
    target_signal   = "IBR";
    read_signal     = read_channel_value( infile, this->interia_by_rad );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read IXX
    target_signal   = "IXX";
    read_signal     = read_channel_value( infile, this->inertia[0] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read IXY
    target_signal   = "IXY";
    read_signal     = read_channel_value( infile, this->inertia[1] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read IXZ
    target_signal   = "IXZ";
    read_signal     = read_channel_value( infile, this->inertia[2] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read IYY
    target_signal   = "IYY";
    read_signal     = read_channel_value( infile, this->inertia[3] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read IYZ
    target_signal   = "IYZ";
    read_signal     = read_channel_value( infile, this->inertia[4] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read IZZ
    target_signal   = "IZZ";
    read_signal     = read_channel_value( infile, this->inertia[5] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read RXX
    target_signal   = "RXX";
    read_signal     = read_channel_value( infile, this->rad_inertia[0] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read RYY
    target_signal   = "RYY";
    read_signal     = read_channel_value( infile, this->rad_inertia[1] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read RZZ
    target_signal   = "RZZ";
    read_signal     = read_channel_value( infile, this->rad_inertia[2] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    if ( this->interia_by_rad )
    {
        this->inertia[0] = this->mass * pow2s( this->rad_inertia[0] );
        this->inertia[3] = this->mass * pow2s( this->rad_inertia[1] );
        this->inertia[5] = this->mass * pow2s( this->rad_inertia[2] );
    }

}