
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

// Include local modules
#include "morison_element_def.hpp"
#include "../inout/reader.hpp"
#include "../math/math_tools.hpp"
#include "../tools.hpp"
#include "../version.hpp"


MorisonElementDef::MorisonElementDef( 
                                        std::ifstream& file,
                                        cusfloat*      cog 
                                    )
{
    // Storage of the center of gravity of the element
    this->cog[0] = cog[0];
    this->cog[1] = cog[1];
    this->cog[2] = cog[2];
    
    // Parse input file to fill the rest of the element properties
    this->_parse_input_file( file );

}


MorisonElementDef::MorisonElementDef( 
                                        const std::string& fopath,
                                        const std::string& finame,
                                        cusfloat*          cog
                                    )
{
    // Storage of the center of gravity of the element
    this->cog[0] = cog[0];
    this->cog[1] = cog[1];
    this->cog[2] = cog[2];


    // Parse input file to fill the rest of the element properties
    this->_parse_input_file( fopath, finame );

}


void MorisonElementDef::_parse_input_file( 
                                                std::ifstream& file 
                                            )
{
    // Read start position of the element
    file >> this->start_pos[0] >> this->start_pos[1] >> this->start_pos[2];

    // Read azimuth angle of the element
    file >> this->azimuth_angle;

    // Read declination angle of the element
    file >> this->declination_angle;

    // Read rotation angle along the direction axis
    file >> this->rotation_angle;

    // Read length of the element
    file >> this->length;

    // Read number of divisons of the element
    file >> this->divisions_np;

    // Read Keulegan-Carpenter characteristic length
    file >> this->kc_length;

    // Read axial area for inertial forces calculation
    file >> this->area_axial_ca;

    // Read axial area for drag calculation
    file >> this->area_axial_cd;

    // Read transversal area for drag calculation
    file >> this->area_transversal_cd;

    // Read vertical area for drag calculation
    file >> this->area_vertical_cd;

    // Read axis added mass coefficient file name
    file >> this->axis_added_mass_coeff_file;

    // Read transversal added mass coefficient file name
    file >> this->transversal_added_mass_coeff_file;

    // Read vertical added mass coefficient file name
    file >> this->vertical_added_mass_coeff_file;

    // Read axis drag coefficient file name
    file >> this->axis_drag_coeff_file;

    // Read transversal drag coefficient file name
    file >> this->transversal_drag_coeff_file;

    // Read vertical drag coefficient file name
    file >> this->vertical_drag_coeff_file;

}

void MorisonElementDef::_parse_input_file( 
                                                const std::string& fopath,
                                                const std::string& target_file
                                            )
{
    // Generate auxiliar variables
    int         line_count      = 0;
    std::string read_signal     = "";
    std::string target_signal   = "";
    std::string _version        = "";

    // Compose file path
    namespace   fs          = std::filesystem;
    fs::path    file_path_  = fs::path( fopath ) / fs::path( target_file );
    std::string file_path   = file_path_.string( );


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

    // Read Start Position X
    target_signal   = "START_POINT_X";
    read_signal     = read_channel_value( infile, this->start_pos[0] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read Start Position Y
    target_signal   = "START_POINT_Y";
    read_signal     = read_channel_value( infile, this->start_pos[1] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read Start Position Z
    target_signal   = "START_POINT_Z";
    read_signal     = read_channel_value( infile, this->start_pos[2] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );


    //////////////////////////////////////////////
    /******* Orientation of the Element *********/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read Azimuth angle
    target_signal   = "AZIMUTH";
    read_signal     = read_channel_value( infile, this->azimuth_angle );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );
    this->azimuth_angle = deg_to_rad( this->azimuth_angle );

    // Read Declination angle
    target_signal   = "DECLINATION";
    read_signal     = read_channel_value( infile, this->declination_angle );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );
    this->declination_angle = deg_to_rad( this->declination_angle );

    // Read Direction axis rotation angle
    target_signal   = "DIR_AXIS_ROT";
    read_signal     = read_channel_value( infile, this->rotation_angle );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );
    this->rotation_angle = deg_to_rad( this->rotation_angle );


    //////////////////////////////////////////////
    /********* Element discretization ************/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read Element length
    target_signal   = "LENGTH";
    read_signal     = read_channel_value( infile, this->length );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read number of divisions
    target_signal   = "NUM_DIVISIONS";
    read_signal     = read_channel_value( infile, this->divisions_np );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read Keulegan-Carpenter characteristic length
    target_signal   = "KC_LENGTH";
    read_signal     = read_channel_value( infile, this->kc_length );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );


    //////////////////////////////////////////////
    /********** Hydrodynamic Areas ***************/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read axial area for inertial forces calculation
    target_signal   = "AXIAL_AREA_CA";
    read_signal     = read_channel_value( infile, this->area_axial_ca );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read axial area for drag calculation
    target_signal   = "AXIAL_AREA_CD";
    read_signal     = read_channel_value( infile, this->area_axial_cd );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read transverse area for drag calculation
    target_signal   = "TRANS_AREA_CD";
    read_signal     = read_channel_value( infile, this->area_transversal_cd );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read vertical area for drag calculation
    target_signal   = "VERT_AREA_CD";
    read_signal     = read_channel_value( infile, this->area_vertical_cd );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );


    //////////////////////////////////////////////
    /******* Hydrodynamic Coefficients ***********/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read axial added-mass coefficient file
    target_signal   = "AXIAL_CA_FILE";
    read_signal     = read_channel_value( infile, this->axis_added_mass_coeff_file );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read transverse added-mass coefficient file
    target_signal   = "TRANS_CA_FILE";
    read_signal     = read_channel_value( infile, this->transversal_added_mass_coeff_file );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read vertical added-mass coefficient file
    target_signal   = "VERT_CA_FILE";
    read_signal     = read_channel_value( infile, this->vertical_added_mass_coeff_file );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read axial drag coefficient file
    target_signal   = "AXIAL_CD_FILE";
    read_signal     = read_channel_value( infile, this->axis_drag_coeff_file );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read transverse drag coefficient file
    target_signal   = "TRANS_CD_FILE";
    read_signal     = read_channel_value( infile, this->transversal_drag_coeff_file );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

    // Read vertical drag coefficient file
    target_signal   = "VERT_CD_FILE";
    read_signal     = read_channel_value( infile, this->vertical_drag_coeff_file );
    CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );


    // Close file
    infile.close( );
}