
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
#include "../../containers/logger.hpp"
#include "../../inout/reader.hpp"
#include "stab_input.hpp"
#include "../../tools.hpp"
#include "../../version.hpp"


void    StabInput::_initialize(
                                    void
                                )
{
    // Compose mesh file path
    std::filesystem::path fopath_ ( this->_fopath );
    std::filesystem::path mesh_   ( "mesh" );
    std::filesystem::path finame_ ( this->mesh_finame );
    std::filesystem::path fipath_ = fopath_ / mesh_ / finame_;

    this->mesh_fipath = fipath_.string( );

    // Convert hydrostatic and gz curves heelings to 
    // radians
    this->heel_gz_rad.resize( this->heel_gz_deg.size( ) );
    for ( std::size_t i=0; i<this->heel_hs_deg.size( ); i++ )
    {
        this->heel_gz_rad[i] = deg_to_rad( this->heel_gz_deg[i] );
    }

    this->heel_hs_rad.resize( this->heel_hs_deg.size( ) );
    for ( std::size_t i=0; i<this->heel_hs_deg.size( ); i++ )
    {
        this->heel_hs_rad[i] = deg_to_rad( this->heel_hs_deg[i] );
    }

}


void StabInput::_read_input_file( 
                                    std::string fopath 
                                )
{
    // Generate auxiliar variables
    int         line_count      = 0;
    std::string read_signal     = "";
    std::string target_signal   = "";
    std::string _version        = "";

    // Generate case.input.dat file path
    std::filesystem::path fopath_( fopath );
    std::filesystem::path finame_( this->_finame );
    std::filesystem::path fipath = fopath_ / finame_;

    // Storage folder and file locations
    this->_fopath = fopath;
    this->_fipath = fipath.string( );

    // Open file unit
    std::ifstream infile;
    infile.open( fipath, std::ios_base::in );
    CHECK_FILE_UNIT_STATUS( infile, fipath );

    // Read file header line
    skip_header( infile, line_count, 1 );

    // Read file version
    target_signal   = "version";
    read_signal     = read_channel_value( infile, _version );
    CHECK_SIGNAL_NAME( read_signal, target_signal, this->_finame, line_count );
    CHECK_INPUT_FILE_VERSION( VERSION_LABEL, _version, fipath );

    //////////////////////////////////////////////
    /************** Site Conditions *************/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read water density
    target_signal   = "RhoW";
    read_signal     = read_channel_value( infile, this->water_density );
    CHECK_SIGNAL_NAME( read_signal, target_signal, this->_finame, line_count );
    
    // Read gravitational acceleration
    target_signal   = "GravAcc";
    read_signal     = read_channel_value( infile, this->grav_acc );
    CHECK_SIGNAL_NAME( read_signal, target_signal, this->_finame, line_count );

    //////////////////////////////////////////////
    /************** Body Definition *************/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read mass
    target_signal   = "BodyMass";
    read_signal     = read_channel_value( infile, this->mass );
    CHECK_SIGNAL_NAME( read_signal, target_signal, this->_finame, line_count );

    // Read COGX
    target_signal   = "COGX";
    read_signal     = read_channel_value( infile, this->cog[0] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, this->_finame, line_count );

    // Read COGY
    target_signal   = "COGY";
    read_signal     = read_channel_value( infile, this->cog[1] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, this->_finame, line_count );

    // Read COGZ
    target_signal   = "COGZ";
    read_signal     = read_channel_value( infile, this->cog[2] );
    CHECK_SIGNAL_NAME( read_signal, target_signal, this->_finame, line_count );

    // Read Mesh file nanme
    target_signal   = "BodyName";
    read_signal     = read_channel_value( infile, this->body_name );
    CHECK_SIGNAL_NAME( read_signal, target_signal, this->_finame, line_count );

    // Read Mesh file nanme
    target_signal   = "MeshFile";
    read_signal     = read_channel_value( infile, this->mesh_finame );
    CHECK_SIGNAL_NAME( read_signal, target_signal, this->_finame, line_count );

    //////////////////////////////////////////////
    /************ Output Channels ***************/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read output flag hydrostatics
    target_signal   = "OutEQ";
    read_signal     = read_channel_value( infile, this->out_eq );
    CHECK_SIGNAL_NAME( read_signal, target_signal, this->_finame, line_count );

    // Read output flag hydrostatics
    target_signal   = "OutGZ";
    read_signal     = read_channel_value( infile, this->out_gz );
    CHECK_SIGNAL_NAME( read_signal, target_signal, this->_finame, line_count );

    // Read output flag hydrostatics
    target_signal   = "OutHS";
    read_signal     = read_channel_value( infile, this->out_hs );
    CHECK_SIGNAL_NAME( read_signal, target_signal, this->_finame, line_count );

    //////////////////////////////////////////////
    /**************** GZ Curves *****************/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read hydrostatics draft
    target_signal   = "DRAFTGZ";
    read_signal     = read_channel_name( infile );
    CHECK_SIGNAL_NAME( read_signal, target_signal, this->_finame, line_count );

    read_list_contraction(
                                infile,
                                line_count,
                                this->_finame,
                                this->draft_gz
                            );

    // Read hydrostatics heels
    target_signal   = "HEELGZ";
    read_signal     = read_channel_name( infile );
    CHECK_SIGNAL_NAME( read_signal, target_signal, this->_finame, line_count );

    read_list_contraction(
                                infile,
                                line_count,
                                this->_finame,
                                this->heel_gz_deg
                            );

    

    //////////////////////////////////////////////
    /*************** Hydrostaics ****************/
    //////////////////////////////////////////////

    // Skip header
    skip_header( infile, line_count, 3 );

    // Read hydrostatics draft
    target_signal   = "DRAFTHS";
    read_signal     = read_channel_name( infile );
    CHECK_SIGNAL_NAME( read_signal, target_signal, this->_finame, line_count );

    read_list_contraction(
                                infile,
                                line_count,
                                this->_finame,
                                this->draft_hs
                            );

    // Read hydrostatics heels
    target_signal   = "HEELHS";
    read_signal     = read_channel_name( infile );
    CHECK_SIGNAL_NAME( read_signal, target_signal, this->_finame, line_count );

    read_list_contraction(
                                infile,
                                line_count,
                                this->_finame,
                                this->heel_hs_deg
                            );

    // Close file unit
    infile.close( );

}


StabInput::StabInput( 
                        std::string fopath
                    )
{
    // Read input file
    this->_read_input_file( fopath );

    // Initialize class to configure input
    this->_initialize( );

}