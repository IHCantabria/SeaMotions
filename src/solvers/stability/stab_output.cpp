
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

// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "stab_output.hpp"
#include "../../tools.hpp"
#include "../../version.hpp"


StabOutput::StabOutput( 
                                            StabInput*      input
                )
{
    // Storage pointer to the input class instance
    this->_input = input;

    // Generate output results file path
    std::filesystem::path _case_fopath( input->get_case_fopath( ) );
    std::filesystem::path _results_finame( std::string( "results.stab.h5" ) );
    std::filesystem::path _results_fipath = _case_fopath / _results_finame;

    this->_results_fipath = _results_fipath.string( );


    /******************************************************/
    /************ Define datasets dimensions **************/
    /******************************************************/

    // Storage GZ heelings dataset dimensions
    this->_ds_gz_h[0]   = input->heel_gz_rad.size( );

    // Storage GZ curves dataset dimensions
    this->_ds_gz[0]     = 2;
    this->_ds_gz[1]     = input->heel_gz_rad.size( );

    // Storage Hydrostatics draft dataset dimensions
    this->_ds_hs_d[0]   = input->draft_hs.size( );

    // Storage Hydrostatics heeling dataset dimensions
    this->_ds_hs_h[0]   = input->heel_hs_rad.size( );

    // Storage Hydrostatic dataset dimensions
    this->_ds_hs_s[0]   = 2;
    this->_ds_hs_s[1]   = input->heel_hs_rad.size( );
    this->_ds_hs_s[2]   = input->draft_hs.size( );

    this->_ds_hs_v[0]   = 2;
    this->_ds_hs_v[1]   = input->heel_hs_rad.size( );
    this->_ds_hs_v[2]   = input->draft_hs.size( );
    this->_ds_hs_v[3]   = 3;

    /******************************************************/
    /**** Create datasets for the required output data ****/
    /******************************************************/

    // Open file unit
    H5::H5File fid( this->_results_fipath.c_str( ), H5F_ACC_TRUNC );

    CREATE_ATTRIBUTE( fid, "version_major", VERSION_MAJOR, H5::PredType::NATIVE_INT );
    CREATE_ATTRIBUTE( fid, "version_minor", VERSION_MINOR, H5::PredType::NATIVE_INT );
    CREATE_ATTRIBUTE( fid, "version_patch", VERSION_PATCH, H5::PredType::NATIVE_INT );

    if ( input->out_hs )
    {
        // Create group for hydrostatics
        H5::Group hs_gp( fid.createGroup( _GN_HYDROSTATS ) );

        // Create dataset for hydrostatics drafts and storage the data
        CREATE_DATASET( 
                            hs_gp,
                            _DN_DRAFT,
                            _DS_1D_NP,
                            this->_ds_hs_d,
                            cusfloat_h5
                        );

        hsize_t offset_hs_drafts[_DS_1D_NP] = { 0 };
        SAVE_DATASET_CHUNK(
                            hs_gp,
                            _DN_DRAFT,
                            _DS_1D_NP,
                            this->_ds_hs_d,
                            this->_ds_hs_d,
                            offset_hs_drafts,
                            this->_input->draft_hs.data( ),
                            cusfloat_h5
                        );

        // Create dataset for hydrostatics heels and storage the data
        CREATE_DATASET( 
                            hs_gp,
                            _DN_HEEL,
                            _DS_1D_NP,
                            this->_ds_hs_h,
                            cusfloat_h5
                        );

        hsize_t offset_hs_heels[_DS_1D_NP] = { 0 };
        SAVE_DATASET_CHUNK(
                            hs_gp,
                            _DN_HEEL,
                            _DS_1D_NP,
                            this->_ds_hs_h,
                            this->_ds_hs_h,
                            offset_hs_heels,
                            this->_input->heel_hs_rad.data( ),
                            cusfloat_h5
                        );

        // Create dataset for waterplane area
        CREATE_DATASET( 
                            hs_gp,
                            _DN_AREA_WL,
                            _DS_HS_S_NP,
                            this->_ds_hs_s,
                            cusfloat_h5
                        );

        // Create dataset for waterplane area first order moment w.r.t X axis
        CREATE_DATASET( 
                            hs_gp,
                            _DN_AREA_MX_WL,
                            _DS_HS_S_NP,
                            this->_ds_hs_s,
                            cusfloat_h5
                        );
        
        // Create dataset for waterplane area first order moment w.r.t Y axis
        CREATE_DATASET( 
                            hs_gp,
                            _DN_AREA_MY_WL,
                            _DS_HS_S_NP,
                            this->_ds_hs_s,
                            cusfloat_h5
                        );

        // Create dataset for waterplane area second order moment w.r.t X axis
        CREATE_DATASET( 
                            hs_gp,
                            _DN_AREA_IXX_WL,
                            _DS_HS_S_NP,
                            this->_ds_hs_s,
                            cusfloat_h5
                        );

        // Create dataset for waterplane area second order moment w.r.t Y axis
        CREATE_DATASET( 
                            hs_gp,
                            _DN_AREA_IYY_WL,
                            _DS_HS_S_NP,
                            this->_ds_hs_s,
                            cusfloat_h5
                        );

        // Create dataset for position of the center of buoyancy
        CREATE_DATASET( 
                            hs_gp,
                            _DN_COB,
                            _DS_HS_V_NP,
                            this->_ds_hs_v,
                            cusfloat_h5
                        );

        // Create dataset for hull displacement
        CREATE_DATASET( 
                            hs_gp,
                            _DN_DISPLACEMENT,
                            _DS_HS_S_NP,
                            this->_ds_hs_s,
                            cusfloat_h5
                        );

        // Create dataset for hull underwater volume
        CREATE_DATASET( 
                            hs_gp,
                            _DN_VOLUME,
                            _DS_HS_S_NP,
                            this->_ds_hs_s,
                            cusfloat_h5
                        );

        // Create dataset for hull underwater volume first order moment w.r.t X axis
        CREATE_DATASET( 
                            hs_gp,
                            _DN_VOLUME_MX,
                            _DS_HS_S_NP,
                            this->_ds_hs_s,
                            cusfloat_h5
                        );

        // Create dataset for hull underwater volume first order moment w.r.t Y axis
        CREATE_DATASET( 
                            hs_gp,
                            _DN_VOLUME_MY,
                            _DS_HS_S_NP,
                            this->_ds_hs_s,
                            cusfloat_h5
                        );
        
    }

    if ( input->out_gz )
    {
        // Create group for gz curves
        H5::Group gz_gp( fid.createGroup( _GN_GZ ) );

        // Create dataset for GZ curves heelings and storage the data
        CREATE_DATASET( 
                            gz_gp,
                            _DN_HEEL,
                            _DS_1D_NP,
                            this->_ds_gz_h,
                            cusfloat_h5
                        );

        hsize_t offset_gz_heels[_DS_1D_NP] = { 0 };
        SAVE_DATASET_CHUNK(
                            gz_gp,
                            _DN_HEEL,
                            _DS_1D_NP,
                            this->_ds_gz_h,
                            this->_ds_gz_h,
                            offset_gz_heels,
                            this->_input->heel_gz_rad.data( ),
                            cusfloat_h5
                        );

        // Create dataset for GZ curves
        CREATE_DATASET( 
                            gz_gp,
                            _DN_GZ,
                            _DS_GZ_NP,
                            this->_ds_gz,
                            cusfloat_h5
                        );
    }

    // Close file unit
    fid.close( );

}
