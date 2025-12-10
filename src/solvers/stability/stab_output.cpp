
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


void    StabOutput::_create_hs_scalar_field( 
                                                H5::Group&          hs_gp,
                                                const char*         field_name
                                            )
{
    // Create dataset for waterplane area
    CREATE_DATASET( 
                        hs_gp,
                        field_name,
                        _DS_HS_S_NP,
                        this->_ds_hs_s,
                        cusfloat_h5
                    );
}


void    StabOutput::_create_hs_vector_field( 
                                                H5::Group&          hs_gp,
                                                const char*         field_name
                                            )
{
    // Create dataset for waterplane area
    CREATE_DATASET( 
                        hs_gp,
                        field_name,
                        _DS_HS_V_NP,
                        this->_ds_hs_v,
                        cusfloat_h5
                    );
}


void    StabOutput::_save_1D_input_data(
                                            H5::Group&              hs_gp,
                                            const char*             field_name,
                                            hsize_t*                dset_dims,
                                            std::vector<cusfloat>&  vec
                                        )
{
    CREATE_DATASET( 
                            hs_gp,
                            field_name,
                            _DS_1D_NP,
                            dset_dims,
                            cusfloat_h5
                    );

    hsize_t offset[_DS_1D_NP] = { 0 };
    SAVE_DATASET_CHUNK(
                            hs_gp,
                            field_name,
                            _DS_1D_NP,
                            dset_dims,
                            dset_dims,
                            offset,
                            vec.data( ),
                            cusfloat_h5
                        );
}


void    StabOutput::save_hydrostatics(
                                            int             axis_id,
                                            InitStabVec&    hydrostats
                                        )
{
    // Open file unit
    H5::H5File  fid( this->_results_fipath.c_str( ), H5F_ACC_RDWR );

    // Take the mesh group
    H5::Group   hs_gp   = fid.openGroup( _GN_HYDROSTATS );

    // Allocate space for body data chunks
    cusfloat*   amb_s   = generate_empty_vector<cusfloat>(
                                                            this->_input->heel_hs_rad.size( )
                                                            *
                                                            this->_input->draft_hs.size( )
                                                        );
    
    cusfloat*   amb_v   = generate_empty_vector<cusfloat>(
                                                            this->_input->heel_hs_rad.size( )
                                                            *
                                                            this->_input->draft_hs.size( )
                                                            *
                                                            3
                                                        );

    // Define list of scalar fields to be stored
    ScalarField scalar_fields[] =   {
                                        { _DN_AREA_WL,          [&]( std::size_t i ){ return hydrostats[i].get_area_wl( );      } },
                                        { _DN_AREA_IXX_WL,      [&]( std::size_t i ){ return hydrostats[i].get_area_wl_ixx( );  } },
                                        { _DN_AREA_IYY_WL,      [&]( std::size_t i ){ return hydrostats[i].get_area_wl_iyy( );  } },
                                        { _DN_AREA_MX_WL,       [&]( std::size_t i ){ return hydrostats[i].get_area_wl_mx( );   } },
                                        { _DN_AREA_MY_WL,       [&]( std::size_t i ){ return hydrostats[i].get_area_wl_my( );   } },
                                        { _DN_BMX,              [&]( std::size_t i ){ return hydrostats[i].get_bmx( );          } },
                                        { _DN_BMY,              [&]( std::size_t i ){ return hydrostats[i].get_bmy( );          } },
                                        { _DN_DISPLACEMENT,     [&]( std::size_t i ){ return hydrostats[i].get_mass( );         } },
                                        { _DN_GMX,              [&]( std::size_t i ){ return hydrostats[i].get_gmx( );          } },
                                        { _DN_GMY,              [&]( std::size_t i ){ return hydrostats[i].get_gmy( );          } },
                                        { _DN_KMX,              [&]( std::size_t i ){ return hydrostats[i].get_kmx( );          } },
                                        { _DN_KMY,              [&]( std::size_t i ){ return hydrostats[i].get_kmy( );          } },
                                        { _DN_VOLUME,           [&]( std::size_t i ){ return hydrostats[i].get_volume( );       } },
                                        { _DN_VOLUME_MX,        [&]( std::size_t i ){ return hydrostats[i].get_volume_mx( );    } },
                                        { _DN_VOLUME_MY,        [&]( std::size_t i ){ return hydrostats[i].get_volume_my( );    } }
                                    };

    // Loop over scalar fields to be stored
    for ( const auto& field : scalar_fields )
    {
        this->_save_hs_scalar_field(
                                        hs_gp,
                                        axis_id,
                                        field.name,
                                        amb_s,
                                        field.getter
                                    );
    }

    // Define list of scalar fields to be stored
    VectorField vector_fields[] =   {
                                        { _DN_COB,          [&]( std::size_t i ){ return hydrostats[i].get_cob( );      } },
                                    };

    // Loop over vector fields to be stored
    for ( const auto& field : vector_fields )
    {
        this->_save_hs_vector_field(
                                        hs_gp,
                                        axis_id,
                                        field.name,
                                        amb_v,
                                        field.getter
                                    );
    }

    // Close file unit
    fid.close( );

    // Deallocate heap memory
    mkl_free( amb_s );
    mkl_free( amb_v );

}


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

        // Create and storage input vector fields
        this->_save_1D_input_data( 
                                        hs_gp,
                                        _DN_DRAFT,
                                        this->_ds_hs_d,
                                        this->_input->draft_hs
                                    );

        this->_save_1D_input_data( 
                                        hs_gp,
                                        _DN_HEEL,
                                        this->_ds_hs_h,
                                        this->_input->heel_hs_deg
                                    );
        
        // Create scalar datasets for hydrostatics in the output file
        this->_create_hs_scalar_field(  hs_gp, _DN_AREA_WL          );
        this->_create_hs_scalar_field(  hs_gp, _DN_AREA_IXX_WL      );
        this->_create_hs_scalar_field(  hs_gp, _DN_AREA_IYY_WL      );
        this->_create_hs_scalar_field(  hs_gp, _DN_AREA_MX_WL       );
        this->_create_hs_scalar_field(  hs_gp, _DN_AREA_MY_WL       );
        this->_create_hs_scalar_field(  hs_gp, _DN_BMX              );
        this->_create_hs_scalar_field(  hs_gp, _DN_BMY              );
        this->_create_hs_scalar_field(  hs_gp, _DN_DISPLACEMENT     );
        this->_create_hs_scalar_field(  hs_gp, _DN_GMX              );
        this->_create_hs_scalar_field(  hs_gp, _DN_GMY              );
        this->_create_hs_scalar_field(  hs_gp, _DN_KMX              );
        this->_create_hs_scalar_field(  hs_gp, _DN_KMY              );
        this->_create_hs_scalar_field(  hs_gp, _DN_VOLUME           );
        this->_create_hs_scalar_field(  hs_gp, _DN_VOLUME_MX        );
        this->_create_hs_scalar_field(  hs_gp, _DN_VOLUME_MY        );
        this->_create_hs_scalar_field(  hs_gp, _DN_VOLUME_MZ        );
        
        // Create scalar datasets for hydrostatics in the output file
        this->_create_hs_vector_field(  hs_gp, _DN_COB              );

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
