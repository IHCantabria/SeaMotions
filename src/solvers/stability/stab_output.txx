
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


// Include local modules
#include "stab_output.hpp"


template<typename T>
void    StabOutput::save_hydrostatics(
                                            int             axis_id,
                                            T&              hydrostats
                                    )
{
    // Open file unit
    H5::H5File  fid( this->_results_fipath.c_str( ), H5F_ACC_RDWR );

    // Take the mesh group
    H5::Group   hs_gp = fid.openGroup( _GN_HYDROSTATS );

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


template<typename Getter>
void StabOutput::_save_hs_scalar_field(
                                            H5::Group&  gp,
                                            int         axis_id,
                                            const char* ds_name,
                                            cusfloat*   buffer,
                                            Getter&&    getter
                                        )
{
    // Define datasets dimensions and offset indexes
    std::size_t heel_np                 = this->_input->heel_hs_rad.size( );
    std::size_t draft_np                = this->_input->draft_hs.size( );

    hsize_t _ds_hs_ch_s[_DS_HS_S_NP]    = { 1, static_cast<hsize_t>( heel_np ), static_cast<hsize_t>( draft_np ) };
    hsize_t offset_s[_DS_HS_S_NP]       = { static_cast<hsize_t>( axis_id ), 0, 0 };

    // Exract hydrostatics scalar field to a buffer
    _extract_hydrostats_scalar( 
                                    heel_np, 
                                    draft_np,
                                    [&]( std::size_t i ){ return getter( i ); },
                                    buffer
                                );

    // Save data into disk
    SAVE_DATASET_CHUNK(
                                    gp,
                                    ds_name,
                                    _DS_HS_S_NP,
                                    this->_ds_hs_s,
                                    _ds_hs_ch_s,
                                    offset_s,
                                    buffer,
                                    cusfloat_h5
                        );
    
}


template<typename Getter>
void StabOutput::_save_hs_vector_field(
                                            H5::Group&  gp,
                                            int         axis_id,
                                            const char* ds_name,
                                            cusfloat*   buffer,
                                            Getter&&    getter
                                        )
{
    // Define datasets dimensions and offset indexes
    std::size_t heel_np                 = this->_input->heel_hs_rad.size( );
    std::size_t draft_np                = this->_input->draft_hs.size( );

    hsize_t _ds_hs_ch_v[_DS_HS_V_NP]    = { 1, static_cast<hsize_t>( heel_np ), static_cast<hsize_t>( draft_np ), 3 };
    hsize_t offset_v[_DS_HS_V_NP]       = { static_cast<hsize_t>( axis_id ), 0, 0, 0 };

    // Exract hydrostatics scalar field to a buffer
    _extract_hydrostats_vector( 
                                    heel_np, 
                                    draft_np,
                                    [&]( std::size_t i ){ return getter( i ); },
                                    buffer
                                );

    // Save data into disk
    SAVE_DATASET_CHUNK(
                                    gp,
                                    ds_name,
                                    _DS_HS_V_NP,
                                    this->_ds_hs_v,
                                    _ds_hs_ch_v,
                                    offset_v,
                                    buffer,
                                    cusfloat_h5
                        );
    
}