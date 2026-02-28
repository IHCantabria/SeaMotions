
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
#include "rad_diff_data.hpp"
#include "../mesh/panel_set_view.hpp"
#include <complex>
#include <type_traits>


template<typename T, typename Config>
template<FieldTypeE field_type, FieldComponentE field_component>
const cut::CusTensor<cusfloat>* RadDiffData<T, Config>::get_field_data( std::size_t heading_index ) const
{
    std::size_t count           = 0;
    std::size_t offset          = 0;
    const cut::CusTensor<T>*    field_data_l    = nullptr;
    cusfloat* out_data = this->_field_data.data( );
    for ( std::size_t j=0; j<this->_size_local; j++ )
    {
        offset          = heading_index * this->panel_data[j].field_points_np;
        field_data_l    = this->panel_data[j].template get_field_data<field_type, field_component>( );
        const T* field_data_ptr = field_data_l->data( );
        for ( std::size_t k=0; k<this->panel_data[j].field_points_np; k++ )
        {
            if constexpr ( std::is_same<T, cuscomplex>::value )
            {
                out_data[ count ] = static_cast<cusfloat>( field_data_ptr[ offset + k ].real( ) );
            }
            else
            {
                out_data[ count ] = static_cast<cusfloat>( field_data_ptr[ offset + k ] );
            }
            count++;
        }
    }

    return &(this->_field_data);
}

template<typename T, typename Config>
std::size_t RadDiffData<T, Config>::get_end_pos( void ) const
{
    return this->_end_pos;
}


template<typename T, typename Config>
std::size_t RadDiffData<T, Config>::get_size_global( void ) const
{
    return this->_size_global;
}


template<typename T, typename Config>
std::size_t RadDiffData<T, Config>::get_size_local( void ) const
{
    return this->_size_local;
}


template<typename T, typename Config>
std::size_t RadDiffData<T, Config>::get_size_local_fp( void ) const
{
    return this->_size_local_fp;
}


template<typename T, typename Config>
std::size_t RadDiffData<T, Config>::get_start_pos( void ) const
{
    return this->_start_pos;
}


template<typename T, typename Config>
RadDiffData<T, Config>::RadDiffData( 
                                        MpiConfig*      mpi_config_,
                                        Mesh*           mesh_,
                                        std::size_t     freqs_np_,
                                        std::size_t     headings_np_,
                                        std::size_t     dofs_np_,
                                        bool            use_mesh_nodes_,
                                        bool            use_waterline_
                                    )
{
    // Load field points coordinates from mesh nodes if specified, otherwise they will be loaded from panel geometry in the PanelData constructor
    if ( use_mesh_nodes_ )
    {
        // Store number of field points
        this->_size_global      = mesh_->nodes_np;
        this->_mpi_config       = mpi_config_;

        // Calculate start and end positions for the current process
        this->_mpi_config->get_1d_bounds(
                                            static_cast<int>( this->_size_global ),
                                            reinterpret_cast<int&>( this->_start_pos ),
                                            reinterpret_cast<int&>( this->_end_pos )
                                        );

        // Determine size of the local data chunk to storage the field points values for the differen degrees of freedom and headings
        this->_size_local       = this->_end_pos - this->_start_pos;
        this->_size_local_fp    = this->_size_local;

        // Allocate local data chunk to storage the field points values for the differen degrees of freedom and headings
        this->_field_data.resize( this->_size_local_fp );
        this->_is_heap          = true;

        // Generate temporary array to store field points coordinates according to the current process distribution along field points
        cusfloat*   field_points_l = generate_empty_vector<cusfloat>( this->_size_local * 3 );
        std::size_t start_pos      = mpi_config_->proc_rank * this->_size_local; 
        for ( std::size_t i=0; i<this->_size_local; i++)
        {
            std::size_t gindex          = start_pos + i;
            field_points_l[ 3*i + 0 ]   = mesh_->x[ gindex ];
            field_points_l[ 3*i + 1 ]   = mesh_->y[ gindex ];
            field_points_l[ 3*i + 2 ]   = mesh_->z[ gindex ];
        }

        // Allocate PanelData
        this->panel_data.reserve( this->_size_local );

        std::size_t field_points_np = 1; // Dummy value since field points will be loaded from mesh nodes
        for ( std::size_t i=0; i<this->_size_local; i++ )
        {
            this->panel_data.emplace_back( 
                                                PanelData<T, Config>(
                                                                        field_points_np,
                                                                        &(field_points_l[ 3*i ]),
                                                                        freqs_np_,
                                                                        headings_np_,
                                                                        dofs_np_
                                                                    ) 
                                        );
        }

        // Delete temporary heap array
        mkl_free( field_points_l );
    }
    else
    {
        // Store number of field points
        this->_size_global      = mesh_->elems_np;
        this->_mpi_config       = mpi_config_;

        // Calculate start and end positions for the current process
        this->_mpi_config->get_1d_bounds(
                                            static_cast<int>( this->_size_global ),
                                            reinterpret_cast<int&>( this->_start_pos ),
                                            reinterpret_cast<int&>( this->_end_pos )
                                        );

        // Determine size of the local panel set for the current process
        this->_size_local       = this->_end_pos - this->_start_pos;

        // Allocate PanelData
        this->panel_data.reserve( this->_size_local );

        std::size_t body_id         = 0;
        std::size_t global_panel_id = 0;
        for ( std::size_t i=0; i<this->_size_local; i++ )
        {
            // Calculate global panel ID
            global_panel_id = this->_start_pos + i;

            // Create new panel data
            this->panel_data.emplace_back( 
                                                PanelData<T, Config>(
                                                                        mesh_->panels[ global_panel_id ],
                                                                        body_id,
                                                                        freqs_np_,
                                                                        headings_np_,
                                                                        dofs_np_,
                                                                        use_waterline_
                                                                    ) 
                                        );
        }

        // Determine size of the local data chunk to storage the field points values for the differen degrees of freedom and headings
        this->_size_local_fp = 0;
        for ( std::size_t i=0; i<this->_size_local; i++ )
        {
            this->_size_local_fp += this->panel_data[i].field_points_np;
        }

        // Allocate local data chunk to storage the field points values for the differen degrees of freedom and headings
        this->_field_data.resize( this->_size_local_fp );
        this->_is_heap    = true;
    }
}


template<typename T, typename Config>
RadDiffData<T, Config>::RadDiffData( 
                                    MpiConfig*      mpi_config_,
                                    MeshGroup*      mesh_gp_,
                                    std::size_t     freqs_np_,
                                    std::size_t     headings_np_,
                                    std::size_t     dofs_np_,
                                    bool            use_waterline_
                                )
{    
    // Get panel set view
    PanelSetView panel_view = make_panel_view( mesh_gp_, use_waterline_ );
    
    // Store number of field points
    this->_size_global      = panel_view.panels_tnp;
    this->_mpi_config       = mpi_config_;

    // Calculate start and end positions for the current process
    this->_mpi_config->get_1d_bounds(
                                        static_cast<int>( this->_size_global ),
                                        reinterpret_cast<int&>( this->_start_pos ),
                                        reinterpret_cast<int&>( this->_end_pos )
                                    );

    // Determine size of the local data chunk to storage the field points values for the differen degrees of freedom and headings
    this->_size_local       = this->_end_pos - this->_start_pos;

    // Allocate PanelData
    this->panel_data.reserve( this->_size_local );

    std::size_t body_id         = 0;
    std::size_t global_panel_id = 0;
    for ( std::size_t i=0; i<this->_size_local; i++ )
    {
        // Calculate global panel ID
        global_panel_id = this->_start_pos + i;

        // Search current body ID
        for ( std::size_t j=body_id; j<static_cast<std::size_t>( mesh_gp_->meshes_np ); j++ )
        {
            if ( global_panel_id < static_cast<std::size_t>( panel_view.panels_cnp[j+1] ) )
            {
                body_id     = j;
                break;
            }
        }

        // Create new panel data
        this->panel_data.emplace_back( 
                                            PanelData<T, Config>(
                                                                    panel_view.panels[ global_panel_id ],
                                                                    body_id,
                                                                    freqs_np_,
                                                                    headings_np_,
                                                                    dofs_np_,
                                                                    use_waterline_
                                                                ) 
                                    );
    }

    // Determine size of the local data chunk to storage the field points values for the differen degrees of freedom and headings
    this->_size_local_fp = 0;
    for ( std::size_t i=0; i<this->_size_local; i++ )
    {
        this->_size_local_fp += this->panel_data[i].field_points_np;
    }

    // Allocate local data chunk to storage the field points values for the differen degrees of freedom and headings
    this->_field_data.resize( this->_size_local_fp );
    this->_is_heap    = true;
}


template<typename T, typename Config>
RadDiffData<T, Config>::~RadDiffData( )
{
}