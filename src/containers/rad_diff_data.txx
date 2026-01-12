
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


template<typename Config>
std::size_t RadDiffData<Config>::get_end_pos( void ) const
{
    return this->_end_pos;
}


template<typename Config>
std::size_t RadDiffData<Config>::get_size_global( void ) const
{
    return this->_size_global;
}


template<typename Config>
std::size_t RadDiffData<Config>::get_size_local( void ) const
{
    return this->_size_local;
}


template<typename Config>
std::size_t RadDiffData<Config>::get_start_pos( void ) const
{
    return this->_start_pos;
}


template<typename Config>
RadDiffData<Config>::RadDiffData( 
                                    MpiConfig*      mpi_config_,
                                    std::size_t     panels_np_,
                                    std::size_t     field_points_np_,
                                    std::size_t     freqs_np_,
                                    std::size_t     headings_np_,
                                    std::size_t     dofs_np_
                                )
{
    // Store number of field points
    this->_size_global      = panels_np_;
    this->_mpi_config       = mpi_config_;

    // Calculate start and end positions for the current process
    this->_mpi_config->get_1d_bounds(
                                        static_cast<int>( panels_np_ ),
                                        reinterpret_cast<int&>( this->_start_pos ),
                                        reinterpret_cast<int&>( this->_end_pos )
                                    );

    this->_size_local       = this->_end_pos - this->_start_pos;

    // Allocate PanelData
    this->panel_data.reserve( this->_size_local );

    for ( std::size_t i=0; i<this->_size_local; i++ )
    {
        this->panel_data.emplace_back( 
                                            PanelData<mode_comp, mode_f, mode_dfdn, mode_dfdc>(
                                                                                                    field_points_np_,
                                                                                                    freqs_np_,
                                                                                                    headings_np_,
                                                                                                    dofs_np_
                                                                                                ) 
                                    );
    }

}


template<typename Config>
RadDiffData<Config>::RadDiffData( 
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
                                            PanelData<Config>(
                                                                    panel_view.panels[ global_panel_id ],
                                                                    body_id,
                                                                    freqs_np_,
                                                                    headings_np_,
                                                                    dofs_np_,
                                                                    use_waterline_
                                                                ) 
                                    );
    }
}