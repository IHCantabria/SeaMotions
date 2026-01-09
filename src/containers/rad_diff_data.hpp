
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

#pragma once

// Include local modules
#include "../math/custensor/custensor.hpp"
#include "../mesh/mesh_group.hpp"
#include "mpi_config.hpp"
#include "panel_data.hpp"


template<int mode_f, int mode_dfdn, int mode_dfdc>
struct RadDiffData
{
private:
    // Declare private variables
    std::size_t                     _end_pos        = 0;        // End position along field points for the current process
    bool                            _is_heap        = false;    // Flag to indicate if memory is allocated on heap
    MpiConfig*                      _mpi_config     = nullptr;  // Pointer to MPI configuration
    std::size_t                     _size_global    = 0;        // Total number of field points all across processes
    std::size_t                     _size_local     = 0;        // Total number of field points for the current process
    std::size_t                     _start_pos      = 0;        // Start position along field points for the current process
    
public:
    // Declare public variables
    std::vector<PanelData<mode_f, mode_dfdn, mode_dfdc>>   panel_data; // Store panel data for radiation and diffraction calculations

    /* Declare class constructors and destructor */
    RadDiffData( ) = default;

    RadDiffData( 
                    MpiConfig*      mpi_config_,
                    std::size_t     panels_np_,
                    std::size_t     field_points_np_,
                    std::size_t     headings_np_,
                    std::size_t     dofs_np_
                );

    RadDiffData( 
                    MpiConfig*      mpi_config_,
                    MeshGroup*      mesh_gp_,
                    std::size_t     headings_np_,
                    std::size_t     dofs_np_,
                    bool            use_waterline_ = false
                );

    /* Declare public methods */
    std::size_t get_end_pos( void ) const;

    std::size_t get_size_global( void ) const;

    std::size_t get_size_local( void ) const;

    std::size_t get_start_pos( void ) const;
    
};


// Include template implementation file
#include "rad_diff_data.txx"