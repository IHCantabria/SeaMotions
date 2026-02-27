
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

// Include general usage libraries
#include <string>
#include <vector>

// Include local modules
#include "../config.hpp"


struct FieldMeshDataConfig
{
    int                     compression_level       = 0;        // Compression level for output files (0-9, where 0 is no compression and 9 is maximum compression). Not used if MPI is enabled.
    std::size_t             dofs_np                 = 0;        // Number of degrees of freedom for the field data (e.g., 1 for scalar fields, 3 for vector fields, etc.)
    std::size_t             freqs_np                = 0;        // Number of frequencies for the field data
    std::vector<cusfloat*>  headings                = {};       // Array of heading values for the field data
    std::size_t             headings_np             = 0;        // Number of headings for the field data
    std::string             mesh_file_path          ;           // Path to the mesh file containing the field points
    std::string             body_name               ;           // Name of the body to extract from the mesh file
    bool                    out_components          = false;    // Flag to indicate if field components data should be outputted
    std::string             out_folder_path         ;           // Path to the output folder where field data files will be stored
    bool                    out_parallel            = false;    // Flag to indicate if output files should be written in parallel using MPI. If true, compression_level is ignored and no compression is applied.
    bool                    out_potential           = false;    // Flag to indicate if potential field data should be outputted
    bool                    out_pressure            = false;    // Flag to indicate if pressure field data should be outputted
    bool                    out_velocity            = false;    // Flag to indicate if velocity field data should be outputted
};