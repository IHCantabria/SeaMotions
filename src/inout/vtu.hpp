
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
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

// Include local modules
#include "../config.hpp"


bool    write_vtu_binary_appended(
                                    const std::string       &filename,
                                    const std::size_t       nodes_np,
                                    cusfloat*               nodes_x,
                                    cusfloat*               nodes_y,
                                    cusfloat*               nodes_z,
                                    const std::size_t       elems_np,
                                    const int               enrl,
                                    int*                    elements,
                                    int*                    elements_type
                                );


bool    write_vtu_ascii(
                                    const std::string       &filename,
                                    const std::size_t       nodes_np,
                                    cusfloat*               nodes_x,
                                    cusfloat*               nodes_y,
                                    cusfloat*               nodes_z,
                                    const std::size_t       elems_np,
                                    const int               enrl,
                                    int*                    elements,
                                    int*                    elements_type
                        );


/**
 * @brief Write a VTU (UnstructuredGrid) file with per-panel pressure cell data.
 *
 * Each panel is represented as an independent polygon cell (no node sharing).
 * The pressure scalar field is stored as CellData.
 * Format: binary-appended, little-endian, block lengths as UInt32.
 *
 * @param filename      Full path of the output VTU file.
 * @param n_nodes       Total number of corner nodes (sum of num_nodes over all panels).
 * @param nodes_x       Flat array of node x-coordinates (length n_nodes).
 * @param nodes_y       Flat array of node y-coordinates (length n_nodes).
 * @param nodes_z       Flat array of node z-coordinates (length n_nodes).
 * @param n_cells       Number of panels / cells.
 * @param connectivity  Node indices per cell, concatenated (length = n_nodes).
 * @param offsets       Cumulative node count per cell (length = n_cells).
 * @param types         VTK cell-type code per cell (5=TRI, 9=QUAD, 7=POLYGON).
 * @param pressure      Pressure value per cell (length = n_cells).
 * @return true on success.
 */
bool    write_vtu_panel_pressure(
                                    const std::string&      filename,
                                    std::size_t             n_nodes,
                                    const cusfloat*         nodes_x,
                                    const cusfloat*         nodes_y,
                                    const cusfloat*         nodes_z,
                                    std::size_t             n_cells,
                                    const int32_t*          connectivity,
                                    const int32_t*          offsets,
                                    const uint8_t*          types,
                                    const cusfloat*         pressure
                                );


/**
 * @brief Write a PVD (ParaView Data) collection file linking VTU time steps.
 *
 * @param filename   Full path of the output PVD file.
 * @param timesteps  Ordered list of (time [s], relative VTU file path) pairs.
 * @return true on success.
 */
bool    write_pvd(
                    const std::string&                                  filename,
                    const std::vector<std::pair<double, std::string>>&  timesteps
                );