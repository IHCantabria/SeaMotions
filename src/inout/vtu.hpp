
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