
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


/**
 * @brief   Function to rotate mesh points w.r.t the origin P = ( 0, 0, 0 )
 * 
* @param    Number of nodes
 * @param   Nodes X coordinates
 * @param   Nodes Y coordinates
 * @param   Nodes Z coordinates
 * @param   Angle of rotation w.r.t X axis going through the origin
 * @param   Angle of rotation w.r.t Y axis going through the origin
 * @param   Angle of rotation w.r.t Z axis going through the origin
 */
void    mesh_points_rotation(
                                    int         nodes_np,
                                    cusfloat*   x,
                                    cusfloat*   y,
                                    cusfloat*   z,
                                    cusfloat    rx,
                                    cusfloat    ry,
                                    cusfloat    rz
                            );


/**
 * @brief   Function to rotate mesh points w.r.t the origin P = ( 0, 0, 0 )
 * 
* @param    Number of nodes
 * @param   Nodes X coordinates
 * @param   Nodes Y coordinates
 * @param   Nodes Z coordinates
 * @param   X reference point where to move the mesh
 * @param   Y reference point where to move the mesh
 * @param   Z reference point where to move the mesh
 * @param   Angle of rotation w.r.t X axis going through the origin
 * @param   Angle of rotation w.r.t Y axis going through the origin
 * @param   Angle of rotation w.r.t Z axis going through the origin
 */
void    mesh_points_rotation_refp(
                                    int         nodes_np,
                                    cusfloat*   x,
                                    cusfloat*   y,
                                    cusfloat*   z,
                                    cusfloat    sref_x,
                                    cusfloat    sref_y,
                                    cusfloat    sref_z,
                                    cusfloat    rx,
                                    cusfloat    ry,
                                    cusfloat    rz
                                );


/**
 * @brief   Function to translate mesh points to a reference point
 * 
 * @param   Number of nodes
 * @param   Nodes X coordinates
 * @param   Nodes Y coordinates
 * @param   Nodes Z coordinates
 * @param   X reference point where to move the mesh
 * @param   Y reference point where to move the mesh
 * @param   Z reference point where to move the mesh
 */
void    mesh_points_translation(
                                    int         nodes_np,
                                    cusfloat*   x,
                                    cusfloat*   y,
                                    cusfloat*   z,
                                    cusfloat    sref_x,
                                    cusfloat    sref_y,
                                    cusfloat    sref_z
                                );