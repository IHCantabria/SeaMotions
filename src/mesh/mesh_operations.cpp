
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
#include "../math/euler_transforms.hpp"
#include "mesh_operations.hpp"


void    mesh_points_rotation(
                                    int         nodes_np,
                                    cusfloat*   x,
                                    cusfloat*   y,
                                    cusfloat*   z,
                                    cusfloat    rx,
                                    cusfloat    ry,
                                    cusfloat    rz
                            )
{
    // Calculate rotation matrix
    cusfloat rot_mat[9];
    euler_ypr( rx, ry, rz, rot_mat );

    // Loop over points to rotate them
    cusfloat xa, ya, za;
    for ( int i=0; i<nodes_np; i++ )
    {
        // Copy to aux variables
        xa = x[i];
        ya = y[i];
        za = z[i];

        // Rotate points using rotation matrix
        x[i] = rot_mat[0] * xa + rot_mat[1] * ya + rot_mat[2] * za;
        x[i] = rot_mat[3] * xa + rot_mat[4] * ya + rot_mat[5] * za;
        x[i] = rot_mat[6] * xa + rot_mat[7] * ya + rot_mat[8] * za;

    }

}


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
                                )
{
    // Move the mesh to have the reference of rotation on
    // P = ( 0, 0, 0 ) so rotations system is working properly
    mesh_points_translation(
                                nodes_np,
                                x,
                                y,
                                z,
                                -sref_x,
                                -sref_y,
                                -sref_z
                            );

    // Rotate mesh around the refernce point. At this 
    // point at P = ( 0, 0, 0 )
    mesh_points_rotation(
                                nodes_np,
                                x,
                                y,
                                z,
                                rx,
                                ry,
                                rz
                        );
    
    // Move back the mesh to have the reference point in 
    // its corresponding position w.r.t to the global 
    // reference system
    mesh_points_translation(
                                nodes_np,
                                x,
                                y,
                                z,
                                sref_x,
                                sref_y,
                                sref_z
                            );

}


void    mesh_points_translation(
                                    int         nodes_np,
                                    cusfloat*   x,
                                    cusfloat*   y,
                                    cusfloat*   z,
                                    cusfloat    sref_x,
                                    cusfloat    sref_y,
                                    cusfloat    sref_z
                                )
{
    for ( int i=0; i<nodes_np; i++ )
    {
        x[i] += sref_x;
        y[i] += sref_y;
        z[i] += sref_z;
    }
}