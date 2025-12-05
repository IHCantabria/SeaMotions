
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


// Include external libraries
#include <iostream>

// Include external scientific libraries
#include <cmath>
#include "mkl.h"

// Include local modules
#include "common.hpp"
#include "../containers/containers.hpp"
#include "dipole.hpp"
#include "../math/math_interface.hpp"
#include "../math/math_tools.hpp"


void    calculate_dipole_potential(
                                            PanelGeom*      panel, 
                                            cusfloat*       field_point, 
                                            int             fp_local_flag, 
                                            cusfloat&       phi
                                )
{
    // Translate field point to panel local coordinates
    cusfloat field_point_local[3];
    if (fp_local_flag == 0)
    {
        cusfloat field_point_local_aux[3];
        sv_sub(3, field_point, panel->center, field_point_local_aux);
        cblas_gemv<cusfloat>(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, panel->global_to_local_mat, 3, field_point_local_aux, 1, 0, field_point_local, 1);
    }
    else if (fp_local_flag == 1)
    {
        field_point_local[0] = field_point[0];
        field_point_local[1] = field_point[1];
        field_point_local[2] = field_point[2];
    }
    else
    {
        throw std::runtime_error("Local flag must have values of: 0 (non-local) and 1 (local)");
    }

    // Calculate distances from each node to the field point in local coordinates
    cusfloat node_fieldp_dx[panel->MAX_PANEL_NODES];
    cusfloat node_fieldp_dy[panel->MAX_PANEL_NODES];
    cusfloat node_fieldp_dz[panel->MAX_PANEL_NODES];
    cusfloat node_fieldp_mod[panel->MAX_PANEL_NODES];
    calculate_distance_node_field(panel, field_point_local, node_fieldp_mod, node_fieldp_dx, node_fieldp_dy, node_fieldp_dz);

    // Calculate distances in between nodes
    cusfloat delta_xi [panel->MAX_PANEL_NODES];
    cusfloat delta_eta [panel->MAX_PANEL_NODES];
    calculate_nodes_distance(panel, delta_xi, delta_eta);

    // Calculate potential kernel
    calculate_dipole_potential_kernel(panel, node_fieldp_mod, node_fieldp_dx, node_fieldp_dy, node_fieldp_dz, delta_xi, delta_eta, phi);
   
}


void    calculate_dipole_potential_kernel(
                                            PanelGeom*      panel, 
                                            cusfloat*       node_fieldp_mod,
                                            cusfloat*       node_fieldp_dx,
                                            cusfloat*       node_fieldp_dy, 
                                            cusfloat*       node_fieldp_dz, 
                                            cusfloat*       delta_xi, 
                                            cusfloat*       delta_eta, 
                                            cusfloat&       phi
                                        )
{
    cusfloat a0i, a0i2, a0j, a0j2, a1i, a1j;
    int j = 0;
    for (int i=0; i<panel->num_nodes; i++)
    {
        // Calculate secondary index
        j = (i+1)%panel->num_nodes;

        // Calcualte potential coefficients
        a0i = delta_eta[i]*(pow2s(node_fieldp_dx[i])+pow2s(node_fieldp_dz[i]));
        a0i2 = delta_xi[i]*node_fieldp_dx[i]*node_fieldp_dy[i];
        a1i = node_fieldp_mod[i]*node_fieldp_dz[i]*delta_xi[i];

        a0j = delta_eta[i]*(pow2s(node_fieldp_dx[j])+pow2s(node_fieldp_dz[j]));
        a0j2 = delta_xi[i]*node_fieldp_dx[j]*node_fieldp_dy[j];
        a1j = node_fieldp_mod[j]*node_fieldp_dz[j]*delta_xi[i];

        phi += atan((a0i-a0i2)/a1i) - atan((a0j-a0j2)/a1j);

    }

}


void    calculate_dipole_velocity(
                                            PanelGeom*      panel, 
                                            cusfloat*       field_point,
                                            int             fp_local_flag, 
                                            cusfloat        *velocity
                                )
{
    // Translate field point to panel local coordinates
    cusfloat field_point_local[3];
    if (fp_local_flag == 0)
    {
        cusfloat field_point_local_aux[3];
        sv_sub(3, field_point, panel->center, field_point_local_aux);
        cblas_gemv<cusfloat>(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, panel->global_to_local_mat, 3, field_point_local_aux, 1, 0, field_point_local, 1);
    }
    else if (fp_local_flag == 1)
    {
        field_point_local[0] = field_point[0];
        field_point_local[1] = field_point[1];
        field_point_local[2] = field_point[2];
    }
    else
    {
        throw std::runtime_error("Local flag must have values of: 0 (non-local) and 1 (local)");
    }

    // Calculate distances from each node to the field point in local coordinates
    cusfloat node_fieldp_dx[panel->MAX_PANEL_NODES];
    cusfloat node_fieldp_dy[panel->MAX_PANEL_NODES];
    cusfloat node_fieldp_dz[panel->MAX_PANEL_NODES];
    cusfloat node_fieldp_mod[panel->MAX_PANEL_NODES];
    calculate_distance_node_field(panel, field_point_local, node_fieldp_mod, node_fieldp_dx, node_fieldp_dy, node_fieldp_dz);

    // Calculate distances in between nodes
    cusfloat delta_xi [panel->MAX_PANEL_NODES];
    cusfloat delta_eta [panel->MAX_PANEL_NODES];
    calculate_nodes_distance(panel, delta_xi, delta_eta);

    // Calculate velocities
    calculate_dipole_velocity_kernel(panel, node_fieldp_mod, node_fieldp_dx, node_fieldp_dy, node_fieldp_dz, delta_xi, delta_eta, velocity);

}


void    calculate_dipole_velocity_kernel(
                                            PanelGeom*      panel, 
                                            cusfloat*       node_fieldp_mod,
                                            cusfloat*       node_fieldp_dx, 
                                            cusfloat*       node_fieldp_dy, 
                                            cusfloat*       node_fieldp_dz,
                                            cusfloat*       delta_xi,
                                            cusfloat*       delta_eta, 
                                            cusfloat*   velocity
                                        )
{
    cusfloat a0i, a0j, a1i, a1j;
    cusfloat da0i, da0j, da1i, da1j;
    cusfloat fi, fj;
    int j = 0;
    for (int i=0; i<panel->num_nodes; i++)
    {
        // Calculate secondary index
        j = (i+1)%panel->num_nodes;

        // Calculate argument for the derivative function
        a0i = delta_eta[i]*(pow2s(node_fieldp_dx[i])+pow2s(node_fieldp_dz[i]));
        a0i -= delta_xi[i]*node_fieldp_dx[i]*node_fieldp_dy[i];
        a1i = node_fieldp_mod[i]*node_fieldp_dz[i]*delta_xi[i];

        a0j = delta_eta[i]*(pow2s(node_fieldp_dx[j])+pow2s(node_fieldp_dz[j]));
        a0j -= delta_xi[i]*node_fieldp_dx[j]*node_fieldp_dy[j];
        a1j = node_fieldp_mod[j]*node_fieldp_dz[j]*delta_xi[i];

        fi = 1/(1+pow2s(a0i/a1i));
        fj = 1/(1+pow2s(a0j/a1j));

        // Calculate X derivative of the argument
        da0i = 2.0*delta_eta[i]*node_fieldp_dx[i] - delta_xi[i]*node_fieldp_dy[i];
        da0j = 2.0*delta_eta[i]*node_fieldp_dx[j] - delta_xi[i]*node_fieldp_dy[j];

        velocity[0] += da0i*a1i/pow2s(a1i)*fi - da0j*a1j/pow2s(a1j)*fj;

        // Calculate Y derivative of the argument
        da0i = -delta_xi[i]*node_fieldp_dx[i];
        da0j = -delta_xi[i]*node_fieldp_dx[j];

        velocity[1] += da0i*a1i/pow2s(a1i)*fi - da0j*a1j/pow2s(a1j)*fj;

        // Calculate Z derivative of the argument
        da0i = 2*node_fieldp_dz[i]*delta_eta[i];
        da1i = node_fieldp_mod[i]*delta_xi[i];

        da0j = 2*node_fieldp_dz[j]*delta_eta[i];
        da1j = node_fieldp_mod[j]*delta_xi[i];

        velocity[2] += (da0i*a1i-a0i*da1i)/pow2s(a1i)*fi - (da0j*a1j-a0j*da1j)/pow2s(a1j)*fj;

    }
}