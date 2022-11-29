
// Include external libraries
#include <iostream>

// Include external scientific libraries
#include <cmath>
#include "mkl.h"

// Include local modules
#include "common.hpp"
#include "../containers.hpp"
#include "dipole.hpp"
#include "../math_interface.hpp"
#include "../math_tools.hpp"


void calculate_dipole_potential(PanelGeom &panel, cusfloat (&field_point)[3], int fp_local_flag, cusfloat &phi)
{
    // Translate field point to panel local coordinates
    cusfloat field_point_local[3];
    if (fp_local_flag == 0)
    {
        cusfloat field_point_local_aux[3];
        sv_sub(3, field_point, panel.center, field_point_local_aux);
        cblas_gemv<cusfloat>(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, panel.global_to_local_mat, 3, field_point_local_aux, 1, 0, field_point_local, 1);
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
    cusfloat node_fieldp_dx[panel.MAX_PANEL_NODES];
    cusfloat node_fieldp_dy[panel.MAX_PANEL_NODES];
    cusfloat node_fieldp_dz[panel.MAX_PANEL_NODES];
    cusfloat node_fieldp_mod[panel.MAX_PANEL_NODES];
    calculate_distance_node_field(panel, field_point_local, node_fieldp_mod, node_fieldp_dx, node_fieldp_dy, node_fieldp_dz);

    // Calculate distances in between nodes
    cusfloat delta_xi [panel.MAX_PANEL_NODES];
    cusfloat delta_eta [panel.MAX_PANEL_NODES];
    calculate_nodes_distance(panel, delta_xi, delta_eta);

    // Calculate potential kernel
    calculate_dipole_potential_kernel(panel, node_fieldp_mod, node_fieldp_dx, node_fieldp_dy, node_fieldp_dz, delta_xi, delta_eta, phi);
   
}


void calculate_dipole_potential_kernel(PanelGeom &panel, cusfloat * node_fieldp_mod, cusfloat* node_fieldp_dx, 
    cusfloat* node_fieldp_dy, cusfloat* node_fieldp_dz, cusfloat* delta_xi, cusfloat* delta_eta, cusfloat &phi)
{
    cusfloat a0i, a0j, a1i, a1j;
    int j = 0;
    for (int i=0; i<panel.num_nodes-1; i++)
    {
        // Calculate secondary index
        j = (i+1)%panel.num_nodes;

        // Calcualte potential coefficients
        a0i = delta_eta[i]*(pow2s(node_fieldp_dx[i])+pow2s(node_fieldp_dz[i]));
        a0i -= delta_xi[i]*node_fieldp_dx[i]*node_fieldp_dy[i];
        a1i = node_fieldp_mod[i]*node_fieldp_dz[i]*delta_xi[i];

        a0j = delta_eta[j]*(pow2s(node_fieldp_dx[j])+pow2s(node_fieldp_dz[j]));
        a0j -= delta_xi[j]*node_fieldp_dx[j]*node_fieldp_dy[j];
        a1j = node_fieldp_mod[j]*node_fieldp_dz[j]*delta_xi[j];

        phi += atan(a0i/a1i) - atan(a0j/a1j);

    }

}


void calculate_dipole_velocity(PanelGeom &panel, cusfloat (&field_point)[3], int fp_local_flag, cusfloat (&velocity)[3])
{
    // Translate field point to panel local coordinates
    cusfloat field_point_local[3];
    if (fp_local_flag == 0)
    {
        cusfloat field_point_local_aux[3];
        sv_sub(3, field_point, panel.center, field_point_local_aux);
        cblas_gemv<cusfloat>(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, panel.global_to_local_mat, 3, field_point_local_aux, 1, 0, field_point_local, 1);
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
    cusfloat node_fieldp_dx[panel.MAX_PANEL_NODES];
    cusfloat node_fieldp_dy[panel.MAX_PANEL_NODES];
    cusfloat node_fieldp_dz[panel.MAX_PANEL_NODES];
    cusfloat node_fieldp_mod[panel.MAX_PANEL_NODES];
    calculate_distance_node_field(panel, field_point_local, node_fieldp_mod, node_fieldp_dx, node_fieldp_dy, node_fieldp_dz);

    // Calculate distances in between nodes
    cusfloat delta_xi [panel.MAX_PANEL_NODES];
    cusfloat delta_eta [panel.MAX_PANEL_NODES];
    calculate_nodes_distance(panel, delta_xi, delta_eta);

    // Calculate velocities
    calculate_dipole_velocity_kernel(panel, node_fieldp_mod, node_fieldp_dx, node_fieldp_dy, node_fieldp_dz, delta_xi, delta_eta, velocity);

}


void calculate_dipole_velocity_kernel(PanelGeom &panel, cusfloat * node_fieldp_mod, cusfloat* node_fieldp_dx, 
    cusfloat* node_fieldp_dy, cusfloat* node_fieldp_dz, cusfloat* delta_xi, cusfloat* delta_eta, cusfloat* velocity,
    int potential_flag)
{
    cusfloat a0i, a0j, a1i, a1j;
    cusfloat da0i, da0j, da1i, da1j;
    cusfloat fi, fj;
    int j = 0;
    for (int i=0; i<panel.num_nodes; i++)
    {
        // Calculate secondary index
        j = (i+1)%panel.num_nodes;

        // Calculate argument for the derivative function
        a0i = delta_eta[i]*(pow2s(node_fieldp_dx[i])+pow2s(node_fieldp_dz[i]));
        a0i -= delta_xi[i]*node_fieldp_dx[i]*node_fieldp_dy[i];
        a1i = node_fieldp_mod[i]*node_fieldp_dz[i]*delta_xi[i];

        a0j = delta_eta[j]*(pow2s(node_fieldp_dx[j])+pow2s(node_fieldp_dz[j]));
        a0j -= delta_xi[j]*node_fieldp_dx[j]*node_fieldp_dy[j];
        a1j = node_fieldp_mod[j]*node_fieldp_dz[j]*delta_xi[j];

        if(potential_flag == 1)
        {
            velocity[4] += atan(a0i/a1i) - atan(a0j/a1j);
        }

        fi = 1/(1+pow2s(a0i/a1i));
        fj = 1/(1+pow2s(a0j/a1j));

        // Calculate X derivative of the argument
        da0i = 2*node_fieldp_dx[i]*delta_eta[i]-delta_xi[i]*node_fieldp_dy[i];
        da1i = node_fieldp_dx[i]*node_fieldp_dz[i]*delta_xi[i]/node_fieldp_mod[i];

        da0j = 2*node_fieldp_dx[j]*delta_eta[j]-delta_xi[j]*node_fieldp_dy[j];
        da1j= node_fieldp_dx[j]*node_fieldp_dz[j]*delta_xi[j]/node_fieldp_mod[j];

        velocity[0] += (da0i*a1i-a0i*da1i)/pow2s(a1i)*fi - (da0j*a1j-a0j*da1j)/pow2s(a1j)*fj;

        // Calculate Y derivative of the argument
        da0i = -delta_xi[i]*node_fieldp_dx[i];
        da1i = node_fieldp_dy[i]*node_fieldp_dz[i]*delta_xi[i]/node_fieldp_mod[i];

        da0j = -delta_xi[j]*node_fieldp_dx[j];
        da1j = node_fieldp_dy[j]*node_fieldp_dz[j]*delta_xi[j]/node_fieldp_mod[j];

        velocity[1] += (da0i*a1i-a0i*da1i)/pow2s(a1i)*fi - (da0j*a1j-a0j*da1j)/pow2s(a1j)*fj;

        // Calculate Z derivative of the argument
        da0i = 2*node_fieldp_dz[i]*delta_eta[i];
        da1i = pow2s(node_fieldp_dz[i])*delta_xi[i]/node_fieldp_mod[i] + node_fieldp_mod[i]*delta_xi[i];

        da0j = 2*node_fieldp_dz[j]*delta_eta[j];
        da1j = pow2s(node_fieldp_dz[j])*delta_xi[j]/node_fieldp_mod[j] + node_fieldp_mod[j]*delta_xi[j];

        velocity[2] += (da0i*a1i-a0i*da1i)/pow2s(a1i)*fi - (da0j*a1j-a0j*da1j)/pow2s(a1j)*fj;

    }
}