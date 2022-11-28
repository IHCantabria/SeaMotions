
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
    cusfloat a0, a1, a2;
    cusfloat b0, b1, b2;
    for (int i=0; i<panel.num_nodes-1; i++)
    {
        a0 = delta_eta[i]*(pow2s(node_fieldp_dx[i])+pow2s(node_fieldp_dz[i]));
        a1 = delta_xi[i]*node_fieldp_dx[i]*node_fieldp_dy[i];
        a2 = node_fieldp_mod[i]*node_fieldp_dz[i]*delta_xi[i];
        b0 = delta_eta[i]*(pow2s(node_fieldp_dx[i+1])+pow2s(node_fieldp_dz[i+1]));
        b1 = delta_xi[i]*node_fieldp_dx[i+1]*node_fieldp_dy[i+1];
        b2 = node_fieldp_mod[i+1]*node_fieldp_dz[i+1]*delta_xi[i];
        phi += atan((a0-a1)/a2) - atan((b0-b1)/b2);
    }
    int n = panel.num_nodes-1;
    a0 = delta_eta[n]*(pow2s(node_fieldp_dx[n])+pow2s(node_fieldp_dz[n]));
    a1 = delta_xi[n]*node_fieldp_dx[n]*node_fieldp_dy[n];
    a2 = node_fieldp_mod[n]*node_fieldp_dz[n]*delta_xi[n];
    b0 = delta_eta[n]*(pow2s(node_fieldp_dx[0])+pow2s(node_fieldp_dz[0]));
    b1 = delta_xi[n]*node_fieldp_dx[0]*node_fieldp_dy[0];
    b2 = node_fieldp_mod[0]*node_fieldp_dz[0]*delta_xi[n];
    phi += atan((a0-a1)/a2) - atan((b0-b1)/b2);
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
    cusfloat* node_fieldp_dy, cusfloat* node_fieldp_dz, cusfloat* delta_xi, cusfloat* delta_eta, cusfloat (&velocity)[3])
{
    cusfloat a0, a1, a2;
    cusfloat b0, b1, b2;
    cusfloat der_arg0, der_arg1;
    cusfloat f0, f1;
    for (int i=0; i<panel.num_nodes-1; i++)
    {
        // Calculate argument for the derivative function
        a0 = delta_eta[i]*(pow2s(node_fieldp_dx[i])+pow2s(node_fieldp_dz[i]));
        a1 = delta_xi[i]*node_fieldp_dx[i]*node_fieldp_dy[i];
        a2 = node_fieldp_mod[i]*node_fieldp_dz[i]*delta_xi[i];
        b0 = delta_eta[i]*(pow2s(node_fieldp_dx[i+1])+pow2s(node_fieldp_dz[i+1]));
        b1 = delta_xi[i]*node_fieldp_dx[i+1]*node_fieldp_dy[i+1];
        b2 = node_fieldp_mod[i+1]*node_fieldp_dz[i+1]*delta_xi[i];
        f0 = 1/(1+pow2s((a0-a1)/a2));
        f1 = 1/(1+pow2s((b0-b1)/b2));

        // Calculate X derivative of the argument
        der_arg0 = (2*delta_eta[i]*node_fieldp_dx[i]-delta_xi[i]*node_fieldp_dy[i])/(node_fieldp_mod[i]*node_fieldp_dz[i]*delta_xi[i]);
        der_arg1 = (2*delta_eta[i+1]*node_fieldp_dx[i+1]-delta_xi[i+1]*node_fieldp_dy[i+1])/(node_fieldp_mod[i+1]*node_fieldp_dz[i+1]*delta_xi[i+1]);
        velocity[0] += der_arg0*f0 - der_arg1*f1;

        // Calculate Y derivative of the argument
        der_arg0 = -delta_xi[i]*node_fieldp_dx[i]/(node_fieldp_mod[i]*node_fieldp_dz[i]*delta_xi[i]);
        der_arg1 = -delta_xi[i+1]*node_fieldp_dx[i+1]/(node_fieldp_mod[i+1]*node_fieldp_dz[i+1]*delta_xi[i+1]);
        velocity[1] += der_arg0*f0 - der_arg1*f1;

        // Calculate Z derivative of the argument
        der_arg0 = 2*pow2s(node_fieldp_dz[i])*delta_eta[i]*node_fieldp_mod[i]*delta_xi[i];
        der_arg0 -= (a0-a1)*node_fieldp_mod[i]*delta_xi[i];
        der_arg0 /= pow2s(node_fieldp_mod[i]*node_fieldp_dz[i]*delta_xi[i]);

        der_arg1 = 2*pow2s(node_fieldp_dz[i+1])*delta_eta[i+1]*node_fieldp_mod[i+1]*delta_xi[i+1];
        der_arg1 -= (b0-b1)*node_fieldp_mod[i+1]*delta_xi[i+1];
        der_arg1 /= pow2s(node_fieldp_mod[i+1]*node_fieldp_dz[i+1]*delta_xi[i+1]);

        velocity[2] += der_arg0*f0 - der_arg1*f1;

    }
    int n = panel.num_nodes-1;
    a0 = delta_eta[n]*(pow2s(node_fieldp_dx[n])+pow2s(node_fieldp_dz[n]));
    a1 = delta_xi[n]*node_fieldp_dx[n]*node_fieldp_dy[n];
    a2 = node_fieldp_mod[n]*node_fieldp_dz[n]*delta_xi[n];
    b0 = delta_eta[0]*(pow2s(node_fieldp_dx[0])+pow2s(node_fieldp_dz[0]));
    b1 = delta_xi[0]*node_fieldp_dx[0]*node_fieldp_dy[0];
    b2 = node_fieldp_mod[0]*node_fieldp_dz[0]*delta_xi[0];
    f0 = 1/(1+pow2s((a0-a1)/a2));
    f1 = 1/(1+pow2s((b0-b1)/b2));

    // Calculate X derivative of the argument
    der_arg0 = (2*delta_eta[n]*node_fieldp_dx[n]-delta_xi[n]*node_fieldp_dy[n])/(node_fieldp_mod[n]*node_fieldp_dz[n]*delta_xi[n]);
    der_arg1 = (2*delta_eta[0]*node_fieldp_dx[0]-delta_xi[0]*node_fieldp_dy[0])/(node_fieldp_mod[0]*node_fieldp_dz[0]*delta_xi[0]);
    velocity[0] += der_arg0*f0 - der_arg1*f1;

    // Calculate Y derivative of the argument
    der_arg0 = -delta_xi[n]*node_fieldp_dx[n]/(node_fieldp_mod[n]*node_fieldp_dz[n]*delta_xi[n]);
    der_arg1 = -delta_xi[0]*node_fieldp_dx[0]/(node_fieldp_mod[0]*node_fieldp_dz[0]*delta_xi[0]);
    velocity[1] += der_arg0*f0 - der_arg1*f1;

    // Calculate Z derivative of the argument
    der_arg0 = 2*pow2s(node_fieldp_dz[n])*delta_eta[n]*node_fieldp_mod[n]*delta_xi[n];
    der_arg0 -= (a0-a1)*node_fieldp_mod[n]*delta_xi[n];
    der_arg0 /= pow2s(node_fieldp_mod[n]*node_fieldp_dz[n]*delta_xi[n]);

    der_arg1 = 2*pow2s(node_fieldp_dz[0])*delta_eta[0]*node_fieldp_mod[0]*delta_xi[0];
    der_arg1 -= (b0-b1)*node_fieldp_mod[0]*delta_xi[0];
    der_arg1 /= pow2s(node_fieldp_mod[0]*node_fieldp_dz[0]*delta_xi[0]);

    velocity[2] += der_arg0*f0 - der_arg1*f1;
}