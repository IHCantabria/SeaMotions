
// Include external libraries
#include <iostream>

// Include external scientific libraries
#include <cmath>
#include "mkl.h"

// Include local modules
#include "../containers.hpp"
#include "dipole.hpp"
#include "../math_interface.hpp"
#include "../math_tools.hpp"


cusfloat dipole_potential(PanelGeom &panel, cusfloat (&field_point)[3], int fp_local_flag)
{
    // Translate field point to panel local coordinates
    if (fp_local_flag == 0)
    {
        cusfloat field_point_local_aux[3];
        cusfloat field_point_local[3];
        sv_sub(3, field_point, panel.center, field_point_local_aux);
        cblas_gemv<cusfloat>(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, panel.global_to_local_mat, 3, field_point_local_aux, 1, 0, field_point_local, 1);
    }

    // Calculate distances from each node to the field point in local coordinates
    cusfloat node_fieldp_dx[panel.MAX_PANEL_NODES];
    cusfloat node_fieldp_dy[panel.MAX_PANEL_NODES];
    cusfloat node_fieldp_dz[panel.MAX_PANEL_NODES];
    cusfloat node_fieldp_mod[panel.MAX_PANEL_NODES];
    cusfloat node_fieldp_vec[3];
    cusfloat node_pos[3];
    for (int i=0; i<panel.num_nodes; i++)
    {
        panel.get_node_local_position(i, node_pos);
        sv_sub(3, field_point_local, node_pos, node_fieldp_vec);
        sv_mod(3, node_fieldp_vec, node_fieldp_mod[i]);

        // Storage vector components for futher use
        node_fieldp_dx[i] = node_fieldp_vec[0];
        node_fieldp_dy[i] = node_fieldp_vec[1];
        node_fieldp_dz[i] = node_fieldp_vec[2];
    }

    // Calculate distances in between nodes
    cusfloat delta_xi [panel.MAX_PANEL_NODES];
    cusfloat delta_eta [panel.MAX_PANEL_NODES];

    for (int i=0; i<panel.num_nodes-1; i++)
    {
        delta_xi[i] = panel.xl[i+1] - panel.xl[i];
        delta_eta[i] = panel.yl[i+1] - panel.yl[i];
    }
    delta_xi[panel.num_nodes-1] = panel.xl[0]-panel.xl[panel.num_nodes-1];
    delta_eta[panel.num_nodes-1] = panel.yl[0]-panel.yl[panel.num_nodes-1];

    std::cout << "XL: "; print_vector(panel.num_nodes, panel.xl, 0);
    std::cout << "YL: "; print_vector(panel.num_nodes, panel.yl, 0);
    std::cout << "dXi: "; print_vector(panel.num_nodes, delta_xi, 0);
    std::cout << "dEta: "; print_vector(panel.num_nodes, delta_eta, 0);

    // Calculate potential
    cusfloat phi = 0.0;
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
    a1 = delta_xi[n]*node_fieldp_dx[n]*node_fieldp_dy[n];
    a2 = node_fieldp_mod[n]*node_fieldp_dz[n]*delta_xi[n];
    b0 = delta_eta[n]*(pow2s(node_fieldp_dx[0])+pow2s(node_fieldp_dz[0]));
    b1 = delta_xi[n]*node_fieldp_dx[0]*node_fieldp_dy[0];
    b2 = node_fieldp_mod[0]*node_fieldp_dz[0]*delta_xi[n];
    phi += atan((a0-a1)/a2) - atan((b0-b1)/b2);

    return phi;
}