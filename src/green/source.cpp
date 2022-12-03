
// Include general usage scientific libraries
#include <cmath>
#include "mkl.h"

// Include local modules
#include "common.hpp"
#include "../config.hpp"
#include "../containers.hpp"
#include "dipole.hpp"
#include "source.hpp"


void calculate_source_potential_newman(PanelGeom &panel, cusfloat (&field_point)[3], int fp_local_flag, cusfloat &phi)
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

    cusfloat sides_len [panel.MAX_PANEL_NODES];
    calculate_sides_len_local(panel, delta_xi, delta_eta, sides_len);

    // Calculate polar angles
    cusfloat polar_angles [panel.MAX_PANEL_NODES];
    calculate_polar_angles(panel, delta_xi, delta_eta, polar_angles);

    // Calculate velocity potential
    int i1 = 0;
    cusfloat r_sum = 0.0;
    cusfloat a, b;
    for (int i=0; i<panel.num_nodes; i++)
    {
        // Calculate next index
        i1 = (i+1)%panel.num_nodes;

        // Calculate potential
        r_sum = node_fieldp_mod[i]+node_fieldp_mod[i1];
        a = node_fieldp_dx[i]*std::sin(polar_angles[i])-node_fieldp_dy[i]*std::cos(polar_angles[i]);
        b = std::log((r_sum+sides_len[i])/(r_sum-sides_len[i]));
        phi += a*b;
    }
    
    cusfloat phi_dipole = 0.0;
    calculate_dipole_potential_kernel(panel, node_fieldp_mod, node_fieldp_dx, node_fieldp_dy, node_fieldp_dz, delta_xi, delta_eta, phi_dipole);
    phi -= node_fieldp_dz[0]*phi_dipole;

}


void calculate_source_velocity_newman(PanelGeom &panel, cusfloat (&field_point)[3], int fp_local_flag, cusfloat (&velocity)[3])
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

    cusfloat sides_len [panel.MAX_PANEL_NODES];
    calculate_sides_len_local(panel, delta_xi, delta_eta, sides_len);

    // Calculate induced velocities
    cusfloat r_sum = 0.0;
    cusfloat b, b0, b1;
    int i1;
    for (int i=0; i<panel.num_nodes; i++)
    {
        // Calculate forward index
        i1 = (i+1)%panel.num_nodes;

        // Calculate potential-like coefficients
        r_sum = node_fieldp_mod[i]+node_fieldp_mod[i1];
        b0 = r_sum + sides_len[i];
        b1 = r_sum - sides_len[i];
        b = std::log(b0/b1);


        // Calculate X derivative coefficients
        velocity[0] -= delta_eta[i]/sides_len[i]*b;

        // Calculate Y derivative coefficients
        velocity[1] += delta_xi[i]/sides_len[i]*b;

    }
    calculate_dipole_potential_kernel(panel, node_fieldp_mod, node_fieldp_dx, node_fieldp_dy, node_fieldp_dz, delta_xi, delta_eta, velocity[2]);

}


void calculate_source_potential_hess(PanelGeom &panel, cusfloat (&field_point)[3], int fp_local_flag, cusfloat &phi)
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

    cusfloat sides_len [panel.MAX_PANEL_NODES];
    calculate_sides_len_local(panel, delta_xi, delta_eta, sides_len);

    // Calcualte potential value
    int j;
    cusfloat a, b, si, sj, r;
    cusfloat cx, cy;
    for (int i=0; i<panel.num_nodes; i++)
    {
        // Calculate forward index
        j = (i+1)%panel.num_nodes;

        // Calculate local parameters
        a = delta_xi[i]/sides_len[i];
        b = delta_eta[i]/sides_len[i];
        si = -(a*node_fieldp_dx[i] + b*node_fieldp_dy[i]);
        sj = -(a*node_fieldp_dx[j] + b*node_fieldp_dy[j]);
        r = node_fieldp_dx[i]*b - node_fieldp_dy[i]*a;

        // Calculate potential
        cy = r*std::abs(node_fieldp_dz[i])*(node_fieldp_mod[i]*sj-node_fieldp_mod[j]*si);
        cx = node_fieldp_mod[i]*node_fieldp_mod[j]*pow2s(r)+pow2s(node_fieldp_dz[i])*si*sj;
        phi += std::atan2(cy, cx);

    }

}


void calculate_source_velocity_hess(PanelGeom &panel, cusfloat (&field_point)[3], int fp_local_flag, cusfloat (&velocity)[3])
{
    // Create axuliar velocity vector to store the velocities in local coordinates
    cusfloat velocity_local[3] = {0.0, 0.0, 0.0};

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

    cusfloat sides_len [panel.MAX_PANEL_NODES];
    calculate_sides_len_local(panel, delta_xi, delta_eta, sides_len);

    // Calculate velocities
    cusfloat eki, eki1, hki, hki1, mi, r_sum, log_r;
    int i1 = 0;
    for (int i=0; i<panel.num_nodes; i++)
    {
        // Calculate forward index
        i1 = (i+1)%panel.num_nodes;

        // Calculate local quatities
        eki = pow2s(node_fieldp_dz[i]) + pow2s(node_fieldp_dx[i]);
        eki1 = pow2s(node_fieldp_dz[i1]) + pow2s(node_fieldp_dx[i1]);
        hki = node_fieldp_dy[i]*node_fieldp_dx[i];
        hki1 = node_fieldp_dy[i1]*node_fieldp_dx[i1];
        r_sum = node_fieldp_mod[i]+node_fieldp_mod[i1];
        log_r = std::log((r_sum-sides_len[i])/(r_sum+sides_len[i]));

        // Calculate X velocity
        velocity_local[0] += delta_eta[i]/sides_len[i]*log_r;

        // Calculate Y velocity
        velocity_local[1] -= delta_xi[i]/sides_len[i]*log_r;

        // Calculate Z velocity
        mi = delta_eta[i]/delta_xi[i];
        velocity_local[2] += std::atan((mi*eki-hki)/(node_fieldp_dz[i]*node_fieldp_mod[i]));
        velocity_local[2] -= std::atan((mi*eki1-hki1)/(node_fieldp_dz[i1]*node_fieldp_mod[i1]));
        
    }

    // Convert back the velocity vector to global coordinates if the field
    // vector was given in that base
    if (fp_local_flag == 0)
    {
        cblas_gemv<cusfloat>(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, panel.local_to_global_mat, 3, velocity_local, 1, 0, velocity, 1);
    }
}