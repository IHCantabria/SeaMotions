
#include "../config.hpp"
#include "../containers/containers.hpp"
#include "../math/math_tools.hpp"


void calculate_distance_node_field(PanelGeom &panel, cusfloat (&field_point_local)[3], cusfloat* node_fieldp_mod,
    cusfloat* node_fieldp_dx, cusfloat* node_fieldp_dy, cusfloat* node_fieldp_dz)
{
    // Calculate distances from each node to the field point in local coordinates
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
}


void calculate_nodes_distance(PanelGeom &panel, cusfloat* delta_xi, cusfloat* delta_eta)
{
    int i1 = 0;
    for (int i=0; i<panel.num_nodes; i++)
    {
        i1 = (i+1)%panel.num_nodes;
        delta_xi[i] = panel.xl[i1] - panel.xl[i];
        delta_eta[i] = panel.yl[i1] - panel.yl[i];
    }
}


void calculate_polar_angles(PanelGeom &panel, cusfloat* delta_xi, cusfloat* delta_eta, cusfloat* polar_angles)
{
    for (int i=0; i<panel.num_nodes; i++)
    {
        polar_angles[i] = std::atan2(delta_eta[i], delta_xi[i]);
    }
}


void calculate_sides_len_local(PanelGeom &panel, cusfloat* delta_xi, cusfloat* delta_eta, cusfloat* sides_len)
{
    for (int i=0; i<panel.num_nodes; i++)
    {
        sides_len[i] = std::sqrt(pow2s(delta_xi[i])+pow2s(delta_eta[i]));
    }
}