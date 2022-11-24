
#ifndef __containers_hpp
#define __containers_hpp


// Include local modules
#include "config.hpp"
#include "math_interface.hpp"
#include "math_tools.hpp"


struct PanelGeom
{
    cusfloat area = 0.0;
    cusfloat center[3];
    cusfloat global_to_local_mat[9];
    cusfloat local_to_global_mat[9];
    int num_nodes = 0;
    static constexpr int MAX_PANEL_NODES = 4;
    cusfloat x[MAX_PANEL_NODES];
    cusfloat xl[MAX_PANEL_NODES];
    cusfloat y[MAX_PANEL_NODES];
    cusfloat yl[MAX_PANEL_NODES];
    cusfloat z[MAX_PANEL_NODES];
    cusfloat zl[MAX_PANEL_NODES];

    // Add method to calculate the geometric propertiess
    void calculate_properties(void)
    {
        // Calculate ceter of the panel
        for (int i=0; i<this->num_nodes; i++)
        {
            this->center[0] += this->x[i];
            this->center[1] += this->y[i];
            this->center[2] += this->z[i];
        }
        this->center[0] /= this->num_nodes;
        this->center[1] /= this->num_nodes;
        this->center[2] /= this->num_nodes;

        // Calculate base vectors
        cusfloat v0[3] = {0.0, 0.0, 0.0};
        cusfloat v1[3] = {0.0, 0.0, 0.0};
        cusfloat v1_aux[3] = {0.0, 0.0, 0.0};
        cusfloat v2[3] = {0.0, 0.0, 0.0};
        if (this->num_nodes == 3)
        {
            // Define first vector of the base
            v0[0] = this->x[2] - this->x[0];
            v0[1] = this->y[2] - this->y[0];
            v0[2] = this->z[2] - this->z[0];

            // Define second vector of the base. This one is auxiliar because
            // the vertexes of the triangle could not conform an orthonormal base
            v1_aux[0] = this->x[1] - this->x[0];
            v1_aux[1] = this->y[1] - this->y[0];
            v1_aux[2] = this->z[1] - this->z[0];
        }
        else if (this->num_nodes == 4)
        {
            // Define first vector of the base
            v0[0] = this->x[2] - this->x[0];
            v0[1] = this->y[2] - this->y[0];
            v0[2] = this->z[2] - this->z[0];

            // Define second vector of the base. This one is auxiliar because
            // the vertexes of the quadrilateral could not conform an orthonormal base
            v1_aux[0] = this->x[3] - this->x[1];
            v1_aux[1] = this->y[3] - this->y[1];
            v1_aux[2] = this->z[3] - this->z[1];
        }
        else
        {
            throw std::runtime_error("Number of vertexes specified is not correct.");
        }

        // Calculate normal vector to the panel
        cross(v0, v1_aux, v2);

        // Calculate base vector 1 in a orthogonal base
        cross(v2, v0, v1);

        // Normalize base vectors
        double v0_mod = cblas_nrm2<cusfloat>(3, v0, 1);
        double v1_mod = cblas_nrm2<cusfloat>(3, v1, 1);
        double v2_mod = cblas_nrm2<cusfloat>(3, v2, 1);
        
        svs_mult(3, v0, 1/v0_mod, v0);
        svs_mult(3, v1, 1/v1_mod, v1);
        svs_mult(3, v2, 1/v2_mod, v2);

        // Calculate area of the panel based on the norm of the cross product
        this->area = v2_mod/2.0;

        // Generate translation matrixes
        for (int i=0; i<3; i++)
        {
            // Generate global to local coordinate sytem matrix
            this->global_to_local_mat[i] = v0[i];
            this->global_to_local_mat[3+i] = v1[i];
            this->global_to_local_mat[6+i] = v2[i];

            // Generate local to global coordinate system matrix
            this->local_to_global_mat[3*i] = v0[i];
            this->local_to_global_mat[3*i+1] = v1[i];
            this->local_to_global_mat[3*i+2] = v2[i];
        }

        // Declare auxiliar vectors to perform the vector rotations
        cusfloat global_pos[3];
        cusfloat local_pos[3];
        // Calculate local coordinates of the panel.
        for (int i=0; i<this->num_nodes; i++)
        {
            // Remove mean point of the panel in order to rotate the panel
            // around a point inside of it
            this->xl[i] = this->x[i] - this->center[0];
            this->yl[i] = this->y[i] - this->center[1];
            this->zl[i] = this->z[i] - this->center[2];

            // Storage the nodes global coordinate positions into the input vector
            global_pos[0] = this->xl[i];
            global_pos[1] = this->yl[i];
            global_pos[2] = this->zl[i];

            // Rotate node position to express the node coordinates in
            // the local coordinate system
            cblas_gemv<cusfloat>(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, this->global_to_local_mat, 3, global_pos, 1, 0, local_pos, 1);

            // Storage of the solution in the local x,y,z storage vectors
            this->xl[i] = local_pos[0];
            this->yl[i] = local_pos[1];
            this->zl[i] = local_pos[2];
        }

    }
};

#endif