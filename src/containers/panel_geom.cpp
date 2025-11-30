
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

// Include general usage libraries
#include <cassert>
#include <fstream>

// Include local modules
#include "../math/math_interface.hpp"
#include "../math/math_tools.hpp"
#include "../math/nonlinear_solvers_math.hpp"
#include "../math/shape_functions.hpp"
#include "../math/topology.hpp"
#include "panel_geom.hpp"


void PanelGeom::calcualte_free_surface_singularity( void )
{
    // Clean container to avoid spurious data
    this->free_surface_log_int = 0.0;

    // Declare local variables
    cusfloat A, B, C, df;
    cusfloat a0, a1;
    cusfloat b0, b1;
    cusfloat delta;
    cusfloat int_0, int_1;
    cusfloat p0[3], p1[3];
    cusfloat p0l[3], p1l[3];
    cusfloat p0t[3], p1t[3];
    cusfloat u[3]   = { 0.0, 0.0, 0.0 };
    cusfloat v[3]   = { 0.0, 0.0, 0.0 };
    cusfloat v_norm = 0.0;
    cusfloat w[3]   = { 0.0, 0.0, 1.0 };
    
    // Loop over number of nodes to calculate the different
    // chunks of the integral
    int i1 = 0;
    for ( int i=0; i<this->num_nodes; i++ )
    {
        // Get forward node index
        i1 = ( i + 1 ) % this->num_nodes;

        // Get points for the study
        p0[0] = this->x[i];
        p0[1] = this->y[i];
        p0[2] = this->z[i];
        p1[0] = this->x[i1];
        p1[1] = this->y[i1];
        p1[2] = this->z[i1];

        // Calculate v base vector
        sv_sub( 3, p1, p0, v );
        v_norm = std::sqrt( pow2s( v[0] ) + pow2s( v[1] ) + pow2s( v[2] ) );
        svs_div( 3, v, v_norm, v );

        // Calculate u vector
        cross( v, w, u );

        // Convert to local coordinates p0 and p1
        sv_sub( 3, p0, this->center, p0t );
        sv_sub( 3, p1, this->center, p1t );

        p0l[0] = u[0] * p0t[0] + u[1] * p0t[1];
        p0l[1] = v[0] * p0t[0] + v[1] * p0t[1];
        p1l[0] = u[0] * p1t[0] + u[1] * p1t[1];
        p1l[1] = v[0] * p1t[0] + v[1] * p1t[1];

        // Calculate alpha angle of cylindrical coordinates
        a0 = std::atan2( p0l[1], p0l[0] );
        a1 = std::atan2( p1l[1], p1l[0] );

        // Calculate A and B factors
        df  = p0l[0] * p1l[1] - p0l[1] * p1l[0];
        A   = ( p1l[1] - p0l[1] ) / df;
        B   = ( p0l[0] - p1l[0] ) / df;
        C   = std::sqrt( pow2s( A ) + pow2s( B ) );

        // Calculate delta angle
        delta = std::atan2( B, A );

        // Calculate beta angles
        b0 = a0 - delta;
        b1 = a1 - delta;

        // Calculate integral value
        int_0 = - std::tan( b0 ) / 2.0 / pow2s( C ) * ( std::log( C ) + std::log( std::cos( b0 ) ) + 1.5 ) + b0/ 2.0 / pow2s( C );
        int_1 = - std::tan( b1 ) / 2.0 / pow2s( C ) * ( std::log( C ) + std::log( std::cos( b1 ) ) + 1.5 ) + b1/ 2.0 / pow2s( C );
        
        this->free_surface_log_int += int_1 - int_0;
    }
}


// Add method to calculate the geometric propertiess
void PanelGeom::calculate_properties( 
                                        cusfloat* cog
                                    )
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
        v1_aux[0] = this->x[1] - this->x[3];
        v1_aux[1] = this->y[1] - this->y[3];
        v1_aux[2] = this->z[1] - this->z[3];
    }
    else
    {
        throw std::runtime_error("Number of vertexes specified is not correct.");
    }

    // Calculate normal vector to the panel
    cross(v1_aux, v0, v2);

    // Calculate base vector 1 in a orthogonal base
    cross(v2, v0, v1);

    // Normalize base vectors
    cusfloat v0_mod = cblas_nrm2<cusfloat>(3, v0, 1);
    cusfloat v1_mod = cblas_nrm2<cusfloat>(3, v1, 1);
    cusfloat v2_mod = cblas_nrm2<cusfloat>(3, v2, 1);
    
    svs_mult(3, v0, 1/v0_mod, v0);
    svs_mult(3, v1, 1/v1_mod, v1);
    svs_mult(3, v2, 1/v2_mod, v2);

    // Calculate area of the panel based on the norm of the cross product
    this->area = v2_mod/2.0;

    // Calcualte the maximum representative distance of the panel
    if (this->num_nodes == 3)
    {
        // Calculate vector that joins the third side
        cusfloat v_star[3];
        v_star[0] = this->x[2] - this->x[1];
        v_star[1] = this->y[2] - this->y[1];
        v_star[2] = this->z[2] - this->z[1];

        // Calculate modulus of the second and third sides
        cusfloat v1_aux_mod = cblas_nrm2<cusfloat>(3, v1_aux, 1);
        cusfloat v_star_mod = cblas_nrm2<cusfloat>(3, v_star, 1);

        // Calculate panel length
        this->length = v0_mod;
        if (v1_aux_mod > this->length)
        {
            this->length = v1_aux_mod;
        }
        
        if (v_star_mod > this->length)
        {
            this->length = v_star_mod;
        }

    }
    else if (this->num_nodes == 4)
    {
        // Calculate modulus of the second diagonal of the panel
        cusfloat v1_aux_mod = cblas_nrm2<cusfloat>(3, v1_aux, 1);

        // Calculate panel length
        this->length = v0_mod;

        if (v1_aux_mod > this->length)
        {
            this->length = v1_aux_mod;
        }

    }

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

    // Storage normal vector to the panel
    copy_vector( 3, v2, this->normal_vec );
    svs_mult( 3, this->normal_vec, -1.0, this->normal_vec );

    // Take the local reference system centre
    this->sysref_centre[0] = this->x[0];
    this->sysref_centre[1] = this->y[0];
    this->sysref_centre[2] = this->z[0];

    // Declare auxiliar vectors to perform the vector rotations
    cusfloat global_pos[3];
    cusfloat local_pos[3];

    // Calculate local coordinates of the panel.
    for (int i=0; i<this->num_nodes; i++)
    {
        // Remove mean point of the panel in order to rotate the panel
        // around a point inside of it
        global_pos[0] = this->x[i] - this->sysref_centre[0];
        global_pos[1] = this->y[i] - this->sysref_centre[1];
        global_pos[2] = this->z[i] - this->sysref_centre[2];

        // Rotate node position to express the node coordinates in
        // the local coordinate system
        cblas_gemv<cusfloat>(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, this->global_to_local_mat, 3, global_pos, 1, 0, local_pos, 1);

        // Storage of the solution in the local x,y,z storage vectors
        this->xl[i] = local_pos[0];
        this->yl[i] = local_pos[1];
        this->zl[i] = local_pos[2];
    }

    for (int i=0; i<this->num_nodes; i++)
    {
        // Remove mean point of the panel in order to rotate the panel
        // around a point inside of it
        global_pos[0] = this->x[i] - this->center[0];
        global_pos[1] = this->y[i] - this->center[1];
        global_pos[2] = this->z[i] - this->center[2];

        // Rotate node position to express the node coordinates in
        // the local coordinate system
        cblas_gemv<cusfloat>(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, this->global_to_local_mat, 3, global_pos, 1, 0, local_pos, 1);

        // Storage of the solution in the local x,y,z storage vectors
        this->xlc[i] = local_pos[0];
        this->ylc[i] = local_pos[1];
        this->zlc[i] = local_pos[2];
    }

    // Calculate normal vector components for rotational modes
    // around the COG of the body
    copy_vector( 3, cog, this->body_cog );

    cusfloat cog_to_panel[3] = { 0.0, 0.0, 0.0 };
    sv_sub( 
                3, 
                this->center, 
                cog, 
                cog_to_panel 
            );
    
    this->normal_vec[3] = (
                                cog_to_panel[1]*this->normal_vec[2]
                                -
                                cog_to_panel[2]*this->normal_vec[1]
                            );
    this->normal_vec[4] = (
                                cog_to_panel[2]*this->normal_vec[0]
                                -
                                cog_to_panel[0]*this->normal_vec[2]
                            );
    this->normal_vec[5] = (
                                cog_to_panel[0]*this->normal_vec[1]
                                -
                                cog_to_panel[1]*this->normal_vec[0]
                            );

}


void PanelGeom::calculate_source_nodes(
                                            int         poly_order,
                                            cusfloat*    
                                        )
{
    // Define sources data
    if ( poly_order == 0 )
    {
        // Allocate heap memory to storage the source points data
        this->_source_positions     = generate_empty_vector<cusfloat>( 3 );
        this->_source_normal_vec    = generate_empty_vector<cusfloat>( 6 );

        // Define source points data
        copy_vector( 3, this->center, this->_source_positions );
        copy_vector( 6, this->normal_vec, this->_source_normal_vec );

    }
    else if ( poly_order > 0 )
    {
        std::cout << "Polynomial orders higher than 0 not implemented yet!" << std::endl;
        throw std::runtime_error( "" );
        // Allocate heap memory to storage the source points data
        // this->source_positions = new

        // Define source points data
    }
    else
    {
        std::cerr << "ERROR - PanelGeom" << std::endl;
        std::cerr << "Polynomial order not valid." << std::endl;
        throw std::runtime_error( "" );
    }

    this->_is_source_nodes = true;
}


void PanelGeom::get_node_local_position( int num_node, cusfloat* node_pos )
{
    node_pos[0] = this->xl[num_node];
    node_pos[1] = this->yl[num_node];
    node_pos[2] = this->zl[num_node];
}


void PanelGeom::get_node_local_position_c( int num_node, cusfloat* node_pos )
{
    node_pos[0] = this->xlc[num_node];
    node_pos[1] = this->ylc[num_node];
    node_pos[2] = this->zlc[num_node];
}


void PanelGeom::get_node_position( int num_node, cusfloat* node_pos )
{
    node_pos[0] = this->x[num_node];
    node_pos[1] = this->y[num_node];
    node_pos[2] = this->z[num_node];
}


void PanelGeom::get_panel_xy_proj( 
                                    PanelGeom* new_panel 
                                )
{
    // Copy number of nodes
    new_panel->num_nodes = this->num_nodes;

    // Copy XY nodes projection
    for ( int i=0; i<this->num_nodes; i++ )
    {
        new_panel->x[i] = this->x[i];
        new_panel->y[i] = this->y[i];
        new_panel->z[i] = 0.0;
    }

    // Calculate panel properties
    new_panel->calculate_properties( this->body_cog );
    
}


void PanelGeom::get_source_nodes_data(
                                        cusfloat*& position,
                                        cusfloat*& normals_vec
                                    )
{
    assert( this->_is_source_nodes && "Thre is no source nodes definition to take from PanelGeom class." );
    position    = this->_source_positions;
    normals_vec = this->_source_normal_vec;
}


int PanelGeom::is_inside( cusfloat* field_point )
{
    // Check if the field point is inside of the bounding box
    // of the panel.
    cusfloat x_min = 1e16;
    cusfloat x_max = -1e16;
    cusfloat y_min = 1e16;
    cusfloat y_max = -1e16;

    for (int i=0; i<this->num_nodes; i++)
    {
        // Check bounding box limits in X axis
        if (this->xl[i]>x_max)
        {
            x_max = this->xl[i];
        }

        if (this->xl[i]<x_min)
        {
            x_min = this->xl[i];
        }

        // Check bounding box limits in Y axis
        if (this->yl[i]>y_max)
        {
            y_max = this->yl[i];
        }

        if (this->yl[i]<y_min)
        {
            y_min = this->yl[i];
        }
    }

    if ((field_point[0]>x_max) || (field_point[0]<x_min) || (field_point[1]>y_max) || (field_point[1]<y_min))
    {
        return 0;
    }

    // Look if the point is also inside of the panel
    cusfloat lam, mu;
    cusfloat rx = x_min-2*std::abs(x_min);
    cusfloat ry = field_point[1];
    cusfloat uiy = 0.0;
    cusfloat uix = 0.0;
    cusfloat vx = field_point[0] - rx;
    int count_cross = 0;
    int j = 0;
    for (int i=0; i<this->num_nodes; i++)
    {
        // Calculate forward index
        j = (i+1)%num_nodes;

        // Calculate side vector coefficients
        uix = this->xl[j] - this->xl[i];
        uiy = this->yl[j] - this->yl[i];

        // Calculate vectors scale
        mu = (ry-this->yl[i])/uiy;
        lam = (this->xl[i]+mu*uix-rx)/vx;

        // Check mu value to control precision
        if (std::abs(mu)<EPS_PRECISION)
        {
            mu = 0.0;
        }

        // Check if there is cross in between the ray traced and 
        // the side of the panel
        if ((lam>0.0) & (lam<1.0) & (mu>=0.0) & (mu<1.0))
        {
            count_cross += 1;
        }

    }
    
    int is_inside = 0;
    if ((count_cross%2)==1)
    {
        is_inside = 1;
    }

    return is_inside;
}


void PanelGeom::local_to_global( 
                                    cusfloat    xi, 
                                    cusfloat    eta,
                                    cusfloat*   global_pos
                                )
{
    // Generate shape functions value container
    cusfloat N[MAX_LIN_NODES]; clear_vector( MAX_LIN_NODES, N );

    // Get shape functions value
    shape_fcn_2d( this->num_nodes, xi, eta, N );

    // Get global coordinates
    // cusfloat x2d[3]         = { 0.0, 0.0, 0.0 };
    this->_x2d[0]                  = cblas_dot<cusfloat>( this->num_nodes, this->xl, 1, N, 1 );
    this->_x2d[1]                  = cblas_dot<cusfloat>( this->num_nodes, this->yl, 1, N, 1 );
    this->_x2d[2]                  = 0.0;

    clear_vector( 3, this->_global_pos );
    clear_vector( 3, global_pos );
    cblas_gemv<cusfloat>(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, this->local_to_global_mat, 3, this->_x2d, 1, 0, this->_global_pos, 1);

    sv_add( 3, this->_global_pos, this->sysref_centre, global_pos );

}


void PanelGeom::local_coords_from_z_proj(
                                            cusfloat    x,
                                            cusfloat    y,
                                            cusfloat&   xi,
                                            cusfloat&   eta
                                        )
{
    // Define function for the optimization
    auto aux_fcn    =   [ x, y, this ]
                        (int* , int*, cusfloat* x_in, cusfloat* x_out )
                        {
                            // Get global coordinates
                            cusfloat global_pos[3] = { 0.0, 0.0, 0.0};
                            this->local_to_global( 
                                                    x_in[0],
                                                    x_in[1],
                                                    global_pos
                                                );

                            // Calculate difference between the target 
                            // and the estimated value
                            x_out[0] = x - global_pos[0];
                            x_out[1] = y - global_pos[1];
                        };

    // Define initial paramets
    cusfloat sol0[2] = { 0.0, 0.0 };

    // Define output vector
    cusfloat sol[2] = { 0.0, 0.0 };
    
    // Use Trust-Region algorithm to estimate the local coordinates
    int status = 0;
    trust_region(
                    aux_fcn,
                    2,
                    2,
                    sol0,
                    sol,
                    status,
                    false
                );

    // Copy output values
    xi  = sol[0];
    eta = sol[1];

}


PanelGeom::PanelGeom(
                        int         npe,
                        cusfloat*   x_in,
                        cusfloat*   y_in,
                        cusfloat*   z_in,
                        int         is_move_f_in,
                        int         type_in,
                        cusfloat*   cog,
                        bool        force_auto_type
                    )
{
    // Get panel vertexes form element list
    this->num_nodes = npe;
    for ( int j=0; j<npe; j++ )
    {
        this->x[j]   = x_in[j];
        this->y[j]   = y_in[j];
        this->z[j]   = z_in[j];
    }

    // Set panel movility
    this->is_move_f = is_move_f_in;

    // Calculate panel properties
    this->calculate_properties( cog );

    // Set panel type. It is interpreted based on the distance to the FS
    this->type = type_in;
    if ( force_auto_type )
    {
        if ( std::abs( this->center[2] ) < FS_SEL_THR )
        {
            this->type = LID_PANEL_CODE;
        }
        else
        {
            this->type = DIFFRAC_PANEL_CODE;
        }
    }

    // Calculate integral properties
    this->calculate_integration_properties<NUM_GP>( );

    // Calculate free surface singularity if necessary
    if ( this->type == LID_PANEL_CODE )
    {
        this->calcualte_free_surface_singularity( );
    }
}


PanelGeom::~PanelGeom(
                        void
                    )
{
    if ( this->_is_source_nodes )
    {
        mkl_free( this->_source_positions );
        mkl_free( this->_source_normal_vec );
    }
}


void PanelGeom::write(
                        std::string finame
                    )
{
    // Open file unit
    std::ofstream outfile( finame );

    // Write number of panel nodes
    outfile << "NumNodes: " << this->num_nodes << std::endl;

    // Write nodes X position
    outfile << "NodesXPos   NodesYPos   NodesZPos" << std::endl;
    for ( int i=0; i<this->num_nodes; i++ )
    {
        outfile << this->x[i] << " " << this->y[i] << " " << this->z[i] << std::endl;
    }

    // Close file unit
    outfile.close( );

}


std::ostream& operator<< ( std::ostream& os, PanelGeom& panel )
{
    os << "PANEL PROPERTIES:" << std::endl;
    os << " - PANEL COORDINATES: " << std::endl;
    for ( int i=0; i<panel.num_nodes; i++ )
    {
        os << "   -> X: " << panel.x[i] << " Y: " << panel.y[i] << " Z: " << panel.z[i] << std::endl;
    }
    os << std::endl;

    os << " - CENTROID:" << std::endl;
    os << "   -> X: " << panel.center[0] << " Y: " << panel.center[1] << " Z: " << panel.center[2] << std::endl;
    os << std::endl;

    os << " - NORMAL VECTOR:" << std::endl;
    os << "   -> X: " << panel.normal_vec[0] << " Y: " << panel.normal_vec[1] << " Z: " << panel.normal_vec[2] << std::endl;
    os << std::endl;

    os << " - x2d:" << std::endl;
    os << "   -> X: " << panel._x2d[0] << " Y: " << panel._x2d[1] << " Z: " << panel._x2d[2] << std::endl;
    os << std::endl;

    os << " - SysCentre:" << std::endl;
    os << "   -> X: " << panel.sysref_centre[0] << " Y: " << panel.sysref_centre[1] << " Z: " << panel.sysref_centre[2] << std::endl;
    os << std::endl;

    os << " - global_pos:" << std::endl;
    os << "   -> X: " << panel._global_pos[0] << " Y: " << panel._global_pos[1] << " Z: " << panel._global_pos[2] << std::endl;
    os << std::endl;

    os << " - Gauss Points Global:" << std::endl;
    os << "   -> X: " << panel.gauss_points_global_x[0] << " Y: " << panel.gauss_points_global_y[0] << " Z: " << panel.gauss_points_global_z[0] << std::endl;
    os << std::endl;

    return os;
}