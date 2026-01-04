
// Include local modules
#include "../math/math_interface.hpp"
#include "source.hpp"
#include "../static_tools.hpp"


template<int mode_f, int mode_dfdr, int mode_dfdz>
void    calculate_source_newman_t(
                                            PanelGeom*      panel, 
                                            cusfloat*       field_point, 
                                            int             fp_local_flag,
                                            int             multipole_flag, 
                                            cusfloat*       velocity,
                                            cusfloat&       phi
                                        )
{
    // Create axuliar velocity vector to store the velocities in local coordinates
    cusfloat velocity_local[3] = {0.0, 0.0, 0.0};

    // Reset potential value
    phi = 0.0;

    // Translate field point to panel local coordinates. Check, if multipole
    // formulation applies
    cusfloat field_point_local[3];
    if (fp_local_flag == 0)
    {
        // Calculate vector from center of the panel to field point
        cusfloat field_point_local_aux[3];
        sv_sub(3, field_point, panel->center, field_point_local_aux);
                
        // Calculate distance from the center of the panel to the field point
        cusfloat r0 = cblas_nrm2<cusfloat>(3, field_point_local_aux, 1);

        // Check if multipole expansion applies
        if ((r0/panel->length > 4.0) && multipole_flag)
        {
            calculate_source_monopole<mode_f, mode_dfdr, mode_dfdz>(
                                                                        panel, 
                                                                        r0, 
                                                                        field_point_local_aux,
                                                                        phi,
                                                                        velocity
                                                                    );
            return;
        }
        else if ( r0/panel->length > 2.5 && multipole_flag )
        {
            calculate_source_multipole<mode_f, mode_dfdr, mode_dfdz>(
                                                                        panel, 
                                                                        r0, 
                                                                        field_point_local_aux,
                                                                        phi,
                                                                        velocity
                                                                    );
            return;
        }

        // Change from global system of coordinates to the local one
        cblas_gemv<cusfloat>(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, panel->global_to_local_mat, 3, field_point_local_aux, 1, 0, field_point_local, 1);
    }
    else if (fp_local_flag == 1)
    {
        // Store field point in the local container. This is done to have compatilibity
        // with the case where the field point is provided in global coordinate system.
        field_point_local[0] = field_point[0];
        field_point_local[1] = field_point[1];
        field_point_local[2] = field_point[2];

        // Calculate distance from the center of the panel to the field point
        cusfloat r0 = cblas_nrm2<cusfloat>(3, field_point_local, 1);

        // Check if multipole expansion applies
        if ((r0/panel->length > 4.0) && multipole_flag)
        {
            calculate_source_monopole<mode_f, mode_dfdr, mode_dfdz>(
                                                                        panel, 
                                                                        r0, 
                                                                        field_point_local, 
                                                                        phi,
                                                                        velocity
                                                                    );
            return;
        }
        if ((r0/panel->length > 2.5) && multipole_flag)
        {
            calculate_source_multipole<mode_f, mode_dfdr, mode_dfdz>(
                                                                        panel, 
                                                                        r0, 
                                                                        field_point_local, 
                                                                        phi,
                                                                        velocity
                                                                    );
            return;
        }
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

    cusfloat sides_len [panel->MAX_PANEL_NODES];
    calculate_sides_len_local(panel, delta_xi, delta_eta, sides_len);

    // Calculate polar angles
    cusfloat polar_angles [panel->MAX_PANEL_NODES];
    calculate_polar_angles(panel, delta_xi, delta_eta, polar_angles);

    // Calculate induced velocities
    cusfloat r_sum = 0.0;
    cusfloat b, b0, b1;
    cusfloat a, c, d;
    int i1;
    for (int i=0; i<panel->num_nodes; i++)
    {
        // Calculate forward index
        i1 = (i+1)%panel->num_nodes;

        // Calculate potential-like coefficients
        r_sum   = node_fieldp_mod[i]+node_fieldp_mod[i1];
        b0      = r_sum + sides_len[i];
        b1      = r_sum - sides_len[i]; b1 = check_zero_eps( b1, ZEROTH_EPS );
        b       = std::log(b0/b1); 

        a       = node_fieldp_dx[i]*std::sin(polar_angles[i])-node_fieldp_dy[i]*std::cos(polar_angles[i]);
        c       = std::log((r_sum+sides_len[i])/b1);
        phi     += a*c;


        // Calculate X derivative coefficients
        velocity_local[0] -= delta_eta[i]/sides_len[i]*b;

        // Calculate Y derivative coefficients
        velocity_local[1] += delta_xi[i]/sides_len[i]*b;

    }

    // Calculate potential
    cusfloat phi_dipole = 0.0;
    calculate_dipole_potential_kernel(panel, node_fieldp_mod, node_fieldp_dx, node_fieldp_dy, node_fieldp_dz, delta_xi, delta_eta, phi_dipole);
    phi -= node_fieldp_dz[0]*phi_dipole;
    phi *= -1.0;

    // Calculate Z local velocity
    if ( std::abs( field_point_local[2] ) < FIELD_POINT_LOCAL_TOL )
    {
        velocity_local[2] = 0.0;
    }
    else
    {
        velocity_local[2] = phi_dipole;
    }

    // Convert back the velocity vector to global coordinates if the field
    // vector was given in that base
    if (fp_local_flag == 0)
    {
        cblas_gemv<cusfloat>(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, panel->local_to_global_mat, 3, velocity_local, 1, 0, velocity, 1);
    }
    else
    {
        copy_vector(3, velocity_local, velocity);
    }

    // Invert sign to have the outgoing velocities for positive sources
    svs_mult( 3, velocity, -1.0, velocity );
}


template<int mode_f, int mode_dfdr, int mode_dfdz>
void    calculate_source_monopole(
                                            PanelGeom*      panel, 
                                            cusfloat        r0, 
                                            cusfloat*       field_point,
                                            cusfloat*       velocity,
                                            cusfloat&       phi
                                        )
{
    // Calculate potential and velocity based on the requested modes
    STATIC_COND( ONLY_FCN,      (phi         = panel->area/r0);                         )
    STATIC_COND( ONLY_FCNDR,    (velocity[0] = field_point[0]*panel->area/pow3s(r0));   )
    STATIC_COND( ONLY_FCNDR,    (velocity[1] = field_point[1]*panel->area/pow3s(r0));   )
    STATIC_COND( ONLY_FCNDZ,    (velocity[2] = field_point[2]*panel->area/pow3s(r0));   )
}


template<int mode_f, int mode_dfdr, int mode_dfdz>
void    calculate_source_multipole(
                                            PanelGeom*      panel, 
                                            cusfloat        r0, 
                                            cusfloat*       field_point,
                                            cusfloat*       velocity,
                                            cusfloat&       phi
                                        )
{
    // Load panel properties into local variables to have a 
    // clearer notation along the function
    cusfloat A      = panel->area;
    cusfloat Mx     = panel->moments_fo[0];
    cusfloat My     = panel->moments_fo[1];
    cusfloat Ixx    = panel->moments_so[0];
    cusfloat Iyy    = panel->moments_so[1];
    cusfloat Ixy    = panel->moments_so[2];
    cusfloat x      = field_point[0];
    cusfloat y      = field_point[1];
    cusfloat z      = field_point[2];

    // Calculate inverse radius powers
    cusfloat w      = 1.0 / r0;
    cusfloat w2     = w * w;
    cusfloat w3     = w2 * w;
    cusfloat w4     = w3 * w;
    cusfloat w7     = w4 * w3;

    // Calculate auxiliar powers of field point coordinates
    STATIC_COND( ONLY_FCNDRZ,   (cusfloat x2     = pow2s( x )); )
    STATIC_COND( ONLY_FCNDRZ,   (cusfloat y2     = pow2s( y )); )
    STATIC_COND( ONLY_FCNDRZ,   (cusfloat z2     = pow2s( z )); )

    // Calculate inverse radius convinations
    STATIC_COND( ONLY_FCN,      (cusfloat wx     = - x * w3); )
    STATIC_COND( ONLY_FCN,      (cusfloat wy     = - y * w3); )
    STATIC_COND( ONLY_FCN,      (cusfloat wz     = - z * w3); )
    STATIC_COND( ONLY_FCN,      (cusfloat wxx    = - w3 + 3.0 * x2 * w5); )
    STATIC_COND( ONLY_FCN,      (cusfloat wxy    = + 3.0 * x * y * w5); )
    STATIC_COND( ONLY_FCN,      (cusfloat wyy    = - w3 + 3.0 * y2 * w5); )

    STATIC_COND( ONLY_FCNDRZ,   (cusfloat wxxx   = 3.0 * x * ( 3.0 * p + 10.0 * x2 ) * w7);     )
    STATIC_COND( ONLY_FCNDRZ,   (cusfloat wxxy   = 3.0 * y * p * w7);                           )
    STATIC_COND( ONLY_FCNDRZ,   (cusfloat wxyy   = 3.0 * x * q * w7);                           )
    STATIC_COND( ONLY_FCNDRZ,   (cusfloat wyyy   = 3.0 * y * ( 3.0 * q + 10.0 * y2 ) * w7);     )
    STATIC_COND( ONLY_FCNDRZ,   (cusfloat wxxz   = 3.0 * z * p * w7);                           )
    STATIC_COND( ONLY_FCNDRZ,   (cusfloat wxyz   = -15.0 * x * y * z * w7);                     )
    STATIC_COND( ONLY_FCNDRZ,   (cusfloat wyyz   = 3.0 * z * q * w7);                           )

    // Calculate potential and velocity based on the requested modes
    STATIC_COND( ONLY_FCN,      (phi             = A * w - ( Mx * wx + My * wy ) + 0.5 * ( Ixx * wxx + 2.0 * Ixy * wxy + Iyy * wyy ) ); )
    STATIC_COND( ONLY_FCN,      (velocity[0]     = - ( A * wx + 0.5 * Ixx * wxxx + Ixy * wxxy + 0.5 * Iyy * wxyy ) );                   )
    STATIC_COND( ONLY_FCN,      (velocity[1]     = - ( A * wy + 0.5 * Ixx * wxxy + Ixy * wxyy + 0.5 * Iyy * wyyy ) );                   )
    STATIC_COND( ONLY_FCN,      (velocity[2]     = - ( A * wz + 0.5 * Ixx * wxxz + Ixy * wxyz + 0.5 * Iyy * wyyz ) );                   )

}