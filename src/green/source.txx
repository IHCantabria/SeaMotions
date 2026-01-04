
// Include local modules
#include "source.hpp"
#include "../static_tools.hpp"


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
    STATIC_COND( ONLY_FCN,      (phi = panel->area/r0);                                 )
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