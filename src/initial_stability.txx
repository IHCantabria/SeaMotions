
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


// Include local modules
#include "initial_stability.hpp"
#include "math/math_tools.hpp"
#include "math/integration.hpp"


template<int NUM_GP, typename T>
cusfloat            InitialStability<NUM_GP, T>::get_area_wl( void ) const
{
    return this->_area_wl;
}


template<int NUM_GP, typename T>
const cusfloat*     InitialStability<NUM_GP, T>::get_area_wl_cog( void ) const
{
    return this->_area_wl_cog;
}


template<int NUM_GP, typename T>
cusfloat            InitialStability<NUM_GP, T>::get_area_wl_mx( void ) const
{
    return this->_area_wl_mx;
}


template<int NUM_GP, typename T>
cusfloat            InitialStability<NUM_GP, T>::get_area_wl_my( void ) const
{
    return this->_area_wl_my;
}


template<int NUM_GP, typename T>
cusfloat            InitialStability<NUM_GP, T>::get_area_wl_ixx( void ) const
{
    return this->_area_wl_ixx;
}


template<int NUM_GP, typename T>
cusfloat            InitialStability<NUM_GP, T>::get_area_wl_ixy( void ) const
{
    return this->_area_wl_ixy;
}


template<int NUM_GP, typename T>
cusfloat            InitialStability<NUM_GP, T>::get_area_wl_iyy( void ) const
{
    return this->_area_wl_iyy;
}


template<int NUM_GP, typename T>
cusfloat            InitialStability<NUM_GP, T>::get_bmx( void ) const
{
    return this->_bmx;
}


template<int NUM_GP, typename T>
cusfloat            InitialStability<NUM_GP, T>::get_bmy( void ) const
{
    return this->_bmy;
}


template<int NUM_GP, typename T>
const cusfloat*     InitialStability<NUM_GP, T>::get_cob( void ) const
{
    return this->_cob;
}


template<int NUM_GP, typename T>
cusfloat            InitialStability<NUM_GP, T>::get_draft( void ) const
{
    return this->_draft;
}


template<int NUM_GP, typename T>
const cusfloat*     InitialStability<NUM_GP, T>::get_eigen_period( void ) const
{
    return this->_eigen_period;
}


template<int NUM_GP, typename T>
cusfloat            InitialStability<NUM_GP, T>::get_gmx( void ) const
{
    return this->_gmx;
}


template<int NUM_GP, typename T>
cusfloat            InitialStability<NUM_GP, T>::get_gmy( void ) const
{
    return this->_gmy;
}


template<int NUM_GP, typename T>
const cusfloat*     InitialStability<NUM_GP, T>::get_gz( void ) const
{
    return this->_gz;
}


template<int NUM_GP, typename T>
const cusfloat*     InitialStability<NUM_GP, T>::get_hydrostiffmat( void ) const
{
    return this->_hydstiffmat;
}


template<int NUM_GP, typename T>
cusfloat            InitialStability<NUM_GP, T>::get_kmx( void ) const
{
    return this->_kmx;
}


template<int NUM_GP, typename T>
cusfloat            InitialStability<NUM_GP, T>::get_kmy( void ) const
{
    return this->_kmy;
}


template<int NUM_GP, typename T>
cusfloat            InitialStability<NUM_GP, T>::get_mass( void ) const
{
    return this->_mass;
}


template<int NUM_GP, typename T>
cusfloat            InitialStability<NUM_GP, T>::get_mch( void ) const
{
    return this->_mch;
}


template<int NUM_GP, typename T>
cusfloat            InitialStability<NUM_GP, T>::get_mct( void ) const
{
    return this->_mct;
}


template<int NUM_GP, typename T>
cusfloat            InitialStability<NUM_GP, T>::get_tpc( void ) const
{
    return this->_tpc;
}


template<int NUM_GP, typename T>
cusfloat            InitialStability<NUM_GP, T>::get_volume( void ) const
{
    return this->_volume;
}


template<int NUM_GP, typename T>
cusfloat            InitialStability<NUM_GP, T>::get_volume_mx( void ) const
{
    return this->_volume_mx;
}


template<int NUM_GP, typename T>
cusfloat            InitialStability<NUM_GP, T>::get_volume_my( void ) const
{
    return this->_volume_my;
}


template<int NUM_GP, typename T>
cusfloat            InitialStability<NUM_GP, T>::get_volume_mz( void ) const
{
    return this->_volume_mz;
}


template<int NUM_GP, typename T>
InitialStability<NUM_GP, T>::InitialStability(
                                                    const cusfloat      rhow_in,
                                                    const cusfloat      grav_acc_in,
                                                    const cusfloat      draft_in,
                                                    const cusfloat*     cog_in,
                                                    const cusfloat*     radii_inertia_in,
                                                    T*                  mesh
                                                )
{
    // Call private class implementation for the recalculation
    // of the hydrostatic values
    this->_recalculate(
                            rhow_in,
                            grav_acc_in,
                            draft_in,
                            cog_in,
                            radii_inertia_in,
                            mesh
                        );
}


template<int NUM_GP, typename T>
void InitialStability<NUM_GP, T>::_process_panel_data(
                                                            PanelGeom*  paneli
                                                        )
{
    // Calculate water plane area
    auto    area_wl_fcn     =   [&]( int ){ return paneli->normal_vec[2]; };

    gauss2d_loop<NUM_GP>(
                            this->_area_wl,
                            area_wl_fcn,
                            paneli
                        );

    // Calculate water plane area first order moments
    auto    area_wl_mx_fcn  =   [&]( int i )
                                { 
                                    return (
                                                paneli->gauss_points_global_y[i]
                                                *
                                                paneli->normal_vec[2]
                                            ); 
                                };

    gauss2d_loop<NUM_GP>(
                            this->_area_wl_mx,
                            area_wl_mx_fcn,
                            paneli
                        );

    auto    area_wl_my_fcn  =   [&]( int i )
                                { 
                                    return (
                                                paneli->gauss_points_global_x[i]
                                                *
                                                paneli->normal_vec[2]
                                            ); 
                                };

    gauss2d_loop<NUM_GP>(
                            this->_area_wl_my,
                            area_wl_my_fcn,
                            paneli
                        );

    // Calculate water plane area second order moments
    auto    area_wl_ixx_fcn  =   [&]( int i )
                                { 
                                    return (
                                                pow2s( paneli->gauss_points_global_y[i] )
                                                *
                                                paneli->normal_vec[2]
                                            ); 
                                };

    gauss2d_loop<NUM_GP>(
                            this->_area_wl_ixx,
                            area_wl_ixx_fcn,
                            paneli
                        );

    auto    area_wl_iyy_fcn  =   [&]( int i )
                                { 
                                    return (
                                                pow2s( paneli->gauss_points_global_x[i] )
                                                *
                                                paneli->normal_vec[2]
                                            ); 
                                };

    gauss2d_loop<NUM_GP>(
                            this->_area_wl_iyy,
                            area_wl_iyy_fcn,
                            paneli
                        );

    // Calculate volume
    auto    volume_fcn      =   [&](int i)
                                { 
                                    return ( 
                                                -paneli->gauss_points_global_z[i] 
                                                * 
                                                paneli->normal_vec[2] 
                                            );
                                };
    
    gauss2d_loop<NUM_GP>(
                                this->_volume,
                                volume_fcn,
                                paneli
                            );

    // Calculate volume moments
    auto    volume_mx_fcn   =   [&](int i)
                                { 
                                    return (
                                                -paneli->gauss_points_global_z[i] 
                                                *  
                                                paneli->gauss_points_global_x[i]
                                                * 
                                                paneli->normal_vec[2] 
                                            ); 
                                };
    
    gauss2d_loop<NUM_GP>(
                                this->_volume_mx,
                                volume_mx_fcn,
                                paneli
                            );

    auto    volume_my_fcn   =   [&](int i)
                                { 
                                    return (
                                                -paneli->gauss_points_global_z[i] 
                                                *  
                                                paneli->gauss_points_global_y[i]
                                                * 
                                                paneli->normal_vec[2] 
                                            ); 
                                };
    
    gauss2d_loop<NUM_GP>(
                                this->_volume_my,
                                volume_my_fcn,
                                paneli
                            );

    auto    volume_mz_fcn   =   [&](int i)
                                { 
                                    return (
                                                -paneli->gauss_points_global_z[i] 
                                                *  
                                                paneli->gauss_points_global_z[i] / 2.0
                                                * 
                                                paneli->normal_vec[2] 
                                            ); 
                                };
    
    gauss2d_loop<NUM_GP>(
                                this->_volume_mz,
                                volume_mz_fcn,
                                paneli
                            );
    
}


/**
 * @brief Print initial stabilty properties on the CLI.
 *
 */
template<int NUM_GP, typename T>
void    InitialStability<NUM_GP, T>::print( 
                                                    void
                                    )
{
    std::cout << std::fixed << std::setprecision( 6 );
    std::cout << std::endl;
    std::cout << "HYDROSTATICS: " << std::endl;
    std::cout << " - Volume properties: "             << std::endl;
    std::cout << "--------------------------------"   << std::endl;
    std::cout << "      + Displaced Volume [m3]:    " << this->_volume << std::endl;
    std::cout << "      + Displacement [kg]:        " << this->_mass << std::endl;
    std::cout << "      + COB [m]:                  " << this->_cob[0] << ", " << this->_cob[1] << ", " << this->_cob[2] << std::endl;
    std::cout << "      + Draft [m]:                " << this->_draft << std::endl;
    std::cout << " - Area properties: " << std::endl;
    std::cout << "--------------------------------"   << std::endl;
    std::cout << "      + WL Area [m2]:             " << this->_area_wl << std::endl;
    std::cout << "      + WL Area Mx [m3]:          " << this->_area_wl_mx << std::endl;
    std::cout << "      + WL Area My [m3]:          " << this->_area_wl_my << std::endl;
    std::cout << "      + WL area Ixx [m4]:         " << this->_area_wl_ixx << std::endl;
    std::cout << "      + WL area Ixy [m4]:         " << this->_area_wl_ixy << std::endl;
    std::cout << "      + WL area Iyy [m4]:         " << this->_area_wl_iyy << std::endl;
    std::cout << " - Mass Properties: "               << std::endl;
    std::cout << "--------------------------------"   << std::endl;
    std::cout << "      + Structural Mass [kg]:     " << this->_mass << std::endl;
    std::cout << "      + COG [m]:                  " << this->_cog[0] << ", " << this->_cog[1] << ", " << this->_cog[2] << std::endl;
    std::cout << " - Stability Properties: "          << std::endl;
    std::cout << "--------------------------------"   << std::endl;
    std::cout << "      + KB [m]:                   " << this->_cob[2] << std::endl;
    std::cout << "      + BMX [m]:                  " << this->_bmx << std::endl;
    std::cout << "      + BMY [m]:                  " << this->_bmy << std::endl;
    std::cout << "      + KMX [m]:                  " << this->_kmx << std::endl;
    std::cout << "      + KMY [m]:                  " << this->_kmy << std::endl;
    std::cout << "      + KG [m]:                   " << this->_cog[2] << std::endl;
    std::cout << "      + GMX [m]:                  " << this->_gmx << std::endl;
    std::cout << "      + GMY [m]:                  " << this->_gmy << std::endl;
    std::cout << " - Hydrostatic Stiffness Matrix: "  << std::endl;
    std::cout << "--------------------------------"   << std::endl;
    print_matrix( 6, 6, this->_hydstiffmat, 6, 0, 1 );
    std::cout << " - Dynamic Properties: "  << std::endl;
    std::cout << "--------------------------------"   << std::endl;
    std::cout << "      + T22 [s]:                  " << this->_eigen_period[0] << std::endl;
    std::cout << "      + T33 [s]:                  " << this->_eigen_period[1] << std::endl;
    std::cout << "      + T44 [s]:                  " << this->_eigen_period[2] << std::endl;
    std::cout << std::endl;

}


template<int NUM_GP, typename T>
void        InitialStability<NUM_GP, T>::_recalculate( 
                                                        const cusfloat      rhow_in,
                                                        const cusfloat      grav_acc_in,
                                                        const cusfloat      draft_in,
                                                        const cusfloat*     cog_in,
                                                        const cusfloat*     radii_inertia_in,
                                                        T*                  mesh_in
                                                    )
{
    // Reset class storage sytem in order to avoid spurious
    // data in case of reuse
    this->_reset( );

    // Storage scalar input arguments
    this->_rhow             = rhow_in;
    this->_grav_acc         = grav_acc_in;
    this->_draft            = draft_in;
    this->_mesh             = mesh_in;

    // Storage vector input arguments
    for ( std::size_t i=0; i<3; i++ )
    {
        this->_cog[i]       = cog_in[i];
        this->_rad_gyr[i]   = radii_inertia_in[i];
    }
    
    // Recalculate geomertrical properties of the 
    // mesh to feed hydrostatics values calculation
    this->_recalculate_geom_props( );

    // Re-calculate hydrostatic properties based
    // on the inputs and the geometrical properties
    // of the mesh
    this->_recalculate_hydro_props( );

}


template<int NUM_GP, typename T>
void    InitialStability<NUM_GP, T>::recalculate( 
                                                    const cusfloat      rhow_in,
                                                    const cusfloat      grav_acc_in,
                                                    const cusfloat      draft_in,
                                                    const cusfloat*     cog_in,
                                                    const cusfloat*     radii_inertia_in,
                                                    T*                  mesh
                                        )
{
    // Call private class implementation for the recalculation
    // of the hydrostatic values
    this->_recalculate(
                            rhow_in,
                            grav_acc_in,
                            draft_in,
                            cog_in,
                            radii_inertia_in,
                            mesh
                        );
    
}


template<int NUM_GP, typename T>
void InitialStability<NUM_GP, T>::_recalculate_geom_props(
                                                                void
                                                            )
{
    // Loop over fully underwater panels to get their properties
    PanelGeom* paneli   = nullptr;
    for ( int i=0; i<this->_mesh->get_elems_np( ); i++ )
    {
        paneli = this->_mesh->get_panel( i );
        if ( 
                paneli->type == 0
                &&
                paneli->location_zone == PANEL_LOC_UW 
            )
        {
            // Process panel data
            this->_process_panel_data( paneli );

        }

    }

}


template<int NUM_GP, typename T>
void    InitialStability<NUM_GP, T>::_recalculate_hydro_props(
                                                                void
                                                            )
{
    // Calcualte mass
    this->_mass             = this->_volume * this->_rhow;

    // Calculate water plane area COG position
    this->_area_wl_cog[0]   = this->_area_wl_my / this->_area_wl;
    this->_area_wl_cog[1]   = this->_area_wl_mx / this->_area_wl;
    this->_area_wl_cog[2]   = 0.0;

    // Calculate underwater volume COG position
    this->_cob[0]           = this->_volume_mx / this->_volume;
    this->_cob[1]           = this->_volume_my / this->_volume;
    this->_cob[2]           = this->_draft + this->_volume_mz / this->_volume;

    // Calculate metacentric radius
    this->_bmx              = this->_area_wl_ixx / this->_volume;
    this->_bmy              = this->_area_wl_iyy / this->_volume;

    // Calculate metacentric height w.r.t the keel
    this->_kmx              = this->_cob[2] + this->_bmx;
    this->_kmy              = this->_cob[2] + this->_bmy;

    // Calculate metacentric height
    this->_gmx              = this->_kmx - this->_cog[2];
    this->_gmy              = this->_kmy - this->_cog[2];

    // Calculate gz arms
    this->_gz[0]            = this->_cob[1] - this->_cog[1];   // Roll
    this->_gz[1]            = this->_cob[0] - this->_cog[0];   // Pitch

    // Calculate tons per cm of inmersion
    this->_tpc              = ( 
                                    this->_rhow 
                                    * 
                                    this->_area_wl 
                                    /
                                    1000            // To tonnes
                                    / 
                                    100             // Per centimiter
                                );

    // Calculate moment to heel / trim 1 deg
    this->_mch              = this->_mass * this->_grav_acc * this->_gmx * std::sin( PI / 180.0 ) / 1000.0;
    this->_mct              = this->_mass * this->_grav_acc * this->_gmy * std::sin( PI / 180.0 ) / 1000.0;

    /* Calculate hydrostatic stiffness */
    // Clean hydrostatic stiffness matrix
    clear_vector( 36, this->_hydstiffmat );

    // K33 - Heave
    this->_hydstiffmat[14] = this->_grav_acc * this->_rhow * this->_area_wl;

    // K34 - K43 - Heave/Roll
    this->_hydstiffmat[15] = this->_grav_acc * this->_rhow * this->_area_wl_my;
    this->_hydstiffmat[20] = this->_grav_acc * this->_rhow * this->_area_wl_my;

    // K35 - K53 - Heave/Pitch
    this->_hydstiffmat[16] = - this->_grav_acc * this->_rhow * this->_area_wl_mx;
    this->_hydstiffmat[26] = - this->_grav_acc * this->_rhow * this->_area_wl_mx;

    // K44 - Roll
    this->_hydstiffmat[21] = this->_grav_acc * this->_rhow * ( this->_area_wl_ixx + this->_volume * ( this->_cob[2] - this->_cog[2] ) );

    // K45 - K54 - Roll/Pitch
    this->_hydstiffmat[22] = this->_grav_acc * this->_rhow * this->_area_wl_ixy;
    this->_hydstiffmat[27] = this->_grav_acc * this->_rhow * this->_area_wl_ixy;

    // K55 - Pitch
    this->_hydstiffmat[28] = this->_grav_acc * this->_rhow * ( this->_area_wl_iyy + this->_volume * ( this->_cob[2] - this->_cog[2] ) );

    // K46 - K56 - Yaw
    this->_hydstiffmat[23] = - this->_grav_acc * this->_rhow * ( this->_cob[0] - this->_cog[0] ) * this->_volume;
    this->_hydstiffmat[29] = - this->_grav_acc * this->_rhow * ( this->_cob[1] - this->_cog[1] ) * this->_volume;

    /* Calculate eigen periods */
    this->_eigen_period[0] = 2.0 * PI * std::sqrt( this->_mass / this->_hydstiffmat[14] );
    this->_eigen_period[1] = 2.0 * PI * std::sqrt( ( this->_mass * pow2s( this->_rad_gyr[0] ) ) / this->_hydstiffmat[21] );
    this->_eigen_period[2] = 2.0 * PI * std::sqrt( ( this->_mass * pow2s( this->_rad_gyr[1] ) ) / this->_hydstiffmat[28] );
}


/**
 * @brief Resets the internal state of the InitialStability object.
 *
 * This method clears or reinitializes any internal variables so the
 * object returns to its default stable state
 */
template<int NUM_GP, typename T>
void    InitialStability<NUM_GP, T>::_reset(
                                                    void
                                )
{
    // Reset scalar variables
    this->_area_wl        = 0.0;  // Water plane area
    this->_area_wl_mx     = 0.0;  // Water plane area first order moemnto around X axis
    this->_area_wl_my     = 0.0;  // Water plane area first order moemnto around Y axis
    this->_area_wl_ixx    = 0.0;  // Water plane area second order moment around X axis
    this->_area_wl_ixy    = 0.0;  // Water plane area second order moment combined X and Y axis
    this->_area_wl_iyy    = 0.0;  // Water plane area second order moment around Y axis
    this->_bmx            = 0.0;  // Metacentric radius aroun X axis
    this->_bmy            = 0.0;  // Metacentric radius aroun Y axis
    this->_draft          = 0.0;  // Draft of the floater. This value is storaged for reference
    this->_gmx            = 0.0;  // Metacentric height aroun X axis
    this->_gmy            = 0.0;  // Metacentric height aroun Y axis
    this->_grav_acc       = 0.0;  // Gravitational acceleration 
    this->_kmx            = 0.0;  // Height of the metacentre over the keel for X axis
    this->_kmy            = 0.0;  // Height of the metacentre over the keel for Y axis
    this->_mass           = 0.0;  // Mass of water displaced by the hull at the required draft
    this->_mch            = 0.0;  // Moment to change heel 1 degree
    this->_mct            = 0.0;  // Moment to change trim 1 degree
    this->_rhow           = 0.0;  // Water density
    this->_tpc            = 0.0;  // Tonnes per cm immersion
    this->_volume         = 0.0;  // Volume of water displaced by the hull at the required draft
    this->_volume_mx      = 0.0;  // Volume of water displaced by the hull at the required draft first order moment X
    this->_volume_my      = 0.0;  // Volume of water displaced by the hull at the required draft first order moment Y
    this->_volume_mz      = 0.0;  // Volume of water displaced by the hull at the required draft first order moment Z

    // Reset to zero vector variables
    clear_vector( 3,  this->_area_wl_cog    );
    clear_vector( 3,  this->_cog            );
    clear_vector( 3,  this->_cob            );
    clear_vector( 3,  this->_eigen_period   );
    clear_vector( 36, this->_hydstiffmat    );

}
