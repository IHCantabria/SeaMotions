
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

#pragma once

// Include general usage libraries
#include <iostream>

// Include local modules
#include "config.hpp"
#include "containers/panel_geom.hpp"


template<int NUM_GP, typename T>
struct InitialStability
{
private:
    /* Declare class private attributes */
    cusfloat    _area_wl        = 0.0;      // Water plane area
    cusfloat    _area_wl_mx     = 0.0;      // Water plane area first order moemnto around X axis
    cusfloat    _area_wl_my     = 0.0;      // Water plane area first order moemnto around Y axis
    cusfloat    _area_wl_cog[3] ;           // Water plane area center of gravity
    cusfloat    _area_wl_ixx    = 0.0;      // Water plane area second order moment around X axis
    cusfloat    _area_wl_ixy    = 0.0;      // Water plane area second order moment combined X and Y axis
    cusfloat    _area_wl_iyy    = 0.0;      // Water plane area second order moment around Y axis
    cusfloat    _bmx            = 0.0;      // Metacentric radius aroun X axis
    cusfloat    _bmy            = 0.0;      // Metacentric radius aroun Y axis
    cusfloat    _cog[3]         ;           // Center of gravity position w.r.t to the keel
    cusfloat    _cob[3]  ;                  // Center of gravity of the volume of water displaced by the hull at the required draft
    cusfloat    _draft          = 0.0;      // Draft of the floater. This value is storaged for reference
    cusfloat    _eigen_period[3];           // Eigen periods: Heave, Roll and Pitch
    cusfloat    _gmx            = 0.0;      // Metacentric height aroun X axis
    cusfloat    _gmy            = 0.0;      // Metacentric height aroun Y axis
    cusfloat    _grav_acc       = 0.0;      // Gravitational acceleration
    cusfloat    _gz[2]          ;           // Restoring lever arms for current loading condition. GZ[0] -> Roll, GZ[1] -> Pitch
    cusfloat    _hydstiffmat[36];           // Hydrostatic stiffness matrix 
    cusfloat    _kmx            = 0.0;      // Height of the metacentre over the keel for X axis
    cusfloat    _kmy            = 0.0;      // Height of the metacentre over the keel for Y axis
    cusfloat    _mass           = 0.0;      // Mass of water displaced by the hull at the required draft
    cusfloat    _mch            = 0.0;      // Moment to change heel 1 degree
    cusfloat    _mct            = 0.0;      // Moment to change trim 1 degree
    T*          _mesh           = nullptr;  // Pointer to mesh object containing the mesh description
    cusfloat    _rad_gyr[3];                // Radius of gyration for the three principal axes: X, Y, Z
    cusfloat    _rhow           = 0.0;      // Water density
    cusfloat    _tpc            = 0.0;      // Tonnes per cm immersion
    cusfloat    _volume         = 0.0;      // Volume of water displaced by the hull at the required draft
    cusfloat    _volume_mx      = 0.0;      // Volume of water displaced by the hull at the required draft first order moment X
    cusfloat    _volume_my      = 0.0;      // Volume of water displaced by the hull at the required draft first order moment Y
    cusfloat    _volume_mz      = 0.0;      // Volume of water displaced by the hull at the required draft first order moment Z

    /* Declare private class methods */

    /**
     * @brief   Calculates geometrical properties over target panel
     */
    void    _process_panel_data(
                                        PanelGeom* paneli
                                );

    /**
     * @brief   Public interface to recalculate hydrostatic values for new inputs without destroying the instance.
     * 
     * @param   Water density
     * @param   Gravitational acceleration
     * @param   Floater draft
     * @param   COG position of the floater relative to the center of coordiates
     * @param   Radius of inertial for the floater. Optional for natural periods prediction
     * @param   Mesh instance pointer to calculate hydrostatic properties at the requested loading condition
     */
    void    _recalculate(
                                        const cusfloat      rhow_in,
                                        const cusfloat      grav_acc_in,
                                        const cusfloat      draft_in,
                                        const cusfloat*     cog_in,
                                        const cusfloat*     radii_inertia_in,
                                        T*                  mesh_in
                        );

    /**
     * @brief   Method to wrap up calculation of geometrical mesh properties
     */
    void    _recalculate_geom_props(
                                        void
                                    );


    /**
     * @brief   Method to calculate hydrostatic properties. It should be executed after geometric properties recalculation
     */
    void    _recalculate_hydro_props(
                                        void
                                    );

    /**
     * @brief   Reset state and values of the class.
     */
    void    _reset( 
                                        void            
                    );

public:
    /* Define class constructors */
    InitialStability( ) = default;

    /**
     * @brief   Calculates hydrostatic properties for the input mesh and loading condition
     * 
     * @param   Water density
     * @param   Gravitational acceleration
     * @param   Floater draft
     * @param   COG position of the floater relative to the center of coordiates
     * @param   Radius of inertial for the floater. Optional for natural periods prediction
     * @param   Mesh instance pointer to calculate hydrostatic properties at the requested loading condition
     */
    InitialStability( 
                        const cusfloat      rhow_in,
                        const cusfloat      grav_acc_in,
                        const cusfloat      draft_in,
                        const cusfloat*     cog_in,
                        const cusfloat*     radii_inertia_in,
                        T*                  mesh_in
                    );

    /* Define class public methods */

    /**
     * @brief   Print out on CLI the main hydrostatic properties calculated during the analysis.
     */
    void                print( void );

    /**
     * @brief   Public interface to recalculate hydrostatic values for new inputs without destroying the instance.
     * 
     * @param   Water density
     * @param   Gravitational acceleration
     * @param   Floater draft
     * @param   COG position of the floater relative to the center of coordiates
     * @param   Radius of inertial for the floater. Optional for natural periods prediction
     * @param   Mesh instance pointer to calculate hydrostatic properties at the requested loading condition
     */
    void                recalculate( 
                                        const cusfloat      rhow_in,
                                        const cusfloat      grav_acc_in,
                                        const cusfloat      draft_in,
                                        const cusfloat*     cog_in,
                                        const cusfloat*     radii_inertia_in,
                                        T*                  mesh_in
                                    );

    /**
     * @brief   Getter method for water plane area
     * 
     * @return  Area in m²
     */
    cusfloat            get_area_wl( void ) const;

    /**
     * @brief   Getter method for water plane area C.O.G
     * 
     * @return  Vector of size 3 containing X, Y, Z coordinates of the water plane area C.O.G
     */
    const cusfloat*     get_area_wl_cog( void ) const;

    /**
     * @brief   Getter method for water plane area first order moment around X axis
     * 
     * @return  Moment of area in m³.
     */
    cusfloat            get_area_wl_mx( void ) const;

    /**
     * @brief   Getter method for water plane area first order moment around Y axis
     * 
     * @return  Moment of area in m³.
     */
    cusfloat            get_area_wl_my( void ) const;

    /**
     * @brief   Getter method for water plane area second order moment around X axis
     * 
     * @return  Inertia of area in m4.
     */
    cusfloat            get_area_wl_ixx( void ) const;

    /**
     * @brief   Getter method for water plane area second order moment combined X and Y axis
     * 
     * @return  Inertia of area in m4.
     */
    cusfloat            get_area_wl_ixy( void ) const;

    /**
     * @brief   Getter method for water plane area second order moment around Y axis
     * 
     * @return  Inertia of area in m4.
     */
    cusfloat            get_area_wl_iyy( void ) const;

    /**
     * @brief   Getter method for metacentric radius around X axis
     * 
     * @return  Floating point value representing the metacentric radius around X axis.
     */
    cusfloat            get_bmx( void ) const;

    /**
     * @brief   Getter method for metacentric radius around Y axis
     * 
     * @return  Floating point value representing the metacentric radius around Y axis.
     */
    cusfloat            get_bmy( void ) const;

    /**
     * @brief   Getter method for centre of buoyancy
     * 
     * @return  Vector of size 3 containing X, Y, Z coordinates of the underwater volume C.O.G
     */
    const cusfloat*     get_cob( void ) const;

    /**
     * @brief   Getter method for the eigen periods
     * 
     * @return  Floater draft of the current load condition
     */
    cusfloat            get_draft( void ) const;

    /**
     * @brief   Getter method for the floater eigen periods for the current loading condition
     * 
     * @return  Vector of size 3 containing: Heave, Roll, Pitch eigen periods measured in seconds
     */
    const cusfloat*     get_eigen_period( void ) const;

    /**
     * @brief   Getter method for metacentric height around X axis
     * 
     * @return  Floating point value representing the metacentric height around X axis.
     */
    cusfloat            get_gmx( void ) const;

    /**
     * @brief   Getter method for metacentric height around Y axis
     * 
     * @return  Floating point value representing the metacentric height around Y axis.
     */
    cusfloat            get_gmy( void ) const;

    /**
     * @brief   Getter method GZ arms for the current loading condition
     * 
     * @return  Vector of size 3 containing the level arms. GZ[0] -> Roll, GZ[1] -> Pitch
     */
    const cusfloat*     get_gz( void ) const;

    /**
     * @brief   Getter method for the hydrostatic stiffness matrix.
     * 
     * @return  Vector of size 36 (6x6 row major) containing each coefficient of the matrix: c00, c01, c02, ..., cn0, cn1, ..., cnn.
     */
    const cusfloat*     get_hydrostiffmat( void ) const;

    /**
     * @brief   Getter method for the total metacentric height (Vertical position of Metacentre above the keel) around X axis
     * 
     * @return  Floating point value representing the total metacentric height around X axis.
     */
    cusfloat            get_kmx( void ) const;

    /**
     * @brief   Getter method for the total metacentric height (Vertical position of Metacentre above the keel) around X axis
     * 
     * @return  Floating point value representing the total metacentric height around Y axis.
     */
    cusfloat            get_kmy( void ) const;

    /**
     * @brief   Getter method for the mass of the underwater volume
     * 
     * @return  Floating point value representing the mass of the underwater volume in kilograms
     */
    cusfloat            get_mass( void ) const;

    /**
     * @brief   Getter method for the moment to change heel 1 degree
     * 
     * @return  Floating point value representing the mch in tonnes·meters
     */
    cusfloat            get_mch( void ) const;

    /**
     * @brief   Getter method for the moment to change trim 1 degree
     * 
     * @return  Floating point value representing the mct in tonnes·meters
     */
    cusfloat            get_mct( void ) const;

    /**
     * @brief   Getter method for the tonnes per 1 cm of inmersion
     * 
     * @return  Floating point value representing the tpc in tonnes
     */
    cusfloat            get_tpc( void ) const;

    /**
     * @brief   Getter method for the underwater volume
     * 
     * @return  Floating point value representing the underwater volume in m³
     */
    cusfloat            get_volume( void ) const;

    /**
     * @brief   Getter method for the underwater volume first order moment around X axis
     * 
     * @return  Floating point value representing the underwater volume in m4.
     */
    cusfloat            get_volume_mx( void ) const;

    /**
     * @brief   Getter method for the underwater volume first order moment around Y axis
     * 
     * @return  Floating point value representing the underwater volume in m4.
     */
    cusfloat            get_volume_my( void ) const;

    /**
     * @brief   Getter method for the underwater volume first order moment around Z axis
     * 
     * @return  Floating point value representing the underwater volume in m4.
     */
    cusfloat            get_volume_mz( void ) const;

};

// Include template class definition
#include "initial_stability.txx"