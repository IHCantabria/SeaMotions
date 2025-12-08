
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
                                        cusfloat    rhow_in,
                                        cusfloat    grav_acc_in,
                                        cusfloat    draft_in,
                                        cusfloat*   cog_in,
                                        cusfloat*   radii_inertia_in,
                                        T*          mesh_in
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
                                        cusfloat    rhow_in,
                                        cusfloat    grav_acc_in,
                                        cusfloat    draft_in,
                                        cusfloat*   cog_in,
                                        cusfloat*   radii_inertia_in,
                                        T*          mesh_in
                    );

    /* Define class public methods */
    void    print( void );

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
    void    recalculate( 
                                        cusfloat    rhow_in,
                                        cusfloat    grav_acc_in,
                                        cusfloat    draft_in,
                                        cusfloat*   cog_in,
                                        cusfloat*   radii_inertia_in,
                                        T*          mesh_in
                            );

};

// Include template class definition
#include "initial_stability.txx"