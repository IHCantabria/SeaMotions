
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

// Include local modules
#include "../../containers/load_condition.hpp"
#include "../../initial_stability.hpp"
#include "../../mesh/rigid_body_mesh.hpp"
#include "stab_input.hpp"


class GZPoint
{
private:
    /* Declare type aliases for the current class */
    using HSInitStab = InitialStability<NUM_GP, RigidBodyMesh>;

    /* Declare class private attributes */
    int                     _axis_id        = 0;            // Axis identifier to set if GZ is calculated about X or Y axis
    cusfloat                _draft          = 0.0;          // Draft at the current loading condition [m]
    cusfloat                _grav_acc       = 0.0;          // Gravitational acceleration [m/s^2]
    cusfloat                _gz[2]          = { 0.0, 0.0 }; // GZ arm value [m]
    cusfloat                _heel           = 0.0;          // Heeling angle [rad]
    HSInitStab              _hs_final_state ;               // Hydrostatics final equilibrium state for the imposed heeling angle
    const LoadCondition*    _load_cond      = nullptr;      // Pointer to load condition definition
    cusfloat                _mass           = 0.0;          // Mass at the current loading condition [kg]
    RigidBodyMesh*          _mesh           = nullptr;      // Pointer to rigid body mesh
    cusfloat                _water_density  = 0.0;          // Water density [kg/m^3]

    /* Declare class private methods */

    /**
     * @brief   Method to find equilibrium position for the current heeling angle
     * 
     * @param   mass_eq     Equilibrium mass at the current loading condition [kg]
     * @param   y0_pos      Initial position defining the heel and trim for the state to calculate [m]
     * 
     */
    cusfloat    _find_vertical_equilibrium( 
                                                cusfloat    mass_eq,
                                                cusfloat*   y0_pos
                                            );

    /**
     * @brief   Method to find equilibrium position for the current heeling angle
     * 
     * @param   mass_eq     Equilibrium mass at the current loading condition [kg]
     * @param   y0_pos      Initial position defining the heel and trim for the state to calculate [m]
     * @param   min_z       Minimum search limit for the vertical equilibrium position [m]
     * @param   max_z       Maximum search limit for the vertical equilibrium position [m]
     */
    cusfloat    _find_vertical_equilibrium( 
                                                cusfloat    mass_eq,
                                                cusfloat*   y0_pos,
                                                cusfloat    min_z,
                                                cusfloat    max_z
                                            );

    /**
     * @brief   Method to update GZ value for the input load condition and heeling angle
     */
    void        _update( 
                                    const LoadCondition*    load_cond_in,
                                    RigidBodyMesh*          mesh_in,
                                    cusfloat                heel_in,
                                    int                     axis_id_in,
                                    cusfloat                water_density_in,
                                    cusfloat                grav_acc_in
                        );
public:
    /* Declare class constructors */

    /**
     * @brief   Default constructor method for GZPoint class
     */
    GZPoint( ) = default;

    /**
     * @brief   Constructor method for GZPoint class
     * 
     * @param   load_cond_in    Pointer to load condition definition
     * @param   heel_in         Heeling angle [rad]
     */
    GZPoint(
                const LoadCondition*    load_cond_in,
                RigidBodyMesh*          mesh_in,
                cusfloat                heel_in,
                int                     axis_id_in,
                cusfloat                water_density_in,
                cusfloat                grav_acc_in
            );


    /**
     * @brief   Getter method for GZ value for the current heeling angle and loading condition and specified axis
     * 
     * @return  Floating point value representing the GZ arm [m]
     */
    cusfloat   get_gz( void ) const;

    /**
     * @brief   Update method to recalculate GZ value for new input parameters
     * 
     * @param   load_cond_in        Pointer to load condition definition
     * @param   mesh_in             Pointer to rigid body mesh
     * @param   heel_in             Heeling angle [rad]
     * @param   axis_id_in          Axis identifier to set if GZ is calculated about X or Y axis
     * @param   water_density_in    Water density [kg/m^3]
     * @param   grav_acc_in         Gravitational acceleration [m/s^2]
     */
    void        update(
                            const LoadCondition*    load_cond_in,
                            RigidBodyMesh*          mesh_in,
                            cusfloat                heel_in,
                            int                     axis_id_in,
                            cusfloat                water_density_in,
                            cusfloat                grav_acc_in
                        );

};