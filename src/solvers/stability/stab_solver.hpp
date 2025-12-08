
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
#include "../../initial_stability.hpp"
#include "stab_input.hpp"
#include "stab_output.hpp"


class StabSolver
{
    /* Class local type aliases */
    using IS = InitialStability<NUM_GP, RigidBodyMesh>;

private:
    /* Define private attributes */
    std::vector<cusfloat>   _gz             ;           // Vector storaging GZ results for an axis
    std::vector<IS>         _hydrostatics   ;           // Vector storaging hydrostatic results for an axis
    StabInput*              _input          = nullptr;  // Pointer to Stability simulations input system
    RigidBodyMesh*          _mesh           = nullptr;  // Pointer to rigid body dynamic mesh
    StabOutput*             _output         = nullptr;  // Pointer to Stability simulations output system
    
    /* Define private methods */

    /**
     * @brief   Method to initialize class attributes
     */
    void    _initialize( void );

public:
    /* Define class constructors */
    StabSolver( 
                        StabInput*  input_in
                );

    ~StabSolver(
                        void
                );

    /* Define class public methods */
    
    /**
     * @brief   Calculate hydrostatic values for different drafts and headings and stores them into results file
     * 
     * This method allows to calculate variables such as water plane area, its centre and its first and second order moments; the 
     * submerged volume and its centre of gravity; metacentric radius, tons per inmmersion centimiter.
     */
    void    calculate_hydrostatics( void );

};