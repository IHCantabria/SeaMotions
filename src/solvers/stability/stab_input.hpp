
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
#include <vector>

// Include local modules
#include "../../config.hpp"
#include "../../containers/load_condition.hpp"


class StabInput
{
private:
    /* Define private attributes */
    std::string  _fopath;                        // Folder path where the case configuration is available and where the results are storaged
    std::string  _finame = "sm_stab.input.dat";  // Case file name for the configuration of the simulation with SeaMotions Stability Solver
    std::string  _fipath;                        // Case file path for the configuration of the simulation with SeaMotions Stability Sovler

    /* Define private methods */
    void    _initialize(
                                void
                        );
    
    void    _read_input_file( 
                                std::string fipath 
                            );
    
public:
    /* Define public attributes */
    std::string                 body_name           ;           // Name of the body to read from the mesh definition.
    std::vector<cusfloat>       draft_hs            ;           // Storage draft points to be computed in hydrostatics module
    cusfloat                    grav_acc            = 0.0;      // Gravitational acceleration
    std::vector<cusfloat>       heel_gz_deg         ;           // Storage heeling points to be computes in GZ curves module. Value storaged in degrees.
    std::vector<cusfloat>       heel_gz_rad         ;           // Storage heeling points to be computes in GZ curves module. Value storaged in radians.
    std::vector<cusfloat>       heel_hs_deg         ;           // Storage heeling points to be computes in hydrostatics module. Value storaged in degrees.
    std::vector<cusfloat>       heel_hs_rad         ;           // Storage heeling points to be computes in hydrostatics module. Value storaged in radians.
    bool                        is_bodies           = false;    // Flag to check that the body definitions where read correctly
    std::vector<LoadCondition>  load_conds          ;           // Storage input loading conditions description
    std::vector<std::string>    load_conds_finame   ;           // Storage input loading conditions file names
    std::vector<std::string>    load_conds_name     ;           // Storage input loading conditions names
    std::string                 mesh_finame         ;           // Mesh file name containing the geometrical description of the body to analyze
    std::string                 mesh_fipath         ;           // Mesh file path containing the geometrical description of the body to analyze
    bool                        out_eq              ;           // Flag to switch the output of the equilibrium algorithm
    bool                        out_gz              ;           // Flag to switch the output of GZ curves
    bool                        out_hs              ;           // Flag to switch the output of hydrostatic values
    std::vector<std::string>    tag_eq              ;           // Storage of equilibrium tag cases names
    cusfloat                    water_density       = 0.0;      // Water density of the stability test location


    /* Define class constructors */
    StabInput( std::string fipath );

    /* Define class public methods */

    /**
     * @brief   Getter for case folder path internal variable
     */
    std::string     get_case_fopath( void ) const;

};