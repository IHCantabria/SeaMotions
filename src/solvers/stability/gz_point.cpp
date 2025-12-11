
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
#include "hydrostatic_force_nlin.hpp"
#include "../../interfaces/hydrostatic_force_interface.hpp"
#include "gz_point.hpp"


cusfloat    GZPoint::_find_vertical_equilibrium( 
                                                    cusfloat    mass_eq,
                                                    cusfloat*   y0_pos
                                                )
{
    // Estimate vertical limits for the equilibrium search
    cusfloat    lz      = this->_mesh->z_max - this->_mesh->z_min;

    // Launch vertical equilibrium search
    cusfloat    sol     = this->_find_vertical_equilibrium( 
                                                                mass_eq,
                                                                y0_pos,
                                                                -lz,
                                                                +lz
                                                            );
    return sol;
}


cusfloat    GZPoint::_find_vertical_equilibrium( 
                                                    cusfloat    mass_eq,
                                                    cusfloat*   y0_pos,
                                                    cusfloat    min_z,
                                                    cusfloat    max_z
                                                )
{
    // Create functor to calculate hydrostatic forces
    HydrostaticForceInterface<NUM_GP>   hydrostat_force_interf( 
                                                                    this->_water_density, 
                                                                    this->_grav_acc 
                                                                );

    HydrostaticForcesNLin               hydrostatic_force(          
                                                                    this->_mesh,
                                                                    y0_pos,
                                                                    &hydrostat_force_interf, 
                                                                    mass_eq * this->_grav_acc 
                                                                );

    // Use bisection method to find equilibrium position
    cusfloat    sol     = 0.0;
    int         info    = 0;
    
    bisection_false( 
                        hydrostatic_force, 
                        min_z, 
                        max_z, 
                        0.01, 
                        1e-5,
                        10000, 
                        false, 
                        sol, 
                        info 
                    );

    hydrostatic_force( sol );

    // Calculate GZ based on hydrostatic moments
    this->_gz[0] = hydrostatic_force.get_last_hydrostatic_forces( )[3] / this->_mass / this->_grav_acc;
    this->_gz[1] = hydrostatic_force.get_last_hydrostatic_forces( )[4] / this->_mass / this->_grav_acc;
    
    return sol;
}


        GZPoint::GZPoint(
                            const LoadCondition*    load_cond_in,
                            RigidBodyMesh*          mesh_in,
                            cusfloat                heel_in,
                            int                     axis_id_in,
                            cusfloat                water_density_in,
                            cusfloat                grav_acc_in
                        )
{
    // Update GZ value
    this->_update(
                        load_cond_in,
                        mesh_in,
                        heel_in,
                        axis_id_in,
                        water_density_in,
                        grav_acc_in
                    );

}


cusfloat   GZPoint::get_gz( void ) const
{
    cusfloat gz = 0.0;
    if ( this->_axis_id == 0 )
    {
        gz = this->_gz[0];
    }
    else if ( this->_axis_id == 1 )
    {
        gz = this->_gz[1];
    }

    return gz;
}


void    GZPoint::update(
                            const LoadCondition*    load_cond_in,
                            RigidBodyMesh*          mesh_in,
                            cusfloat                heel_in,
                            int                     axis_id_in,
                            cusfloat                water_density_in,
                            cusfloat                grav_acc_in
                        )
{
    // Update GZ value
    this->_update(
                        load_cond_in,
                        mesh_in,
                        heel_in,
                        axis_id_in,
                        water_density_in,
                        grav_acc_in
                    );
}


void    GZPoint::_update( 
                                    const LoadCondition*    load_cond_in,
                                    RigidBodyMesh*          mesh_in,
                                    cusfloat                heel_in,
                                    int                     axis_id_in,
                                    cusfloat                water_density_in,
                                    cusfloat                grav_acc_in
                        )
{
    // Store input parameters
    this->_axis_id          = axis_id_in;
    this->_grav_acc         = grav_acc_in;
    this->_load_cond        = load_cond_in;
    this->_heel             = heel_in;
    this->_mesh             = mesh_in;
    this->_water_density    = water_density_in;

    // Calculate mass at the current loading condition if any
    this->_mass     = this->_load_cond->mass;
    this->_draft    = this->_load_cond->draft;
    if (  !this->_load_cond->use_mass )
    {
        InitialStability<NUM_GP, RigidBodyMesh> init_stab( 
                                                            this->_water_density,
                                                            this->_grav_acc,
                                                            this->_load_cond->draft,
                                                            this->_load_cond->cog,
                                                            this->_load_cond->rad_inertia,
                                                            this->_mesh
                                                        );
        this->_mass = init_stab.get_mass( );

    }
    else
    {
        // Calculate equivalent draft for the given mass
        cusfloat yh0_pos[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

        this->_draft = std::abs( this->_find_vertical_equilibrium( this->_mass, yh0_pos ) );
    }

    // Define initial position vector
    cusfloat    y0_pos[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    if ( this->_axis_id == 0 )
    {
        y0_pos[3] = this->_heel;
    }
    else if ( this->_axis_id == 1 )
    {
        y0_pos[4] = this->_heel;
    }

    // Find vertical equilibrium for the given heel angle
    this->_find_vertical_equilibrium( 
                                        this->_mass, 
                                        y0_pos,
                                        -1.5 * this->_draft,
                                        -0.5 * this->_draft
                                    );

    this->_hs_final_state   = HSInitStab( 
                                            this->_water_density,
                                            this->_grav_acc,
                                            this->_draft,
                                            this->_mesh->cog,
                                            this->_load_cond->rad_inertia,
                                            this->_mesh
                                        );

}