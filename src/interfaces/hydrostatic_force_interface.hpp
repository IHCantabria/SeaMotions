
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
#include "../containers/panel_geom.hpp"
#include "../math/integration.hpp"
#include "../math/math_tools.hpp"


template<int NGP>
struct HydrostaticForceInterface
{
private:
    /* Define class private attriutes */
    cusfloat    _grav_acc   = 0.0;
    cusfloat    _rhow       = 0.0;

public:
    /* Define class public attributes */
    cusfloat    fx[ NGP*NGP ];
    cusfloat    fy[ NGP*NGP ];
    cusfloat    fz[ NGP*NGP ];
    cusfloat    mx[ NGP*NGP ];
    cusfloat    my[ NGP*NGP ];
    cusfloat    mz[ NGP*NGP ];

    /* Define class constructor */
    HydrostaticForceInterface( 
                                    cusfloat water_density_in, 
                                    cusfloat grav_acc_in 
                                )
    {
        this->_rhow     = water_density_in;
        this->_grav_acc = grav_acc_in;
    }

    /* Define operator overloading */
    void    operator() (
                            PanelGeom*  panel,
                            cusfloat*   global_force
                        )
    {
        // Declare local variables to assist during
        // loop iteration
        cusfloat _fm = 0.0;
        cusfloat _f[3];
        cusfloat _m[3];
        cusfloat _r[3];

        for ( int i=0; i< ( NGP * NGP ); i++ )
        {
            if( panel->gauss_points_global_z[i] < 0.0 )
            {
                // Calculate force modulus
                _fm     = - this->_rhow * this->_grav_acc * panel->gauss_points_global_z[i];

                // Calculate force components using panel normal
                _f[0]   = _fm * panel->normal_vec[0];
                _f[1]   = _fm * panel->normal_vec[1];
                _f[2]   = _fm * panel->normal_vec[2];

                // Calculate distance from the body cog to the points in the panel
                _r[0]   = ( panel->gauss_points_global_x[i] - panel->body_cog[0] );
                _r[1]   = ( panel->gauss_points_global_y[i] - panel->body_cog[1] );
                _r[2]   = ( panel->gauss_points_global_z[i] - panel->body_cog[2] );
                
                // Calculate moment due to forces at the integration points
                cross( _r, _f, _m );

                // Sum contributions into global forces
                this->fx[i] = _f[0];
                this->fy[i] = _f[1];
                this->fz[i] = _f[2];
                this->mx[i] = _m[0];
                this->my[i] = _m[1];
                this->mz[i] = _m[2];

            }

        }

        // Integrate all forces and moments
        gauss2d_loop<NGP>( 
                                global_force[0], 
                                [&](int i){ return this->fx[i]; },
                                panel 
                            );
        
        gauss2d_loop<NGP>( 
                                global_force[1], 
                                [&](int i){ return this->fy[i]; },
                                panel 
                            );

        gauss2d_loop<NGP>( 
                                global_force[2],
                                [&](int i){ return this->fz[i]; },
                                panel 
                            );

        gauss2d_loop<NGP>( 
                                global_force[3], 
                                [&](int i){ return this->mx[i]; },
                                panel 
                            );

        gauss2d_loop<NGP>( 
                                global_force[4], 
                                [&](int i){ return this->my[i]; },
                                panel 
                            );

        gauss2d_loop<NGP>( 
                                global_force[5], 
                                [&](int i){ return this->mz[i]; },
                                panel 
                            );

    }

};