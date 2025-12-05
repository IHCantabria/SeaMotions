
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
#include "global_static_matrix.hpp"


void    calculate_global_hydstiffness(
                                                Input*          input,
                                                Hydrostatics**  hydrostatics,
                                                cusfloat*       hydstiffness
                                    )
{
    int index   = 0;
    for ( int i=0; i<input->bodies_np; i++ )
    {
        for ( int j=0; j<input->dofs_np; j++ )
        {
            for ( int k=0; k<input->dofs_np; k++ )
            {
                index               =   (
                                            i * input->bodies_np *  pow2s( input->dofs_np ) 
                                            + 
                                            i * input->dofs_np 
                                            + 
                                            j * input->dofs_np * input->bodies_np
                                            +
                                            k
                                        );
                hydstiffness[index] = hydrostatics[i]->hydstiffmat[j*input->dofs_np+k];
            }
        }
    }
}



void    calculate_global_structural_mass(
                                                Input*          input,
                                                cusfloat*       structural_mass_p0
                                        )
{
    // Allocate space for the body ith structural mass
    // matrix
    cusfloat*   body_mass   = generate_empty_vector<cusfloat>( pow2s( input->dofs_np ) );
    int index = 0;
    for ( int i=0; i<input->bodies_np; i++ )
    {
        // Clear body matrix to not get spurious data from 
        // the previous body definition
        clear_vector( pow2s( input->dofs_np ), body_mass );

        // Define body mass matrix
        body_mass[0] = input->bodies[i]->mass;  // Surge
        body_mass[7] = input->bodies[i]->mass;  // Sway
        body_mass[14] = input->bodies[i]->mass; // Heave

        if ( input->bodies[i]->interia_by_rad )
        {
            body_mass[21] = input->bodies[i]->mass * pow2s( input->bodies[i]->rad_inertia[0] ); // Roll
            body_mass[28] = input->bodies[i]->mass * pow2s( input->bodies[i]->rad_inertia[1] ); // Pitch
            body_mass[35] = input->bodies[i]->mass * pow2s( input->bodies[i]->rad_inertia[2] ); // Yaw
        }
        else
        {
            body_mass[21] = input->bodies[i]->inertia[0]; // Roll
            body_mass[22] = input->bodies[i]->inertia[1]; // Roll - Pitch
            body_mass[23] = input->bodies[i]->inertia[2]; // Roll - Yaw
            body_mass[27] = input->bodies[i]->inertia[1]; // Pitch - Roll
            body_mass[28] = input->bodies[i]->inertia[3]; // Pitch
            body_mass[29] = input->bodies[i]->inertia[4]; // Pitch - Yaw
            body_mass[33] = input->bodies[i]->inertia[2]; // Yaw - Roll
            body_mass[34] = input->bodies[i]->inertia[4]; // Yaw - Pitch-
            body_mass[35] = input->bodies[i]->inertia[5]; // Yaw
        }

        // Copy ith body structural mass matrix to the global matrix assembly
        for ( int j=0; j<input->dofs_np; j++ )
        {
            for ( int k=0; k<input->dofs_np; k++ )
            {
                index                       =   (
                                                    i * input->bodies_np *  pow2s( input->dofs_np ) 
                                                    + 
                                                    i * input->dofs_np 
                                                    + 
                                                    j * input->dofs_np * input->bodies_np
                                                    +
                                                    k
                                                );
                structural_mass_p0[index]   = body_mass[j*input->dofs_np+k];
            }
        }
    }

    // Deallocate heap memory space for variables in the
    // current function
    mkl_free( body_mass );
}