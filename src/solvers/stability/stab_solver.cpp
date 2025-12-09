
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
#include "../../mesh/rigid_body_mesh.hpp"
#include "stab_solver.hpp"


void    StabSolver::calculate_hydrostatics(
                                            void
                                        )
{
    // Define local variables
    std::size_t draft_np    = this->_input->draft_hs.size( );
    std::size_t heel_np     = this->_input->heel_hs_rad.size( );
    cusfloat    _pos[6]     = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    // Check if output is enabled
    if ( this->_input->out_hs )
    {
        // Loop over axes
        for ( std::size_t i=0; i<2; i++ )
        {
            // Loop over heelings
            for ( std::size_t j=0; j<heel_np; j++ )
            {
                // Loop over drafts
                for ( std::size_t k=0; k<draft_np; k++ )
                {
                    // Clear last state to avoid spurious results
                    clear_vector<cusfloat, 6>( _pos );
                    
                    // Update postiion vector
                    _pos[2]     = -this->_input->draft_hs[k];
                    _pos[3+i]   = this->_input->heel_hs_rad[j];

                    // Move mesh to the prescribed postion
                    this->_mesh->move( _pos[0], _pos[1], _pos[2], _pos[3], _pos[4], _pos[5] );

                    // Recalculate mesh intersection
                    this->_mesh->check_underwater_panels( );

                    // Calculate hydrostatic values
                    this->_hydrostatics[j*draft_np+k].recalculate(
                                                                    this->_input->water_density,
                                                                    this->_input->grav_acc,
                                                                    this->_input->draft_hs[k],
                                                                    this->_mesh->cog,
                                                                    this->_input->rad_gyr,
                                                                    this->_mesh
                                                                );
                }
            }

            // Storage axis data
            this->_output->save_hydrostatics( i, this->_hydrostatics );

        }
    }
}


void    StabSolver::_initialize(
                                            void
                                )
{
    // Allocate space for the interal storage system for hydrostatic calculation
    this->_hydrostatics.resize( this->_input->draft_hs.size( ) * this->_input->heel_hs_rad.size( ) );

    // Alloate space for gz results
    this->_gz.resize( this->_input->heel_gz_rad.size( ) );

}

        StabSolver::StabSolver(
                                            StabInput*  input_in
                                        )
{
    // Storage input arguments
    this->_input    = input_in;

    // Create rigid body dynamics mesh
    this->_mesh     = new RigidBodyMesh( 
                                            this->_input->mesh_fipath, 
                                            this->_input->body_name, 
                                            this->_input->cog, 
                                            false, 
                                            0, 
                                            0.0
                                        );

    // Create output system
    this->_output   = new StabOutput( this->_input );

    // Initialize class attributes
    this->_initialize( );

}


        StabSolver::~StabSolver(
                                            void
                               )
{
    // Delete rigid body dynamic mesh
    delete this->_mesh;

    // Delete output system
    delete this->_output;
}