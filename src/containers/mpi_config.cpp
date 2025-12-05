
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

// Include general usage libraries
#include <iostream>

// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "../config.hpp"
#include "mpi_config.hpp"


void    MpiConfig::get_1d_bounds(
                                    int     np,
                                    int&    start_pos,
                                    int&    end_pos
                                )
{
    if ( this->_is_parallel )
    {
        // Check if there is more processors than data
        if ( this->procs_total > np )
        {
            std::cerr << "There is more processes available than ";
            std::cerr << "data to distribute over them." << std::endl;
            throw std::runtime_error( "" );
        }

        // Divide data in chunks
        int chunk_size = static_cast<int>( 
                                            std::ceil( 
                                                        static_cast<cusfloat>( np ) 
                                                        / 
                                                        static_cast<cusfloat>( this->procs_total ) 
                                                    ) 
                                        ) ;

        // Set interval bounds
        start_pos   = this->proc_rank * chunk_size;
        end_pos     = ( this->proc_rank + 1 ) * chunk_size;

        // Set limits to the upper bound
        end_pos     = ( end_pos > np ) ? np : end_pos;
    }
    else
    {
        start_pos   = 0;
        end_pos     = np;
    }
}


bool    MpiConfig::is_root( void )
{
    return this->proc_root == this->proc_rank;
}


MpiConfig::MpiConfig( void )
{
    // Get total number of processors
    int _procs_total = 0;
    MPI_Comm_size(
                    MPI_COMM_WORLD,
                    &_procs_total
                );

    // Get current process rank
    int _proc_rank = 0;
    MPI_Comm_rank(
                    MPI_COMM_WORLD,
                    &_proc_rank
                );

    // Storage MPI configuration
    this->mpi_comm      = MPI_COMM_WORLD;
    this->proc_rank     = _proc_rank;
    this->procs_total   = _procs_total;
    this->proc_root     = 0;
}


MpiConfig::MpiConfig(
                        int proc_rank_in,
                        int procs_total_in,
                        int proc_root_in,
                        MPI_Comm mpi_comm_in
                    )
{
    this->mpi_comm      = mpi_comm_in;
    this->proc_rank     = proc_rank_in;
    this->procs_total   = procs_total_in;
    this->proc_root     = proc_root_in;
}


void MpiConfig::set_parallel( void )
{
    this->_is_parallel = true;
}


void MpiConfig::set_serial( void )
{
    this->_is_parallel = false;
}