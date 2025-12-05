
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
#include "mpi.h"


struct MpiConfig
{
private:
    // Define class private attributes
    bool        _is_parallel    = true;

public:
    // Define class attributes
    MPI_Comm    mpi_comm;
    int         proc_rank   = 0;
    int         procs_total = 0;
    int         proc_root   = 0;

    // Define class constructors and destructor
    MpiConfig( void );

    MpiConfig(
                    int proc_rank_in,
                    int procs_total_in,
                    int proc_root_in,
                    MPI_Comm mpi_comm_in
                );

    // Define class methods
    void    get_1d_bounds(
                            int     np,
                            int&    start_pos,
                            int&    end_pos
                        );

    bool    is_root(    
                            void 
                    );

    void    set_parallel(
                            void
                        );

    void    set_serial( 
                            void 
                        );

};
