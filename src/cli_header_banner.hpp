
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
#include <string>

// Include local modules
#include "config.hpp"
#include "version.hpp"


template<bool is_mpi>
inline void cli_header_banner( std::string case_fopath, std::string solver )
{
    /*****************************************/
    /******** Print Header Section ***********/
    /*****************************************/
    const int width         = 80;
    std::string title       = "SeaMotions";
    std::string solver_hd   = "Solver: " + solver;
    std::string version     = "Version " + VERSION_LABEL;

    // Define lambda funciton to center strings
    auto center =   [ & ]
                    ( const std::string& s ) 
                    {
                        int padding = ( width - s.size( ) ) / 2;
                        return std::string( padding, ' ' ) + s;
                    };

    // Check for mpi configuration and root process if necessary
    bool is_root = false;
    if constexpr( is_mpi )
    {
        #include "mpi.h"

        // Get current process rank
        int proc_rank = 0;
        MPI_Comm_rank(
                        MPI_COMM_WORLD,
                        &proc_rank
                    );

        // Check if current process is the root
        is_root = ( proc_rank == MPI_ROOT_PROC_ID );

    }
    else
    {
        is_root = true;

    }

    if ( is_root )
    {
        std::cout << std::endl;
        std::cout << std::string( width , '=') << "\n"
                << center( title )   << "\n"
                << center( solver_hd )   << "\n"
                << center( version ) << "\n"
                << std::string( width , '=') << "\n\n";

        
        // Add licensing banners
        std::cout << std::endl;
        std::cout << "This program is free software: you can redistribute it and/or modify it" << std::endl;
        std::cout << "under the terms of the GNU General Public License v3.0 as published by" << std::endl;
        std::cout << "the Free Software Foundation." << std::endl;
        std::cout << std::endl;
        std::cout << "This program comes WITHOUT ANY WARRANTY; without even the implied" << std::endl;
        std::cout << "warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE." << std::endl;
        std::cout << "See the GNU GPL v3.0 for more details." << std::endl;
        std::cout << "https://www.gnu.org/licenses/gpl-3.0.html" << std::endl << std::endl << std::endl;

        std::cout << " -> Case Path: " << case_fopath << "\n\n";
        std::cout << std::flush;

    }

}