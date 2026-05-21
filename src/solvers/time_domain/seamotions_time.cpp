
/*
 * Copyright (c) 2025 Sergio Fernández Ruano / IHCantabria
 *
 * This file is part of SeaMotionsTimeDev.
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
#include <string>

// Include local modules
#include "../../cli_header_banner.hpp"
#include "../../config.hpp"
#include "../../containers/mpi_config.hpp"
#include "../../containers/mpi_timer.hpp"
#include "../../tools/cli_options.hpp"
#include "../../version.hpp"
#include "input_t.hpp"
#include "time_solver.hpp"

// MPI
#include "mpi.h"


int main( int argc, char* argv[] )
{
    /*****************************************/
    /****** Parse Command-line Options *******/
    /*****************************************/
    CliOptions cli_options = parse_cli_options( argc, argv );
    const std::string program_name = ( argc > 0 && argv[0] ) ? argv[0] : "seamotions_time";

    if ( cli_options.show_help )
    {
        print_cli_usage( std::cout, program_name );
        return 0;
    }

    if ( !cli_options.error.empty( ) )
    {
        std::cerr << cli_options.error << std::endl;
        print_cli_usage( std::cerr, program_name );
        return 1;
    }

    if ( cli_options.case_path.empty( ) )
    {
        print_cli_usage( std::cerr, program_name );
        return 1;
    }

    const std::string case_fopath( cli_options.case_path );

    /*****************************************/
    /****** Initialize MPI environment *******/
    /*****************************************/
    MPI_Init( NULL, NULL );

    int procs_total = 0;
    MPI_Comm_size( MPI_COMM_WORLD, &procs_total );

    int proc_rank = 0;
    MPI_Comm_rank( MPI_COMM_WORLD, &proc_rank );

    // Time-domain solver runs on a single process (no domain decomposition)
    MpiConfig mpi_config( proc_rank, procs_total, MPI_ROOT_PROC_ID, MPI_COMM_WORLD );
    MPI_Barrier( MPI_COMM_WORLD );


    /*****************************************/
    /******** Print Header Section ***********/
    /*****************************************/
    if ( mpi_config.is_root( ) )
    {
        cli_header_banner<true>( case_fopath, "Time" );
    }


    /*****************************************/
    /************ Read Input data ************/
    /*****************************************/
    InputT input;
    input.load( case_fopath );


    /*****************************************/
    /******* Initialize Time Solver **********/
    /*****************************************/
    MpiTimer case_timer;

    TimeSolver<NUM_GP, NUM_GP> time_solver( &input, &mpi_config );


    /*****************************************/
    /******** Run Time-domain Solution *******/
    /*****************************************/
    time_solver.run( );


    /*****************************************/
    /******** Close MPI environment **********/
    /*****************************************/
    case_timer.stop( );
    MPI_Finalize( );


    /*****************************************/
    /**** Close program and final actions ****/
    /*****************************************/
    if ( mpi_config.is_root( ) )
    {
        std::cout << std::endl << std::endl;
        std::cout << " -> Seamotions (Time) finished!" << std::endl;
        std::cout << " ---> Elapsed wall time for calculation [s]: " << case_timer << std::endl;
    }

    return 0;
}
