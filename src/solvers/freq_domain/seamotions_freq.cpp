
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
#include "mpi.h"
#include <string>

// Include local modules
#include "../../containers/mpi_timer.hpp"
#include "frequency_solver.hpp"
#include "../../cli_header_banner.hpp"
#include "../../tools/cli_options.hpp"
#include "../../inout/input.hpp"
#include "../../tools.hpp"
#include "../../version.hpp"



int main( int argc, char* argv[] )
{
    // Read command line arguments
    CliOptions cli_options = parse_cli_options( argc, argv );
    const std::string program_name = ( argc > 0 && argv[0] ) ? argv[0] : "seamotions_freq";

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

    std::string case_fopath( cli_options.case_path );

    /*****************************************/
    /****** Initialize MPI environment *******/
    /*****************************************/
    MPI_Init( NULL, NULL );

    // Get total number of processors
    int procs_total = 0;
    MPI_Comm_size(
                    MPI_COMM_WORLD,
                    &procs_total
                );

    // Get current process rank
    int proc_rank = 0;
    MPI_Comm_rank(
                    MPI_COMM_WORLD,
                    &proc_rank
                );

    // Create MPI Configuration system
    MpiConfig mpi_config( proc_rank, procs_total, MPI_ROOT_PROC_ID, MPI_COMM_WORLD );
    MPI_Barrier( MPI_COMM_WORLD );


    /*****************************************/
    /******** Print Header Section ***********/
    /*****************************************/
    cli_header_banner<true>( case_fopath, "Frequency" );


    /*****************************************/
    /************ Read Input data ************/
    /*****************************************/
    Input input;
    input.load( case_fopath );
    input.apply_cli_options( cli_options );
    

    /*****************************************/
    /****** Initialize Frequency Solver ******/
    /*****************************************/
    MpiTimer case_timer;

    FrequencySolver<NUM_GP, PF_ON> freq_solver( &input, &mpi_config );
    
    /*****************************************/
    /******** First Order Solution ***********/
    /*****************************************/
    freq_solver.calculate_first_order( );

    /*****************************************/
    /******** First Order Solution ***********/
    /*****************************************/
    freq_solver.calculate_second_order( );

    /*****************************************/
    /********* Close MPI environment *********/
    /*****************************************/
    case_timer.stop( );
    MPI_Finalize( );

    /*****************************************/
    /**** Close program and final actions ****/
    /*****************************************/

    // Print Elapsed time
    if ( mpi_config.is_root( ) )
    {
        std::cout << std::endl << std::endl;
        std::cout << " -> Seamotions (Frequency) finished!" << std::endl;
        std::cout << " ---> Elapsed wall time for calculation [s]: " << case_timer << std::endl;
    }

    return 0;
}