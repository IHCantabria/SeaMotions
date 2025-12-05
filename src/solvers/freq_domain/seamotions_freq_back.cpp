
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
#include "mpi.h"
#include <string>

// Include local modules
#include "../../containers/mpi_config.hpp"
#include "freq_solver_tools.hpp"
#include "../../hydrostatics.hpp"
#include "../../inout/output.hpp"
#include "../../inout/reader.hpp"
#include "../../tools.hpp"



int main( int argc, char* argv[] )
{
    // Read command line arguments
    if (!check_num_cmd_args(argc, 1))
    {
        return 1;
    }

    std::string case_fopath(argv[1]);

    /*****************************************/
    /************ Read Input data ************/
    /*****************************************/
    Input* input = read_input_files( case_fopath );

    /*****************************************/
    /****** Initialize MPI environment *******/
    /*****************************************/

    // Init MPI
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

    // Declare container to hold MPI configuration
    MpiConfig mpi_config( 
                            proc_rank,
                            procs_total,
                            0,
                            MPI_COMM_WORLD
                        );


    /*****************************************/
    /**** Launch Hydrostatics Calculation ****/
    /*****************************************/
    MPI_Barrier( MPI_COMM_WORLD );
    double hydrostat_tstart     = MPI_Wtime( );
    double case_tstart          = MPI_Wtime( );
    
    Hydrostatics** hydrostatics = new Hydrostatics*[input->bodies_np];
    for ( int i=0; i<input->bodies_np; i++ )
    {
        hydrostatics[i] =   new Hydrostatics( 
                                                input->bodies[i]->mesh,
                                                input->water_density,
                                                input->grav_acc,
                                                input->bodies[i]->mass,
                                                input->bodies[i]->cog,
                                                input->bodies[i]->rad_inertia,
                                                &mpi_config
                                            );
    }
    MPI_Barrier( MPI_COMM_WORLD );
    double hydrostat_tend = MPI_Wtime( );

    if ( mpi_config.is_root( ) )
    {
        std::cout << "Execution time [s]: " << ( hydrostat_tend - hydrostat_tstart ) << std::endl;
    }

    /*****************************************/
    /********* Launch Output System **********/
    /*****************************************/
    Output* output = nullptr;
    
    if ( mpi_config.is_root( ) )
    {
        output = new Output( input );
    }

    /*****************************************/
    /****** Storage Initial Parameters *******/
    /*****************************************/

    if ( mpi_config.is_root( ) )
    {
        // Storage frequency set
        output->save_frequencies( input->freqs );

        // Storage headings set
        output->save_headings( input->heads.data( ) );

        // Storage structural mass
        if ( input->out_struct_mass )
        {
            output->save_structural_mass( );
        }

        // Storage hydrostatic stiffness matrix
        if ( input->out_hydstiff )
        {
            output->save_hydstiffness( hydrostatics );
        }

        // Storage mesh
        if (  input->out_mesh )
        {
            output->save_mesh( );
        }

    }

    /*****************************************/
    /***** Calculate Source Distribution *****/
    /*****************************************/
    calculate_freq_domain_coeffs(
                                    &mpi_config,
                                    input, 
                                    hydrostatics,
                                    output 
                                );

    /*****************************************/
    /********* Close MPI environment *********/
    /*****************************************/
    double case_tend = MPI_Wtime( );
    MPI_Finalize( );

    /*****************************************/
    /**** Delete heap memory allocations *****/
    /*****************************************/
    delete input;
    delete output;

    for ( int i=0; i<input->bodies_np; i++ )
    {
        delete hydrostatics[i];
    }
    delete [] hydrostatics;


    if ( mpi_config.is_root( ) )
    {
        std::cout << "Elapsed wall time for calculation [s]: " << case_tend - case_tstart << std::endl;
        std::cout << "Seamotions Freqcuency ended!" << std::endl;
    }

    return 0;
}