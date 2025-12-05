
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
#include <iomanip>

// Include local module
#include "../../containers/logger.hpp"
#include "../../containers/mpi_timer.hpp"
#include "diffraction.hpp"
#include "frequency_solver.hpp"
#include "froude_krylov.hpp"
#include "global_static_matrix.hpp"
#include "hydromechanics.hpp"
#include "../../hydrostatics.hpp"
#include "../../math/math_tools.hpp"
#include "../../mesh/mesh.hpp"
#include "raos.hpp"


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::calculate_first_order( void )
{
    // Loop over frequencies
    for ( int i=0; i<this->input->angfreqs_np; i++ )
    {
        LOG_TASK_SS( "FO - Calculating period: " << std::setw( 10 ) << std::fixed << std::setprecision( 3 ) << angfreq_to_period( this->input->angfreqs[i] ) << " s" )
        // Start time measurement for the current frequency
        MpiTimer freq_timer;

        // Solve radiation-diffraction problem for the current frequency
        this->kernel->solve( input->angfreqs[i] );
        this->kernel->update_results( sim_data );

        // Calculate added mass and damping
        calculate_hydromechanic_coeffs_lin( 
                                                this->input,
                                                this->mpi_config,
                                                this->mesh_gp,
                                                this->sim_data->panels_potential,
                                                this->input->angfreqs[i],
                                                this->sim_data->added_mass,
                                                this->sim_data->damping_rad,
                                                this->sim_data->panels_pressure
                                            );

        // Calculate diffraction forces
        calculate_diffraction_forces_lin(
                                                this->input,
                                                this->mpi_config,
                                                this->mesh_gp,
                                                this->sim_data->panels_potential,
                                                this->input->angfreqs[i],
                                                this->sim_data->wave_diffrac,
                                                this->sim_data->panels_pressure
                                        );
        
        // Calculate Froude-Krylov first order forces
        calculate_froude_krylov_fo(
                                                this->input,
                                                this->mpi_config,
                                                this->mesh_gp,
                                                this->input->angfreqs[i],
                                                this->sim_data->froude_krylov
                                    );
        
        // Reduce data into the root processor variablesç
        REDUCE_FO_ROOT( added_mass,    hydmech,  cusfloat   )
        REDUCE_FO_ROOT( damping_rad,   hydmech,  cusfloat   )
        REDUCE_FO_ROOT( froude_krylov, wave_exc, cuscomplex )
        REDUCE_FO_ROOT( wave_diffrac,  wave_exc, cuscomplex )
        
        // Generate total wave excitation forces
        if ( this->mpi_config->is_root( ) )
        {
            sv_add(
                        this->sim_data->wave_exc_np,
                        this->sim_data->wave_diffrac_p0,
                        this->sim_data->froude_krylov_p0,
                        this->sim_data->wave_exc_p0
                    );
        }
        
        // Calculate RAOs
        if ( this->mpi_config->is_root( ) )
        {
            calculate_raos(
                                this->input,
                                this->sim_data->structural_mass_p0,
                                this->sim_data->added_mass_p0,
                                this->sim_data->damping_rad_p0,
                                this->sim_data->hydrostiff_p0,
                                this->sim_data->wave_diffrac_p0,
                                this->sim_data->froude_krylov_p0,
                                this->input->angfreqs[i],
                                this->sim_data->raos
                            );
        }

        
        // Storage results
        // MpiTimer storage_timer;       
        if ( this->input->out_sources )
        {
            MPI_Reduce(
                            this->sim_data->intensities,
                            this->sim_data->intensities_p0,
                            this->kernel->size( ) * ( this->input->dofs_np + this->input->heads_np ),
                            mpi_cuscomplex,
                            MPI_SUM,
                            this->mpi_config->proc_root,
                            MPI_COMM_WORLD
                        );
        }

        if ( this->input->out_potential )
        {
            MPI_Reduce(
                            this->sim_data->panels_potential,
                            this->sim_data->panels_potential_p0,
                            this->kernel->size( ) * ( this->input->dofs_np + this->input->heads_np ),
                            mpi_cuscomplex,
                            MPI_SUM,
                            this->mpi_config->proc_root,
                            MPI_COMM_WORLD
                        );
        }

        if ( this->input->out_pressure )
        {
            MPI_Reduce(
                            this->sim_data->panels_pressure,
                            this->sim_data->panels_pressure_p0,
                            this->kernel->size( ) * ( this->input->dofs_np + this->input->heads_np ),
                            mpi_cuscomplex,
                            MPI_SUM,
                            this->mpi_config->proc_root,
                            MPI_COMM_WORLD
                        );
        }
        
        if ( this->mpi_config->is_root( ) )
        {
            if ( this->input->out_hydmech )
            {
                this->output->save_hydromechanics_format(
                                                            i,
                                                            _DN_ADDED_MASS,
                                                            this->sim_data->added_mass_p0
                                                        );

                this->output->save_hydromechanics_format(
                                                            i,
                                                            _DN_DAMPING_RAD,
                                                            sim_data->damping_rad_p0
                                                        );
            }

            if ( this->input->out_diffrac )
            {
                this->output->save_wave_exciting_format(
                                                            i,
                                                            _DN_DIFFRAC,
                                                            this->sim_data->wave_diffrac_p0
                                                        );
            }

            if ( this->input->out_fk )
            {
                this->output->save_wave_exciting_format(
                                                            i,
                                                            _DN_FK,
                                                            this->sim_data->froude_krylov_p0
                                                        );
            }

            if ( this->input->out_wex )
            {
                this->output->save_wave_exciting_format(
                                                            i,
                                                            _DN_WEX,
                                                            this->sim_data->wave_exc_p0
                                                        );
            }

            if ( this->input->out_raos )
            {
                this->output->save_wave_exciting_format(
                                                            i,
                                                            _DN_RAO,
                                                            this->sim_data->raos
                                                        );
            }

            // Storage sources
            if ( this->input->out_sources )
            {
                this->output->save_fields_data( 
                                                    i,
                                                    _DN_SRC_INT,
                                                    this->sim_data->intensities_p0
                                                );
            }

            // Storage panels potential
            if ( this->input->out_potential )
            {
                this->output->save_fields_data( 
                                                    i,
                                                    _DN_POT_INT,
                                                    this->sim_data->panels_potential_p0
                                                );
            }

            // Storage panels pressure
            if ( this->input->out_pressure )
            {
                this->output->save_fields_data( 
                                                    i,
                                                    _DN_PRESS_INT,
                                                    this->sim_data->panels_pressure_p0
                                                );
            }

        }
        
        // Print out execution times
        LOG_TASK_TIME( freq_timer )
        
    }
}


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::_calculate_hydrostatics( void )
{
    LOG_TASK( "Calculating hydrostatics..." )
    MpiTimer hydro_timer;

    this->hydrostatics = new Hydrostatics*[this->input->bodies_np];
    for ( int i=0; i<input->bodies_np; i++ )
    {
        this->hydrostatics[i]   =   new Hydrostatics( 
                                                        this->input->bodies[i]->mesh,
                                                        this->input->water_density,
                                                        this->input->grav_acc,
                                                        this->input->bodies[i]->mass,
                                                        this->input->bodies[i]->cog,
                                                        this->input->bodies[i]->rad_inertia,
                                                        this->mpi_config
                                                    );
    }

    LOG_TASK_TIME( hydro_timer )

}


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::_calculate_global_static_matrixes( void )
{
    // Calculate global structural mass
    if ( this->mpi_config->is_root( ) )
    {
        calculate_global_structural_mass(
                                            input,
                                            sim_data->structural_mass_p0
                                        );
    }


    // Calculate global stiffness matrix
    if ( this->mpi_config->is_root( ) )
    {
        calculate_global_hydstiffness(
                                            input,
                                            hydrostatics,
                                            sim_data->hydrostiff_p0
                                        );
    }
}


template<std::size_t N, int mode_pf>
FrequencySolver<N, mode_pf>::FrequencySolver( Input* input_in, MpiConfig* mpi_config_in )
{
    // Storage input sytem pointer
    this->input = input_in;

    // Create MPI environment configuration
    this->mpi_config = mpi_config_in;

    // Calculate hydrostatics
    this->_calculate_hydrostatics( );

    // Create output system
    this->_initialize_output_system( );
    
    // Create mesh group
    this->_initialize_mesh_groups( );

    // Allocate space for the simulation data
    this->sim_data  = new SimulationData(
                                            this->input,
                                            this->mpi_config,
                                            this->mesh_gp->meshes_np,
                                            this->input->dofs_np,
                                            this->input->heads_np,
                                            this->mesh_gp->source_nodes_tnp
                                        );

    // Generate formulation kernel
    this->_generate_formulation_kernel( );

    // Calculate mass and hydrostatic stiffness matrixes
    // used for the calculation of the RAOs.
    this->_calculate_global_static_matrixes( );

}


template<std::size_t N, int mode_pf>
FrequencySolver<N, mode_pf>::~FrequencySolver( void )
{
    // Delete simulation data
    delete this->sim_data;

    // Delete mesh groups
    delete this->mesh_gp;
    delete this->mesh_fs_qtf_gp;

    // Delete hydrostatic objects
    for ( int i=0; i<input->bodies_np; i++ )
    {
        delete this->hydrostatics[i];
    }
    delete [] this->hydrostatics;

    // Delete output system
    delete this->output;

}


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::_generate_formulation_kernel( void )
{
    LOG_TASK( "Generating kernel..." )
    MpiTimer kernel_timer;
    this->kernel = new FormulationKernelBackend<NUM_GP, PF_OFF>( this->input, this->mpi_config, this->mesh_gp );
    LOG_TASK_TIME( kernel_timer )
}


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::_initialize_mesh_groups( void )
{
    LOG_TASK( "Initialize mesh groups..." )
    MpiTimer mesh_timer;
    
    // Group all meshes in a vector
    Mesh** all_meshes = new Mesh*[input->bodies_np];
    for ( int i=0; i<input->bodies_np; i++ )
    {
        all_meshes[i] = this->input->bodies[i]->mesh;
    }

    // Create new mesh from the meshes of all objects
    this->mesh_gp       = new MeshGroup( 
                                            all_meshes,
                                            input->bodies_np,
                                            input->is_wl_points
                                        );

    Mesh**      fs_mesh = nullptr;
    if ( this->input->is_fs_qtf )
    {
        // Create vector with the QTF FS mesh to feed the
        // MeshGroup object constructor
        fs_mesh         = new Mesh*[1];
        fs_mesh[0]      = this->input->bodies[0]->mesh_fs_qtf;

        // Create partition circle water line if required
        if ( this->input->out_qtf_so_model == 2 )
        {
            fs_mesh[0]->detect_pc_points( this->input->wl_det_prec );
        }

        // Create mesh group
        this->mesh_fs_qtf_gp    = new MeshGroup(
                                                    fs_mesh,
                                                    1,
                                                    false
                                                );
    }

    this->mesh_gp->define_mirror_panels( );

    // Delete allocated heap memory
    delete []   all_meshes;
    delete []   fs_mesh;

    LOG_TASK_TIME( mesh_timer )

}


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::_initialize_output_system( void )
{
    LOG_TASK( "Creating output system..." )
    MpiTimer output_timer;

    // Create output system
    if ( this->mpi_config->is_root( ) )
    {
        this->output = new Output( this->input );
    }

    // Storage initial parameters
    if ( this->mpi_config->is_root( ) )
    {
        // Storage frequency set
        this->output->save_frequencies( this->input->freqs );

        // Storage headings set
        this->output->save_headings( this->input->heads.data( ) );

        // Storage structural mass
        if ( this->input->out_struct_mass )
        {
            this->output->save_structural_mass( );
        }

        // Storage hydrostatic stiffness matrix
        if ( this->input->out_hydstiff )
        {
            this->output->save_hydstiffness( this->hydrostatics );
        }

        // Storage mesh
        if ( this->input->out_mesh )
        {
            this->output->save_mesh( );
        }

    }

    LOG_TASK_TIME( output_timer )

}