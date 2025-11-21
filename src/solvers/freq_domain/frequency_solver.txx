
// Include local module
#include "../../containers/mpi_timer.hpp"
#include "diffraction.hpp"
#include "frequency_solver.hpp"
#include "froude_krylov.hpp"
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
        INFO( "FO -> Calculating period: " << angfreq_to_period( this->input->angfreqs[i] ) << " s\n" )

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
        if ( mpi_config->is_root( ) )
        {
            sv_add(
                        this->sim_data->wave_exc_np,
                        this->sim_data->wave_diffrac_p0,
                        this->sim_data->froude_krylov_p0,
                        this->sim_data->wave_exc_p0
                    );
        }
        
        // Calculate RAOs
        if ( mpi_config->is_root( ) )
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
        MpiTimer storage_timer;       
        if ( input->out_pressure )
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
                                                    this->sim_data->intensities
                                                );
            }

            // Storage panels potential
            if ( this->input->out_potential )
            {
                this->output->save_fields_data( 
                                                    i,
                                                    _DN_POT_INT,
                                                    this->sim_data->panels_potential
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
        storage_timer.stop( );
        freq_timer.stop( );

        if ( this->mpi_config->is_root( ) )
        {
            std::cout << "Execution time [s]: " << freq_timer;
            std::cout << " - Disk Storage time [s]: " << storage_timer << std::endl;
        }
    }
}


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::calculate_hydrostatics( void )
{
    INFO( "Calculating hydrostatics..." )
    MpiTimer hydro_timer;

    this->hydrostatics = new Hydrostatics*[this->input->bodies_np];
    for ( int i=0; i<input->bodies_np; i++ )
    {
        hydrostatics[i] =   new Hydrostatics( 
                                                this->input->bodies[i]->mesh,
                                                this->input->water_density,
                                                this->input->grav_acc,
                                                this->input->bodies[i]->mass,
                                                this->input->bodies[i]->cog,
                                                this->input->bodies[i]->rad_inertia,
                                                this->mpi_config
                                            );
    }

    hydro_timer.stop( );
    INFO( " --> Done! Elapsed Time: " << hydro_timer << "\n" )
}


template<std::size_t N, int mode_pf>
FrequencySolver<N, mode_pf>::FrequencySolver( Input* input_in )
{
    // Storage input sytem pointer
    this->input = input_in;

    // Create MPI environment configuration
    this->mpi_config = new MpiConfig( );

    // Calculate hydrostatics
    this->calculate_hydrostatics( );

    // Create output system
    this->initialize_output_system( );
    
    // Create mesh group
    this->initialize_mesh_groups( );


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
    this->kernel = new FormulationKernelBackend<NUM_GP, PF_OFF>( this->input, this->mpi_config, this->mesh_gp );

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

    // Delete MPI config
    delete this->mpi_config;

}


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::initialize_mesh_groups( void )
{
    INFO( "Initialize mesh groups..." )
    MpiTimer mesh_timer;
    
    // Group all meshes in a vector
    Mesh** all_meshes = new Mesh*[input->bodies_np];
    for ( int i=0; i<input->bodies_np; i++ )
    {
        all_meshes[i] = this->input->bodies[i]->mesh;
    }

    // Create new mesh from the meshes of all objects
    this->mesh_gp   = new MeshGroup( 
                                        all_meshes,
                                        input->bodies_np,
                                        input->is_wl_points
                                    );

    Mesh**      fs_mesh         = nullptr;
    if ( this->input->is_fs_qtf )
    {
        // Create vector with the QTF FS mesh to feed the
        // MeshGroup object constructor
        fs_mesh         = new Mesh*[1];
        fs_mesh[0]      = this->input->bodies[0]->mesh_fs_qtf;

        // Create partition circle water line if required
        if ( input->out_qtf_so_model == 2 )
        {
            fs_mesh[0]->detect_pc_points( input->wl_det_prec );
        }

        // Create mesh group
        this->mesh_fs_qtf_gp    = new MeshGroup(
                                                    fs_mesh,
                                                    1,
                                                    false
                                                );
    }

    if ( input->is_log_sin_ana )
    {
        mesh_gp->define_mirror_panels( );
    }

    // Delete allocated heap memory
    delete []   all_meshes;
    delete []   fs_mesh;

    mesh_timer.stop( );
    INFO( " --> Done! Time Elapsed: " << mesh_timer << "\n" )

}


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::initialize_output_system( void )
{
    INFO( "Creating output system..." )
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

    INFO( " --> Done!" << "\n" ) // Output system created and initialized
}