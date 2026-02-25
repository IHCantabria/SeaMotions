
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
#include <sstream>

// Include local modules
#include "field_mesh_data.hpp"


inline std::string compose_field_name( const std::string& field_base_name, const cusfloat heading )
{
    std::stringstream ss;
    ss << std::fixed << std::setprecision( 2 ) << heading;
    return field_base_name + "_" + ss.str( );
}


template<typename T, typename ModeComp>
void FieldMeshData<T, ModeComp>::add_step( cusfloat freq )
{
    // Increment internal step counter
    this->_step++;
    this->_exporter->append_time( freq );

    // Retrive field data for the current step and append it to the exporter interface
    std::string field_name  = "";
    cusfloat*   field_data  = nullptr;

    for ( std::size_t i=0; i<this->_rdd->headings_np; i++ )
    {
        if ( this->_config.out_potential )
        {
            // Add total potential field dataset
            field_name  = compose_field_name( "potential", this->_rdd->headings[i] );
            field_data  = this->get_field_data_scalar< FieldTypeE::POTENTIAL, FieldComponentE::TOTAL >( i );
            if ( this->_mpi_config->is_root( ) )
                this->_exporter->append_step( field_name, field_data );

            if ( this->_config.out_components )
            {
                // Add potential field components dataset
                field_name = compose_field_name( "potential_incident", this->_rdd->headings[i] );
                field_data = this->get_field_data_scalar< FieldTypeE::POTENTIAL, FieldComponentE::INCIDENT >( i );
                if ( this->_mpi_config->is_root( ) )
                    this->_exporter->append_step( field_name, field_data );

                field_name = compose_field_name( "potential_scatter", this->_rdd->headings[i] );
                field_data = this->get_field_data_scalar< FieldTypeE::POTENTIAL, FieldComponentE::SCATTERED >( i );
                if ( this->_mpi_config->is_root( ) )
                    this->_exporter->append_step( field_name, field_data );

                field_name = compose_field_name( "potential_radiation", this->_rdd->headings[i] );
                field_data = this->get_field_data_scalar< FieldTypeE::POTENTIAL, FieldComponentE::RADIATED >( i );
                if ( this->_mpi_config->is_root( ) )
                    this->_exporter->append_step( field_name, field_data );
                
            }
        }

        if ( this->_config.out_pressure )
        {
            // Add total pressure field dataset
            field_name  = compose_field_name( "pressure", this->_rdd->headings[i] );
            field_data  = this->get_field_data_scalar< FieldTypeE::PRESSURE, FieldComponentE::TOTAL >( i );
            if ( this->_mpi_config->is_root( ) )
                this->_exporter->append_step( field_name, field_data );

            if ( this->_config.out_components )
            {
                // Add pressure field components dataset
                field_name = compose_field_name( "pressure_incident", this->_rdd->headings[i] );
                field_data = this->get_field_data_scalar< FieldTypeE::PRESSURE, FieldComponentE::INCIDENT >( i );
                if ( this->_mpi_config->is_root( ) )
                    this->_exporter->append_step( field_name, field_data );

                field_name = compose_field_name( "pressure_scatter", this->_rdd->headings[i] );
                field_data = this->get_field_data_scalar< FieldTypeE::PRESSURE, FieldComponentE::SCATTERED >( i );
                if ( this->_mpi_config->is_root( ) )
                    this->_exporter->append_step( field_name, field_data );

                field_name = compose_field_name( "pressure_radiation", this->_rdd->headings[i] );
                field_data = this->get_field_data_scalar< FieldTypeE::PRESSURE, FieldComponentE::RADIATED >( i );
                if ( this->_mpi_config->is_root( ) )
                    this->_exporter->append_step( field_name, field_data );
                
            }
        }

        if ( this->_config.out_velocity )
        {
            // Add total pressure field dataset
            field_name  = compose_field_name( "velocity", this->_rdd->headings[i] );
            field_data  = this->get_field_data_vector< FieldTypeE::VELOCITY, FieldComponentE::TOTAL >( i );
            if ( this->_mpi_config->is_root( ) )
                this->_exporter->append_step( field_name, field_data );

            if ( this->_config.out_components )
            {
                // Add pressure field components dataset
                field_name = compose_field_name( "velocity_incident", this->_rdd->headings[i] );
                field_data = this->get_field_data_vector< FieldTypeE::VELOCITY, FieldComponentE::INCIDENT >( i );
                if ( this->_mpi_config->is_root( ) )
                    this->_exporter->append_step( field_name, field_data );

                field_name = compose_field_name( "velocity_scatter", this->_rdd->headings[i] );
                field_data = this->get_field_data_vector< FieldTypeE::VELOCITY, FieldComponentE::SCATTERED >( i );
                if ( this->_mpi_config->is_root( ) )
                    this->_exporter->append_step( field_name, field_data );

                field_name = compose_field_name( "velocity_radiation", this->_rdd->headings[i] );
                field_data = this->get_field_data_vector< FieldTypeE::VELOCITY, FieldComponentE::RADIATED >( i );
                if ( this->_mpi_config->is_root( ) )
                    this->_exporter->append_step( field_name, field_data );
                
            }
        }
    }

    // Flush XDMF metadata after appending the current step data to ensure it is updated in case of an unexpected interruption
    if ( this->_mpi_config->is_root( ) )
    {
        this->_exporter->write_xdmf( );
    }
}

template<typename T, typename ModeComp>
void FieldMeshData<T, ModeComp>::_configure_exporter( void )
{
    std::size_t hdnp = 1;

    for ( std::size_t i=0; i<this->_rdd->headings_np; i++ )
    {
        if ( this->_config.out_potential )
        {
            // Add total potential field dataset
            this->_exporter->add_field( compose_field_name( "potential", this->_rdd->headings[i] ), hdnp, 1 );

            if ( this->_config.out_components )
            {
                // Add potential field components dataset
                this->_exporter->add_field( compose_field_name( "potential_incident", this->_rdd->headings[i] ),  hdnp, 1 );
                this->_exporter->add_field( compose_field_name( "potential_scatter", this->_rdd->headings[i] ),   hdnp, 1 );
                this->_exporter->add_field( compose_field_name( "potential_radiation", this->_rdd->headings[i] ), hdnp, 1 );
            }
        }

        if ( this->_config.out_potential )
        {
            // Add total potential field dataset
            this->_exporter->add_field( compose_field_name( "pressure", this->_rdd->headings[i] ), hdnp, 1 );

            if ( this->_config.out_components )
            {
                // Add potential field components dataset
                this->_exporter->add_field( compose_field_name( "pressure_incident", this->_rdd->headings[i] ),  hdnp, 1 );
                this->_exporter->add_field( compose_field_name( "pressure_scatter", this->_rdd->headings[i] ),   hdnp, 1 );
                this->_exporter->add_field( compose_field_name( "pressure_radiation", this->_rdd->headings[i] ), hdnp, 1 );
            }
        }

        if ( this->_config.out_velocity )
        {
            // Add total potential field dataset
            this->_exporter->add_field( compose_field_name( "velocity", this->_rdd->headings[i] ), hdnp, 3 );

            if ( this->_config.out_components )
            {
                // Add potential field components dataset
                this->_exporter->add_field( compose_field_name( "velocity_incident", this->_rdd->headings[i] ),  hdnp, 3 );
                this->_exporter->add_field( compose_field_name( "velocity_scatter", this->_rdd->headings[i] ),   hdnp, 3 );
                this->_exporter->add_field( compose_field_name( "velocity_radiation", this->_rdd->headings[i] ), hdnp, 3 );
            }
        }
    }
}


template<typename T, typename ModeComp>
template<FieldTypeE field_type, FieldComponentE field_comp>
const cusfloat* FieldMeshData<T, ModeComp>::get_field_data_scalar( std::size_t heading_index ) const
{
    // Get field data from radiation and diffraction data container
    cusfloat* field_temp = this->_rdd->get_field_data<field_type, field_comp>( heading_index );

    // MPI Gather field data to root process
    MPI_Gather( 
                    field_temp, 
                    this->_rdd->get_size_local_fp( ), 
                    mpi_cusfloat, 
                    this->_data_scalar_r, 
                    this->_rdd->get_size_local_fp( ), 
                    mpi_cusfloat, 
                    0, 
                    this->_mpi_config->mpi_comm 
                );

    MPI_Barrier( this->_mpi_config->mpi_comm );

    return this->_data_scalar_r;
}


template<typename T, typename ModeComp>
template<FieldTypeE field_type, FieldComponentE field_comp>
const cusfloat* FieldMeshData<T, ModeComp>::get_field_data_vector( std::size_t heading_index ) const
{
    // Get field data from radiation and diffraction data container
    std::size_t ndim = 0;
    if constexpr( field_type == FieldTypeE::VELOCITY )
    {
        ndim = 3;
        cusfloat* vel_x = this->_rdd->get_field_data<FieldTypeE::VELOCITY_X, field_comp>( heading_index );
        cusfloat* vel_y = this->_rdd->get_field_data<FieldTypeE::VELOCITY_Y, field_comp>( heading_index );
        cusfloat* vel_z = this->_rdd->get_field_data<FieldTypeE::VELOCITY_Z, field_comp>( heading_index );

        for ( std::size_t i=0; i<this->_rdd->get_size_local_fp( ); i++ )
        {
            this->_data_vector[ 3*i ]     = vel_x[i];
            this->_data_vector[ 3*i + 1 ] = vel_y[i];
            this->_data_vector[ 3*i + 2 ] = vel_z[i];
        }

    }

    // MPI Gather field data to root process
    MPI_Gather( 
                    this->_data_vector, 
                    this->_rdd->get_size_local_fp( ) * ndim, 
                    mpi_cusfloat, 
                    this->_data_vector_r, 
                    this->_rdd->get_size_local_fp( ) * ndim, 
                    mpi_cusfloat, 
                    0, 
                    this->_mpi_config->mpi_comm 
                );

    MPI_Barrier( this->_mpi_config->mpi_comm );

    return this->_data_vector_r;
}


template<typename T, typename ModeComp>
FieldMeshData<T, ModeComp>::FieldMeshData( 
                                            FieldMeshDataConfig&    config,
                                            MpiConfig*              mpi_config
                                        )
{
    // Get a local copy of the configuration
    this->_config           = config;
    this->_mpi_config       = mpi_config;

    // Create mesh pointer
    this->_mesh             = new Mesh( config.mesh_file_path, config.body_name );
    this->_is_heap_mesh     = true;

    // Create radiation and diffraction data pointer
    this->_rdd              = new RadDiffData<T, ModeComp>( 
                                                                mpi_config,
                                                                this->_mesh,
                                                                config.freqs_np,
                                                                config.headings_np,
                                                                config.dofs_np,
                                                                true
                                                            );
    this->_is_heap_rdd      = true;

    // Allocate memory for root data if parallel output is enabled
    this->_data_vector  = generate_empty_vector<cusfloat>( this->_rdd->get_size_local_fp( ) * 3 );
    if ( mpi_config->is_root( ) )
    {
        this->_data_scalar_r = generate_empty_vector<cusfloat>( this->_mesh->nodes_np );
        this->_data_vector_r = generate_empty_vector<cusfloat>( this->_mesh->nodes_np * 3 );
    }
    this->_is_data_heap = true;

    // Generate HDF5 exporter interface
    if ( mpi_config->is_root( ) )
    {
        this->_exporter         = new HDF5TimeSeriesExporter( 
                                                                this->_config.out_folder_path, 
                                                                this->_mesh 
                                                            );
        this->_is_heap_exporter = true;

    }

    // Configure exporter output settings
    this->_configure_exporter( );
    
}


template<typename T, typename ModeComp>
FieldMeshData<T, ModeComp>::~FieldMeshData( void )
{
    // Delete mesh pointer if allocated on heap
    if ( this->_is_heap_mesh && this->_mesh != nullptr )
    {
        delete this->_mesh;
        this->_mesh = nullptr;
    }

    // Delete radiation and diffraction data pointer if allocated on heap
    if ( this->_is_heap_rdd && this->_rdd != nullptr )
    {
        delete this->_rdd;
        this->_rdd = nullptr;
    }

    // Delete root scalar data array if allocated on heap
    if ( this->_mpi_config->is_root( ) && this->_is_data_heap && this->_data_scalar_r != nullptr )
    {
        mkl_free( this->_data_scalar_r );
        this->_data_scalar_r = nullptr;
    }

    // Delete root vector data array if allocated on heap
    mkl_free( this->_data_vector );
    this->_data_vector = nullptr;
    if ( this->_mpi_config->is_root( ) && this->_is_data_heap && this->_data_vector_r != nullptr )
    {
        mkl_free( this->_data_vector_r );
        this->_data_vector_r = nullptr;
    }

    // Delete HDF5 exporter interface if allocated on heap
    if ( this->_is_heap_exporter && this->_exporter != nullptr )
    {
        delete this->_exporter;
        this->_exporter = nullptr;
    }

}