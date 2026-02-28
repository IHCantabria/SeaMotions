
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
#include "../inout/hdf5_time_series_exporter.hpp"


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
    std::string                 field_name  = "";
    cut::CusTensor<cusfloat>*   field_data  = nullptr;

    const std::size_t heads_np = static_cast<std::size_t>( this->_input->heads_np );
    for ( std::size_t i=0; i<heads_np; i++ )
    {
        if ( this->_field_points_def->out_potential )
        {
            // Add total potential field dataset
            field_name  = compose_field_name( "potential", this->_input->heads[i] );
            field_data  = this->template get_field_data_scalar< FieldTypeE::POTENTIAL, FieldComponentE::TOTAL >( i );
            if ( this->_mpi_config->is_root( ) )
                this->_exporter->append_step( field_name, field_data );

            if ( this->_field_points_def->out_components )
            {
                // Add potential field components dataset
                field_name = compose_field_name( "potential_incident", this->_input->heads[i] );
                field_data = this->template get_field_data_scalar< FieldTypeE::POTENTIAL, FieldComponentE::INCIDENT >( i );
                if ( this->_mpi_config->is_root( ) )
                    this->_exporter->append_step( field_name, field_data );

                field_name = compose_field_name( "potential_scatter", this->_input->heads[i] );
                field_data = this->template get_field_data_scalar< FieldTypeE::POTENTIAL, FieldComponentE::SCATTERED >( i );
                if ( this->_mpi_config->is_root( ) )
                    this->_exporter->append_step( field_name, field_data );

                field_name = compose_field_name( "potential_radiation", this->_input->heads[i] );
                field_data = this->template get_field_data_scalar< FieldTypeE::POTENTIAL, FieldComponentE::RADIATED >( i );
                if ( this->_mpi_config->is_root( ) )
                    this->_exporter->append_step( field_name, field_data );
                
            }
        }

        if ( this->_field_points_def->out_pressure )
        {
            // Add total pressure field dataset
            field_name  = compose_field_name( "pressure", this->_input->heads[i] );
            field_data  = this->template get_field_data_scalar< FieldTypeE::PRESSURE, FieldComponentE::TOTAL >( i );
            if ( this->_mpi_config->is_root( ) )
                this->_exporter->append_step( field_name, field_data );

            if ( this->_field_points_def->out_components )
            {
                // Add pressure field components dataset
                field_name = compose_field_name( "pressure_incident", this->_input->heads[i] );
                field_data = this->template get_field_data_scalar< FieldTypeE::PRESSURE, FieldComponentE::INCIDENT >( i );
                if ( this->_mpi_config->is_root( ) )
                    this->_exporter->append_step( field_name, field_data );

                field_name = compose_field_name( "pressure_scatter", this->_input->heads[i] );
                field_data = this->template get_field_data_scalar< FieldTypeE::PRESSURE, FieldComponentE::SCATTERED >( i );
                if ( this->_mpi_config->is_root( ) )
                    this->_exporter->append_step( field_name, field_data );

                field_name = compose_field_name( "pressure_radiation", this->_input->heads[i] );
                field_data = this->template get_field_data_scalar< FieldTypeE::PRESSURE, FieldComponentE::RADIATED >( i );
                if ( this->_mpi_config->is_root( ) )
                    this->_exporter->append_step( field_name, field_data );
                
            }
        }

        if ( this->_field_points_def->out_velocity )
        {
            // Add total pressure field dataset
            field_name  = compose_field_name( "velocity", this->_input->heads[i] );
            field_data  = this->template get_field_data_vector< FieldTypeE::VELOCITY, FieldComponentE::TOTAL >( i );
            if ( this->_mpi_config->is_root( ) )
                this->_exporter->append_step( field_name, field_data );

            if ( this->_field_points_def->out_components )
            {
                // Add pressure field components dataset
                field_name = compose_field_name( "velocity_incident", this->_input->heads[i] );
                field_data = this->template get_field_data_vector< FieldTypeE::VELOCITY, FieldComponentE::INCIDENT >( i );
                if ( this->_mpi_config->is_root( ) )
                    this->_exporter->append_step( field_name, field_data );

                field_name = compose_field_name( "velocity_scatter", this->_input->heads[i] );
                field_data = this->template get_field_data_vector< FieldTypeE::VELOCITY, FieldComponentE::SCATTERED >( i );
                if ( this->_mpi_config->is_root( ) )
                    this->_exporter->append_step( field_name, field_data );

                field_name = compose_field_name( "velocity_radiation", this->_input->heads[i] );
                field_data = this->template get_field_data_vector< FieldTypeE::VELOCITY, FieldComponentE::RADIATED >( i );
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

    const std::size_t heads_np = static_cast<std::size_t>( this->_input->heads_np );
    for ( std::size_t i=0; i<heads_np; i++ )
    {
        if ( this->_field_points_def->out_potential )
        {
            // Add total potential field dataset
            this->_exporter->add_field( compose_field_name( "potential", this->_input->heads[i] ), hdnp, 1 );

            if ( this->_field_points_def->out_components )
            {
                // Add potential field components dataset
                this->_exporter->add_field( compose_field_name( "potential_incident", this->_input->heads[i] ),  hdnp, 1 );
                this->_exporter->add_field( compose_field_name( "potential_scatter", this->_input->heads[i] ),   hdnp, 1 );
                this->_exporter->add_field( compose_field_name( "potential_radiation", this->_input->heads[i] ), hdnp, 1 );
            }
        }

        if ( this->_field_points_def->out_pressure )
        {
            // Add total potential field dataset
            this->_exporter->add_field( compose_field_name( "pressure", this->_input->heads[i] ), hdnp, 1 );

            if ( this->_field_points_def->out_components )
            {
                // Add potential field components dataset
                this->_exporter->add_field( compose_field_name( "pressure_incident", this->_input->heads[i] ),  hdnp, 1 );
                this->_exporter->add_field( compose_field_name( "pressure_scatter", this->_input->heads[i] ),   hdnp, 1 );
                this->_exporter->add_field( compose_field_name( "pressure_radiation", this->_input->heads[i] ), hdnp, 1 );
            }
        }

        if ( this->_field_points_def->out_velocity )
        {
            // Add total potential field dataset
            this->_exporter->add_field( compose_field_name( "velocity", this->_input->heads[i] ), hdnp, 3 );

            if ( this->_field_points_def->out_components )
            {
                // Add potential field components dataset
                this->_exporter->add_field( compose_field_name( "velocity_incident", this->_input->heads[i] ),  hdnp, 3 );
                this->_exporter->add_field( compose_field_name( "velocity_scatter", this->_input->heads[i] ),   hdnp, 3 );
                this->_exporter->add_field( compose_field_name( "velocity_radiation", this->_input->heads[i] ), hdnp, 3 );
            }
        }
    }
}


template<typename T, typename ModeComp>
template<FieldTypeE field_type, FieldComponentE field_comp>
cut::CusTensor<cusfloat>* FieldMeshData<T, ModeComp>::get_field_data_scalar( std::size_t heading_index )
{
    // Get field data from radiation and diffraction data container
    const cut::CusTensor<cusfloat>* field_temp = this->_rdd->template get_field_data<field_type, field_comp>( heading_index );

    // MPI Gather field data to root process
    cusfloat* recv_buf = this->_mpi_config->is_root( ) ? this->_data_scalar_r.data( ) : nullptr;
    MPI_Gather( 
                    field_temp->data( ), 
                    this->_rdd->get_size_local_fp( ), 
                    mpi_cusfloat, 
                    recv_buf, 
                    this->_rdd->get_size_local_fp( ), 
                    mpi_cusfloat, 
                    0, 
                    this->_mpi_config->mpi_comm 
                );

    MPI_Barrier( this->_mpi_config->mpi_comm );

    return &(this->_data_scalar_r);
}


template<typename T, typename ModeComp>
template<FieldTypeE field_type, FieldComponentE field_comp>
cut::CusTensor<cusfloat>* FieldMeshData<T, ModeComp>::get_field_data_vector( std::size_t heading_index )
{
    // Get field data from radiation and diffraction data container
    std::size_t ndim = 0;
    if constexpr( field_type == FieldTypeE::VELOCITY )
    {
        // Set field dimenstions for velocity vector field
        ndim = 3;
        
        // Get velocity X field data
        const cut::CusTensor<cusfloat>* vel_x = this->_rdd->template get_field_data<FieldTypeE::VELOCITY_X, field_comp>( heading_index );
        for ( std::size_t i=0; i<this->_rdd->get_size_local_fp( ); i++ )
        {
            this->_data_vector[ 3*i ]     = vel_x->data( )[i];
        }

        // Get velocity Y field data
        const cut::CusTensor<cusfloat>* vel_y = this->_rdd->template get_field_data<FieldTypeE::VELOCITY_Y, field_comp>( heading_index );
        for ( std::size_t i=0; i<this->_rdd->get_size_local_fp( ); i++ )
        {
            this->_data_vector[ 3*i + 1 ] = vel_y->data( )[i];
        }

        // Get velocity Z field data
        const cut::CusTensor<cusfloat>* vel_z = this->_rdd->template get_field_data<FieldTypeE::VELOCITY_Z, field_comp>( heading_index );
        for ( std::size_t i=0; i<this->_rdd->get_size_local_fp( ); i++ )
        {
            this->_data_vector[ 3*i + 2 ] = vel_z->data( )[i];
        }

    }

    // MPI Gather field data to root process
    cusfloat* recv_buf = this->_mpi_config->is_root( ) ? this->_data_vector_r.data( ) : nullptr;
    MPI_Gather( 
                    this->_data_vector.data( ), 
                    this->_rdd->get_size_local_fp( ) * ndim, 
                    mpi_cusfloat, 
                    recv_buf, 
                    this->_rdd->get_size_local_fp( ) * ndim, 
                    mpi_cusfloat, 
                    0, 
                    this->_mpi_config->mpi_comm 
                );

    MPI_Barrier( this->_mpi_config->mpi_comm );

    return &(this->_data_vector_r);
}


template<typename T, typename ModeComp>
RadDiffData<T, ModeComp>* FieldMeshData<T, ModeComp>::get_rdd( void ) const
{
    return this->_rdd;
}


template<typename T, typename ModeComp>
FieldMeshData<T, ModeComp>::FieldMeshData( 
                                            Input*          input,
                                            FieldPointsDef* field_points_def,
                                            std::string     out_folder_path,
                                            MpiConfig*      mpi_config 
                                        )
{
    // Get a local copy of the configuration
    this->_input            = input;
    this->_field_points_def = field_points_def;
    this->_mpi_config       = mpi_config;

    // Compose mesh file path
    namespace fs            = std::filesystem;
    fs::path mesh_fopath    = fs::path( this->_input->folder_path ) / fs::path( std::string( MESH_FOLDER_NAME ) );
    fs::path mesh_fipath    = mesh_fopath / fs::path( field_points_def->mesh_finame );

    // Create mesh pointer
    this->_mesh             = new Mesh( 
                                            mesh_fipath.string( ), 
                                            field_points_def->mesh_body_name,
                                            PanelTypeE::FIELD_POINT
                                        );
    this->_is_heap_mesh     = true;

    // Create radiation and diffraction data pointer
    this->_rdd              = new RadDiffData<T, ModeComp>( 
                                                                mpi_config,
                                                                this->_mesh,
                                                                input->angfreqs_np,
                                                                input->heads_np,
                                                                input->dofs_np,
                                                                true
                                                            );
    this->_is_heap_rdd      = true;

    // Allocate memory for root data if parallel output is enabled
    this->_data_vector.resize( this->_rdd->get_size_local_fp( ) * 3 );
    if ( mpi_config->is_root( ) )
    {
        this->_data_scalar_r.resize( this->_mesh->nodes_np );
        this->_data_vector_r.resize( this->_mesh->nodes_np * 3 );
    }

    // Generate HDF5 exporter interface
    if ( mpi_config->is_root( ) )
    {
        this->_exporter         = new HDF5TimeSeriesExporter( 
                                                                out_folder_path, 
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

    // Delete HDF5 exporter interface if allocated on heap
    if ( this->_is_heap_exporter && this->_exporter != nullptr )
    {
        delete this->_exporter;
        this->_exporter = nullptr;
    }

}