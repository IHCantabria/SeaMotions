
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


#include "hdf5_time_series_exporter_mpi.hpp"


HDF5TimeSeriesExporterMPI::HDF5TimeSeriesExporterMPI(
                                                        const   std::string&    folder_path,
                                                                Mesh*           mesh,
                                                                MpiConfig*      mpi_config
                                                    )
    : HDF5TimeSeriesExporterBase( folder_path, mesh ),
      _mpi_config( mpi_config )
{
    int start_pos = 0;
    int end_pos   = 0;

    this->_mpi_config->get_1d_bounds(
                                        mesh->nodes_np,
                                        start_pos,
                                        end_pos
                                    );

    this->_global_np   = static_cast<std::size_t>( mesh->nodes_np );
    this->_local_np    = static_cast<std::size_t>( end_pos - start_pos );
    this->_node_offset = static_cast<std::size_t>( start_pos );

    this->write_mesh();
    this->initialize_time_series();
}


void HDF5TimeSeriesExporterMPI::add_field(
                                            std::string             field_name,
                                            hsize_t                 samples_multiplier,
                                            hsize_t                 components_np
                                        )
{
    const auto write_dims = this->register_field_metadata( field_name, samples_multiplier, components_np );

    // Prepare collective file access via MPI-IO
    hid_t plist = H5Pcreate( H5P_FILE_ACCESS );
    H5Pset_fapl_mpio(
                        plist, 
                        this->_mpi_config->mpi_comm, 
                        MPI_INFO_NULL
                    );

    // Open shared field file without truncating previous content
    hid_t file  = H5Fopen(
                                this->_field_file.c_str( ), 
                                H5F_ACC_RDWR,
                                plist
                            );

    // Property list no longer needed
    H5Pclose( plist );

    // Ensure the container group exists (safe for all ranks)
    ensureGroup( file, "/Fields" );

    // Dataset shapes (time dimension starts empty/unlimited)
    std::vector<hsize_t> initial_dims   = { 0, write_dims[1], write_dims[2] };
    std::vector<hsize_t> max_dims       = { H5S_UNLIMITED, write_dims[1], write_dims[2] };
    std::vector<hsize_t> chunk_dims     = write_dims;

    // Create extendable dataset collectively
    create_hdf5_extendable_dataset_mpi(
                                            file,
                                            std::string( "/Fields/" ) + field_name,
                                            cusfloat_h5,
                                            initial_dims,
                                            max_dims,
                                            chunk_dims
                                        );

    // Close file before proceeding
    H5Fclose( file );

    // Keep all ranks in sync before continuing
    MPI_Barrier( this->_mpi_config->mpi_comm );
}


void HDF5TimeSeriesExporterMPI::write_mesh( void )
{
    if ( this->_mpi_config->is_root() )
    {
        this->write_mesh_payload();
    }

    MPI_Barrier( this->_mpi_config->mpi_comm );
}


void HDF5TimeSeriesExporterMPI::initialize_time_series()
{
    hid_t plist = H5Pcreate( H5P_FILE_ACCESS );
    H5Pset_fapl_mpio(
                        plist,
                        this->_mpi_config->mpi_comm,
                        MPI_INFO_NULL
                    );

    hid_t file = H5Fcreate(
                                this->_field_file.c_str(),
                                H5F_ACC_TRUNC,
                                H5P_DEFAULT,
                                plist
                            );
    H5Pclose( plist );

    ensureGroup( file, "/Fields" );
    H5Fclose( file );

    MPI_Barrier( this->_mpi_config->mpi_comm );
}


void HDF5TimeSeriesExporterMPI::append_time( cusfloat append_value )
{
    if ( this->_mpi_config->is_root() )
    {
        this->record_time_value( append_value );
    }

    MPI_Barrier( this->_mpi_config->mpi_comm );
}


void HDF5TimeSeriesExporterMPI::append_step(
                                                std::string                 field_name,
                                                cut::CusTensor<cusfloat>*   field_data
                                            )
{
    hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist, this->_mpi_config->mpi_comm, MPI_INFO_NULL);

    hid_t file = H5Fopen(this->_field_file.c_str(), H5F_ACC_RDWR, plist);
    H5Pclose(plist);

    append_to_hdf5_dataset_mpi(
                                    file, 
                                    "/Fields/" + field_name, 
                                    field_data->data( ),
                                    this->_local_np,
                                    this->_node_offset,
                                    this->_field_dims[field_name][2]
                                );

    H5Fclose(file);

    MPI_Barrier( this->_mpi_config->mpi_comm );
}


void HDF5TimeSeriesExporterMPI::write_xdmf( void )
{
    if ( this->_mpi_config->is_root() )
    {
        this->write_xdmf_payload();
    }

    MPI_Barrier( this->_mpi_config->mpi_comm );
}


void HDF5TimeSeriesExporterMPI::ensureGroup(hid_t file, const std::string& path)
{
    if (H5Oexists_by_name(file, path.c_str(), H5P_DEFAULT) > 0)
    {
        return;
    }

    hid_t gid = H5Gcreate(file, path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(gid);
}
