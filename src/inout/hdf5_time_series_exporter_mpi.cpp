
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

#include <fstream>
#include <stdexcept>


HDF5TimeSeriesExporterMPI::HDF5TimeSeriesExporterMPI(
                                                        const   std::string&    folder_path,
                                                                Mesh*           mesh,
                                                                MpiConfig*      mpi_config
                                                    )
{
    // Store mesh info
    this->_folder_path  = folder_path;
    this->_mesh         = mesh;
    this->_mpi_config   = mpi_config;

    // Compose file paths
    this->_mesh_file    = folder_path + "/mesh.h5";
    this->_field_file   = folder_path + "/fields.h5";
    this->_xdmf_file    = folder_path + "/data.xdmf2";

    // Define global and local node counts
    int start_pos = 0;
    int end_pos   = 0;
    
    this->_mpi_config->get_1d_bounds( 
                                        mesh->nodes_np,
                                        start_pos,
                                        end_pos
                                    );  

    this->_global_np    = static_cast<std::size_t>( mesh->nodes_np );
    this->_local_np     = static_cast<std::size_t>( end_pos - start_pos );
    this->_node_offset  = static_cast<std::size_t>( start_pos );

    // Write mesh on initialization (static mesh)
    this->write_mesh( );

    // Initialize time series datasets (extendable along time dimension)
    this->initialize_time_series( );
}


void HDF5TimeSeriesExporterMPI::add_field(
                                            std::string             field_name,
                                            hsize_t                 comps_np
                                        )
{
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
    std::vector<hsize_t> initial_dims   = {0, this->_global_np, comps_np};
    std::vector<hsize_t> max_dims       = {H5S_UNLIMITED, this->_global_np, comps_np};
    std::vector<hsize_t> chunk_dims     = { 1, this->_global_np, comps_np };

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

    // Record field metadata for later hyperslab writes
    std::vector<hsize_t> write_dims = {1, this->_global_np, comps_np};
    this->_field_names.push_back( field_name );
    this->_field_dims[field_name] = write_dims;

    // Keep all ranks in sync before continuing
    MPI_Barrier( this->_mpi_config->mpi_comm );
}


void HDF5TimeSeriesExporterMPI::write_mesh( void )
{
    if ( this->_mpi_config->is_root( ) )
    {
        hid_t file = H5Fcreate(
                                    this->_mesh_file.c_str( ), 
                                    H5F_ACC_TRUNC,
                                    H5P_DEFAULT, 
                                    H5P_DEFAULT
                                );

        hid_t meshGroup = H5Gcreate(
                                        file, 
                                        "/Mesh",
                                        H5P_DEFAULT,
                                        H5P_DEFAULT,
                                        H5P_DEFAULT
                                    );

        std::size_t nodes_np = static_cast<std::size_t>( this->_mesh->nodes_np );
        cusfloat* node_view_data = generate_empty_vector<cusfloat>( nodes_np * 3 );
        for ( std::size_t i=0; i<nodes_np; i++ )
        {
            node_view_data[3*i + 0] = this->_mesh->x[i];
            node_view_data[3*i + 1] = this->_mesh->y[i];
            node_view_data[3*i + 2] = this->_mesh->z[i];
        }

        write_hdf5_dataset_compressed(
                                        file,
                                        "/Mesh/Points",
                                        cusfloat_h5,
                                        {static_cast<hsize_t>(nodes_np), 3},
                                        node_view_data,
                                        0
                                    );

        mkl_free( node_view_data );

        std::vector<int> mixed;
        int ernl = this->_mesh->enrl;
        for ( size_t i = 0; i < static_cast<size_t>(this->_mesh->elems_np); i++ )
        {
            if ( this->_mesh->elems_type[i] == 5 )
            {
                mixed.push_back(4);
                mixed.push_back(this->_mesh->elems[ernl*i + 1]);
                mixed.push_back(this->_mesh->elems[ernl*i + 2]);
                mixed.push_back(this->_mesh->elems[ernl*i + 3]);
                this->_tri_np++;
            }
            else if ( this->_mesh->elems_type[i] == 9 )
            {
                mixed.push_back(5);
                mixed.push_back(this->_mesh->elems[ernl*i + 1]);
                mixed.push_back(this->_mesh->elems[ernl*i + 2]);
                mixed.push_back(this->_mesh->elems[ernl*i + 3]);
                mixed.push_back(this->_mesh->elems[ernl*i + 4]);
                this->_quads_np++;
            }
        }

        write_hdf5_dataset_compressed(
                                        file,
                                        "/Mesh/MixedCells",
                                        int_h5,
                                        {static_cast<hsize_t>(mixed.size())},
                                        mixed.data( ),
                                        0
                                    );

        H5Gclose( meshGroup );
        H5Fclose( file );
    }

    MPI_Barrier( this->_mpi_config->mpi_comm );
}


void HDF5TimeSeriesExporterMPI::initialize_time_series()
{
    if ( this->_mpi_config->is_root( ) )
    {
        hid_t file = H5Fcreate(
                                    this->_field_file.c_str( ), 
                                    H5F_ACC_TRUNC,
                                    H5P_DEFAULT, 
                                    H5P_DEFAULT
                                );

        hid_t gid = H5Gcreate( file, "/Fields", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
        H5Gclose( gid );

        H5Fclose( file );
    }

    MPI_Barrier( this->_mpi_config->mpi_comm );
}


void HDF5TimeSeriesExporterMPI::append_time( cusfloat append_value )
{
    if ( this->_mpi_config->is_root( ) )
    {
        this->_time.push_back( append_value );
        this->_steps++;
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

    MPI_Barrier(this->_mpi_config->mpi_comm);
}


void HDF5TimeSeriesExporterMPI::write_xdmf( void )
{
    if ( !( this->_mpi_config->is_root( ) ) ) return;

    // Define auxiliary variables
    std::string attr_type;

    // Open file unit for XDMF writing (rank 0 only in MPI context)
    std::ofstream x( this->_xdmf_file );

    // Write XDMF header and domain
    x << "<?xml version=\"1.0\" ?>\n";
    x << "<Xdmf Version=\"2.0\">\n";
    x << "  <Domain>\n";
    x << "    <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">\n";

    for (size_t t = 0; t < this->_steps; ++t)
    {
        x << "      <Grid Name=\"Step" << t << "\" GridType=\"Uniform\">\n";
        x << "        <Time Value=\"" << this->_time[t] << "\"/>\n";

        if (t == 0)
        {
            // -------- MIXED TOPOLOGY --------
            x << "        <Topology TopologyType=\"Mixed\" NumberOfElements=\""
            << this->_mesh->elems_np << "\" "
            << "BaseOffset=\"0\">\n";
            x << "          <DataItem Format=\"HDF\" Dimensions=\""
            << (this->_tri_np*4 + this->_quads_np*5) << "\">\n";
            x << "            mesh.h5:/Mesh/MixedCells\n";
            x << "          </DataItem>\n";
            x << "        </Topology>\n";

            // -------- GEOMETRY --------
            x << "        <Geometry GeometryType=\"XYZ\">\n";
            x << "          <DataItem Format=\"HDF\" Dimensions=\"" << this->_mesh->nodes_np << " 3\">\n";
            x << "            mesh.h5:/Mesh/Points\n";
            x << "          </DataItem>\n";
            x << "        </Geometry>\n";
        }
        else
        {
            x << "        <Topology Reference=\"/Xdmf/Domain/Grid/Grid[1]/Topology\"/>\n";
            x << "        <Geometry Reference=\"/Xdmf/Domain/Grid/Grid[1]/Geometry\"/>\n";
        }

        // Hyperslab selection in XDMF
        for ( auto &c : this->_field_dims )
        {
            // Determine attribute type based on number of components
            if ( c.second[2] == 1 )
            {
                attr_type = "Scalar";
            }
            else if ( c.second[2] == 3 )
            {
                attr_type = "Vector";
            }
            else
            {
                throw std::runtime_error( "Unsupported number of components in field '" + c.first + "'. Only scalar (1) and vector (3) fields are supported." );
            }
            
            // Write attribute section
            x << "        <Attribute Name=\"" << c.first << "\" AttributeType=\"" << attr_type << "\" Center=\"Node\">\n";
            x << "          <DataItem ItemType=\"HyperSlab\" Dimensions=\"" << this->_mesh->nodes_np << " " << c.second[2] << "\">\n";
            x << "            <DataItem Dimensions=\"3 3\" Format=\"XML\">\n";
            x << "              " << t << " 0 0   1 1 1   1 " << this->_mesh->nodes_np << " " << c.second[2] << "\n";
            x << "            </DataItem>\n";
            x << "            <DataItem Format=\"HDF\" Dimensions=\"" << this->_steps << " " << this->_mesh->nodes_np << " " << c.second[2] << "\">\n";
            x << "              fields.h5:/Fields/" + c.first + "\n";
            x << "            </DataItem>\n";
            x << "          </DataItem>\n";
            x << "        </Attribute>\n";
        }

        x << "      </Grid>\n";
    }

    x << "    </Grid>\n";
    x << "  </Domain>\n";
    x << "</Xdmf>\n";
}


void HDF5TimeSeriesExporterMPI::ensureGroup(hid_t file, const std::string& path)
{
    hid_t gid = H5Gopen(file, path.c_str(), H5P_DEFAULT);
    if (gid < 0)
        gid = H5Gcreate(file, path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(gid);
}


void HDF5TimeSeriesExporterMPI::append_field_dataset(
                                                        hid_t               file,
                                                  const std::string&        path,
                                                  const cusfloat*           data,
                                                        hsize_t             components
                                                    )
{
    hid_t dset = H5Dopen(file, path.c_str(), H5P_DEFAULT);

    hid_t filespace = H5Dget_space(dset);
    hsize_t dims[3];
    H5Sget_simple_extent_dims(filespace, dims, nullptr);
    dims[0] += 1;
    H5Dset_extent(dset, dims);

    H5Sclose(filespace);
    filespace = H5Dget_space(dset);
    hsize_t start[3] = {dims[0] - 1, this->_node_offset, 0};
    hsize_t count[3] = {1, this->_local_np, components};

    hid_t memspace = H5I_INVALID_HID;

    if (this->_local_np == 0)
    {
        hsize_t dummy_dims[3] = {1, 1, components};
        memspace = H5Screate_simple(3, dummy_dims, nullptr);
        H5Sselect_none(filespace);
        H5Sselect_none(memspace);
    }
    else
    {
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, nullptr, count, nullptr);
        memspace = H5Screate_simple(3, count, nullptr);
    }

    hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);

    H5Dwrite(dset, cusfloat_h5, memspace, filespace, dxpl, data);

    H5Pclose(dxpl);
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Dclose(dset);
}
