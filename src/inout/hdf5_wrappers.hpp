
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
#include "hdf5.h"
#include <vector>


inline void append_to_hdf5_dataset(
                                                hid_t                   file,
                                        const   std::string&            path,
                                                hid_t                   dtype,
                                        const   std::vector<hsize_t>&   count,
                                        const   void*                   data
                                    )
{
    // Open dataset
    hid_t dset = H5Dopen(file, path.c_str(), H5P_DEFAULT);

    // Get current dataset dimensions
    hid_t filespace = H5Dget_space(dset);
    std::vector<hsize_t> dims(count.size());
    H5Sget_simple_extent_dims(filespace, dims.data(), nullptr);

    // Extend dataset by one timestep
    dims[0] += 1;
    H5Dset_extent(dset, dims.data());

    // Close and reopen filespace to reflect the new dataset dimensions
    H5Sclose(filespace);
    filespace = H5Dget_space(dset);

    // Select hyperslab for the new timestep
    std::vector<hsize_t> offset(count.size(), 0);
    offset[0] = dims[0] - 1;

    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset.data(), nullptr,
                        count.data(), nullptr);

    // Create memory space for the chunk
    hid_t memspace = H5Screate_simple(count.size(), count.data(), nullptr);

    // Write data to the selected hyperslab
    H5Dwrite(dset, dtype, memspace, filespace, H5P_DEFAULT, data);

    // Close memory space, file space and dataset
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Dclose(dset);
}


inline void create_hdf5_dataset_simple( 
                                                    hid_t    loc_id, 
                                            const   char*    name, 
                                                    hsize_t  rank, 
                                            const   hsize_t* dims, 
                                                    hid_t    data_type 
                                    )
{
    // Create simple dataspace to define dataset dimensions
    hid_t dataspace = H5Screate_simple( 
                                            static_cast<int>( rank ), 
                                            dims, 
                                            nullptr 
                                        );

    // Create dataset
    hid_t dataset   = H5Dcreate2( 
                                    loc_id, 
                                    name, 
                                    data_type, 
                                    dataspace, 
                                    H5P_DEFAULT, 
                                    H5P_DEFAULT, 
                                    H5P_DEFAULT 
                                );

    // Close dataset and dataspace
    H5Dclose( dataset );
    H5Sclose( dataspace );
}


inline void create_hdf5_extendable_dataset(
                                                    hid_t                   file,
                                            const   std::string&            path,
                                                    hid_t                   dtype,
                                            const   std::vector<hsize_t>&   initial,
                                            const   std::vector<hsize_t>&   maxdims
                                        )
{
    hid_t space = H5Screate_simple(initial.size(), initial.data(), maxdims.data());

    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);

    std::vector<hsize_t> chunk = initial;
    chunk[0] = 1; // one timestep per chunk
    for (auto& c : chunk)
        if (c == 0) c = 1;

    H5Pset_chunk(plist, chunk.size(), chunk.data());
    H5Pset_deflate(plist, 4);

    hid_t dset = H5Dcreate(file, path.c_str(), dtype, space,
                            H5P_DEFAULT, plist, H5P_DEFAULT);

    H5Dclose(dset);
    H5Pclose(plist);
    H5Sclose(space);
}


inline void write_hdf5_dataset_chunk(
                                                hid_t     loc_id,
                                        const   char*     name,
                                                hsize_t   rank,
                                        const   hsize_t*  dataset_shape,
                                        const   hsize_t*  chunk_shape,
                                        const   hsize_t*  offset,
                                        const   void*     data,
                                                hid_t     data_type
                                    )
{
    // Open dataset
    static_cast<void>( dataset_shape );
    hid_t dataset    = H5Dopen2( loc_id, name, H5P_DEFAULT );

    // Select hyperslab in the file
    hid_t file_space = H5Dget_space( dataset );
    H5Sselect_hyperslab( file_space, H5S_SELECT_SET, offset, nullptr, chunk_shape, nullptr );

    // Create memory space for the chunk
    hid_t mem_space = H5Screate_simple( static_cast<int>( rank ), chunk_shape, nullptr );

    // Write data to the selected hyperslab
    H5Dwrite( dataset, data_type, mem_space, file_space, H5P_DEFAULT, data );

    // Close memory space, file space and dataset
    H5Sclose( mem_space );
    H5Sclose( file_space );
    H5Dclose( dataset );
}


inline void write_hdf5_dataset_compressed(
                                                    hid_t                   file,
                                            const   std::string&            path,
                                                    hid_t                   dtype,
                                            const   std::vector<hsize_t>&   dims,
                                            const   void*                   data,
                                                    int                     compression = 4
                                        )
{
    // Cretae dataspace
    hid_t space = H5Screate_simple(dims.size(), dims.data(), nullptr);

    // Create property list for dataset creation and set chunking and compression
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist, dims.size(), dims.data());
    H5Pset_deflate(plist, compression);

    // Create dataset
    hid_t dset = H5Dcreate(file, path.c_str(), dtype, space,
                            H5P_DEFAULT, plist, H5P_DEFAULT);

    // Write data
    H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    // Close dataset, property list and dataspace
    H5Dclose(dset);
    H5Pclose(plist);
    H5Sclose(space);
}


template<typename T>
inline void write_hdf5_dataset_mpi( 
                                            hid_t                   file,
                                    const   std::string&            path,
                                            hid_t                   dtype,
                                    const   std::vector<hsize_t>&   dims,
                                    const   std::vector<hsize_t>&   offset,
                                    const   std::vector<hsize_t>&   chunk_dims,
                                    const   T*                      data
                                )
{
    // Create global dataspace
    hid_t space     = H5Screate_simple(
                                            dims.size( ), 
                                            dims.data( ), 
                                            nullptr
                                        );

    // Create datasets for points and connectivity
    hid_t dset      = H5Dcreate(
                                    file, 
                                    path.c_str( ), 
                                    dtype,
                                    space, 
                                    H5P_DEFAULT, 
                                    H5P_DEFAULT, 
                                    H5P_DEFAULT
                                );

    // Select hyperslab for this rank
    // hsize_t start[2]    = {0, 0};
    // hsize_t count[2]    = {this->_mesh->nodes_np, 3};

    H5Sselect_hyperslab( 
                            space, 
                            H5S_SELECT_SET,
                            offset.data( ), 
                            nullptr,
                            chunk_dims.data( ), 
                            nullptr
                        );

    // Select memory space
    hid_t memspace      = H5Screate_simple( dims.size( ), chunk_dims.data(), nullptr );

    hid_t dxpl          = H5Pcreate( H5P_DATASET_XFER );
    H5Pset_dxpl_mpio( dxpl, H5FD_MPIO_COLLECTIVE );

    // Write data to the selected hyperslab
    H5Dwrite( 
                dset, 
                dtype, 
                memspace, 
                space, 
                dxpl, 
                data 
            );
    
    // Close dataset, property list and dataspace
    H5Pclose( dxpl      );
    H5Sclose( memspace  );
    H5Dclose( dset      );
    H5Sclose( space     );

}


template<typename T>
inline void write_hdf5_scalar_attribute( 
                                                    hid_t   loc_id, 
                                            const   char*   name, 
                                                    hid_t   data_type, 
                                                    T       value 
                                        )
{
    // Create simple dataspace for scalar attribute
    hid_t dataspace = H5Screate( H5S_SCALAR );

    // Create attribute and write value
    hid_t attr      = H5Acreate2( loc_id, name, data_type, dataspace, H5P_DEFAULT, H5P_DEFAULT );
    H5Awrite( attr, data_type, &value );

    // Close attribute and dataspace
    H5Aclose( attr );
    H5Sclose( dataspace );
}