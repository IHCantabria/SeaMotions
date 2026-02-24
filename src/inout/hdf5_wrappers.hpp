
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

// Include local modules
#include "../config.hpp"


/**
 * @brief Append a contiguous block of samples to an extendable dataset.
 *
 * @tparam T   Plain-old-data type matching the dataset layout.
 * @param file Open HDF5 file identifier that owns the dataset.
 * @param path Absolute dataset path (e.g., "/Fields/Pressure").
 * @param dtype Memory datatype used to interpret @p data.
 * @param count Shape of the block to append; first entry corresponds to the time dimension.
 * @param data  Pointer to the samples to be written.
 */
template<typename T>
inline void append_to_hdf5_dataset(
                                                hid_t                   file,
                                        const   std::string&            path,
                                                hid_t                   dtype,
                                        const   std::vector<hsize_t>&   count,
                                        const   T*                      data
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


/**
 * @brief Append a timestep to a dataset using collective MPI-IO semantics.
 *
 * @tparam T        POD type of the values written from each rank.
 * @param file      Shared file identifier opened with an MPI-enabled access property list.
 * @param path      Dataset path relative to @p file.
 * @param data      Pointer to the local rank contribution.
 * @param local_np  Number of nodes owned by the current rank.
 * @param offset    Global index of the first node owned by the current rank.
 * @param components Number of components per node (1 for scalar fields, 3 for vectors, ...).
 */
template<typename T>
void append_to_hdf5_dataset_mpi(
                                            hid_t           file,
                                    const   std::string&    path,
                                    const   T*              data,
                                    const   hsize_t         local_np,
                                    const   hsize_t         offset,
                                    const   hsize_t         components
                                )
{
    // Open dataset
    hid_t dset = H5Dopen( 
                            file, 
                            path.c_str(), 
                            H5P_DEFAULT 
                        );

    // Get current dataset dimensions
    hid_t filespace = H5Dget_space( dset );

    // Extend dataset by one timestep
    hsize_t dims[3];
    H5Sget_simple_extent_dims( filespace, dims, nullptr );
    dims[0] += 1;
    H5Dset_extent( dset, dims );

    // Close and reopen filespace to reflect the new dataset dimensions
    H5Sclose( filespace );
    filespace = H5Dget_space( dset );

    // Select hyperslab for the new timestep
    hsize_t start[3] = { dims[0] - 1,      offset,          0 };
    hsize_t count[3] = {           1,    local_np, components };
    hid_t   memspace = H5I_INVALID_HID;

    if ( local_np == 0 )
    {
        hsize_t dummy_dims[3] = { 1, 1, components };
        memspace = H5Screate_simple(3, dummy_dims, nullptr);
        H5Sselect_none(filespace);
        H5Sselect_none(memspace);
    }
    else
    {
        H5Sselect_hyperslab(
                                filespace, 
                                H5S_SELECT_SET, 
                                start, 
                                nullptr, 
                                count, 
                                nullptr
                            );

        memspace    = H5Screate_simple( 3, count, nullptr );
    }

    // Write data to the selected hyperslab collectively
    hid_t dxpl  = H5Pcreate( H5P_DATASET_XFER );

    H5Pset_dxpl_mpio( 
                        dxpl, 
                        H5FD_MPIO_COLLECTIVE 
                    );

    H5Dwrite( 
                dset, 
                cusfloat_h5, 
                memspace, 
                filespace, 
                dxpl, 
                data 
            );

    // Close property list, memory space, file space and dataset
    H5Pclose( dxpl      );
    H5Sclose( memspace  );
    H5Sclose( filespace );
    H5Dclose( dset      );

}


/**
 * @brief Create a fixed-size dataset with a simple dataspace.
 *
 * @param loc_id    Identifier of the parent group or file.
 * @param name      Dataset name within @p loc_id.
 * @param rank      Number of dimensions in @p dims.
 * @param dims      Array describing the extent of each dimension.
 * @param data_type HDF5 datatype used for storage.
 */
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


/**
 * @brief Create an extendable chunked dataset with optional compression.
 *
 * @param file     File identifier that will host the dataset.
 * @param path     Full dataset path.
 * @param dtype    Storage datatype.
 * @param initial  Initial dataspace dimensions (first entry usually zero for time).
 * @param maxdims  Maximum allowed extent per dimension (use H5S_UNLIMITED for unlimited growth).
 */
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


/**
 * @brief Create an extendable dataset intended to be accessed collectively through MPI-IO.
 *
 * @param file       Shared file identifier opened with MPI capabilities.
 * @param path       Dataset path relative to the file root.
 * @param dtype      Storage datatype.
 * @param dset_dims  Initial dataspace dimensions.
 * @param max_dims   Maximum extent for each dimension.
 * @param chunk_dims Chunking strategy used by the storage layout.
 */
inline void create_hdf5_extendable_dataset_mpi(
                                                            hid_t                   file,
                                                    const   std::string&            path,
                                                            hid_t                   dtype,
                                                    const   std::vector<hsize_t>&   dset_dims,
                                                    const   std::vector<hsize_t>&   max_dims,
                                                    const   std::vector<hsize_t>&   chunk_dims
                                                )
{
    // Create simpledataspace to define dataset dimensions
    hid_t   space   = H5Screate_simple(
                                            dset_dims.size( ), 
                                            dset_dims.data( ), 
                                            max_dims.data( )
                                        );

    // Create property list for dataset creation and set chunking
    hid_t plist = H5Pcreate( H5P_DATASET_CREATE );
    H5Pset_chunk( plist, chunk_dims.size( ), chunk_dims.data( ) );

    // Create dataset
    hid_t dset = H5Dcreate(
                                file, 
                                path.c_str(), 
                                dtype,
                                space, 
                                H5P_DEFAULT, 
                                plist, 
                                H5P_DEFAULT
                            );

    // Close dataset, property list and dataspace
    H5Pclose(plist);
    H5Sclose(space);
    H5Dclose(dset);
}


/**
 * @brief Write a contiguous chunk into an existing dataset.
 *
 * @param loc_id        Parent group or file that contains the dataset.
 * @param name          Dataset name relative to @p loc_id.
 * @param rank          Number of dimensions in the hyperslab.
 * @param dataset_shape (Unused) Full dataset shape, provided for symmetry.
 * @param chunk_shape   Size of the block being written.
 * @param offset        Starting indices of the hyperslab.
 * @param data          Pointer to the data to copy into the file.
 * @param data_type     Memory datatype describing @p data.
 */
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


/**
 * @brief Create and populate a chunked dataset with optional zlib compression.
 *
 * @param file        File identifier where the dataset will be stored.
 * @param path        Full dataset path.
 * @param dtype       Storage datatype.
 * @param dims        Dataspace dimensions.
 * @param data        Pointer to the full dataset contents.
 * @param compression Deflate level (0 disables compression).
 */
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


/**
 * @brief Create a dataset and populate it through a single collective write.
 *
 * @tparam T         POD type stored in the dataset.
 * @param file       File identifier opened with an MPI-IO access property list.
 * @param path       Dataset path relative to the file root.
 * @param dtype      Storage datatype.
 * @param dims       Global dimensions of the dataset.
 * @param offset     Starting indices of the subsection handled by the current rank.
 * @param chunk_dims Size of the contiguous block provided by the current rank.
 * @param data       Pointer to the local data buffer.
 */
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


/**
 * @brief Attach a scalar attribute to a dataset or group.
 *
 * @tparam T       C++ scalar type compatible with @p data_type.
 * @param loc_id   Identifier of the object receiving the attribute.
 * @param name     Attribute name.
 * @param data_type HDF5 datatype describing @p value.
 * @param value    Attribute payload.
 */
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