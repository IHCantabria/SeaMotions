
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
#include <mpi.h>
#include <hdf5.h>
#include <vector>
#include <string>
#include <unordered_map>
#include <cstddef>

// Include local module
#include "../../src/config.hpp"
#include "../../src/containers/mpi_config.hpp"
#include "../../src/inout/hdf5_wrappers.hpp"
#include "../../src/math/custensor/custensor.hpp"
#include "../../src/mesh/mesh.hpp"


/**
 * @brief MPI-enabled exporter that writes SeaMotions meshes and distributed field time series to HDF5/XDMF.
 * @details The class orchestrates collective filesystem access so each rank contributes its portion of every
 *          field, while rank 0 keeps the timeline and companion XDMF description consistent with the files.
 */
class HDF5TimeSeriesExporterMPI
{
private:
    /** @brief Convenience alias for storing dataset dimensions per field name. */
    using dict_sst = std::unordered_map<std::string, std::vector<hsize_t>>;

    /*** Declare private variables ***/
    std::vector<std::string>    _field_names    ;           // Names of fields added (for XDMF reference)
    dict_sst                    _field_dims     ;           // Dimensions of each field (for XDMF reference)
    std::string                 _field_file     ;           // Path to HDF5 file for time‑series fields
    std::string                 _folder_path    ;           // Base folder for output files (mesh, fields, xdmf)
    std::size_t                 _global_np      = 0;        // Number of global nodes. Total number of nodes in the mesh (for global dataset dimensions)
    std::size_t                 _local_np       = 0;        // Number of local nodes. Those living in this processor memory (for hyperslab selection)
    Mesh*                       _mesh           = nullptr;  // Pointer to mesh object (assumed to be initialized and containing mesh info)
    std::string                 _mesh_file      ;           // Path to HDF5 file for mesh (static, written once)
    MpiConfig*                  _mpi_config     = nullptr;  // Pointer to MPI configuration (assumed to be initialized with rank, size, comm)
    std::size_t                 _quads_np       = 0;        // Number of quads (for connectivity dataset)
    std::size_t                 _steps          = 0;        // Number of time steps written (for XDMF reference)
    std::vector<cusfloat>       _time           ;           // Vector to store time values for each step (for XDMF reference)
    std::size_t                 _tri_np         = 0;        // Number of triangles (for connectivity dataset)
    std::string                 _xdmf_file      ;           // Path to XDMF file (rank 0 only, references mesh and field files)

    std::size_t                 _node_offset    = 0;        // Starting global node index assigned to this rank
public:
    /**
     * @brief Construct an MPI exporter bound to a rank-partitioned mesh.
     * @param folder_path Output directory that will receive mesh, fields, and XDMF files.
     * @param mesh Pointer to the global mesh instance (coordinates accessible on rank 0).
     * @param mpi_config Pointer to the MPI configuration exposing communicator utilities.
     */
    HDF5TimeSeriesExporterMPI(
                                const   std::string&    folder_path,
                                        Mesh*           mesh,
                                        MpiConfig*      mpi_config
                            );

    /**
     * @brief Declare a field dataset to be grown along the time dimension.
     * @param field_name Logical dataset name stored under /Fields.
     * @param comps_np Number of components per node (1 for scalars, 3 for vectors, ...).
     */
    void add_field(
                        std::string             field_name,
                        hsize_t                 comps_np
                    );

    /**
     * @brief Serialize the static mesh geometry and mixed connectivity (root rank only writes).
     */
    void write_mesh( void );

    /**
     * @brief Create (or truncate) the extendable field container shared by all ranks.
     */
    void initialize_time_series();

    /**
     * @brief Record the physical time associated with the next timestep.
     * @param append_value Time value stored only by the root rank for XDMF emission.
     */
    void append_time( cusfloat append_value );

    /**
     * @brief Append the local rank contribution for a specific field at the current timestep.
     * @param field_name Dataset identifier previously registered through add_field().
     * @param field_data Tensor owning the local node values laid out contiguously.
     */
    void append_step(
                        std::string                 field_name,
                        cut::CusTensor<cusfloat>*   field_data
                    );

    /**
     * @brief Write the XDMF companion file referencing mesh and field data (root rank only).
     */
    void write_xdmf( void );

private:
    /**
     * @brief Ensure the requested group exists, creating it when missing.
     * @param file Open HDF5 file where the group lives.
     * @param path Absolute group path inside the file.
     */
    static void ensureGroup(hid_t file, const std::string& path);

    /**
     * @brief Legacy helper that appends a field hyperslab owned by the current rank.
     * @param file Open HDF5 file.
     * @param path Dataset path relative to the root.
     * @param data Pointer to the local buffer (node-major layout).
     * @param components Number of components per node.
     */
    void append_field_dataset(hid_t file,
                              const std::string& path,
                              const cusfloat* data,
                              hsize_t components);
};