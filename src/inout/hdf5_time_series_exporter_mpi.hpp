
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

#include <mpi.h>
#include <hdf5.h>
#include <cstddef>
#include <string>
#include <vector>

#include "hdf5_time_series_exporter_interface.hpp"
#include "../../src/containers/mpi_config.hpp"
#include "../../src/inout/hdf5_wrappers.hpp"
#include "../../src/math/custensor/custensor.hpp"
#include "../../src/mesh/mesh.hpp"


/**
 * @brief MPI-enabled exporter that writes SeaMotions meshes and distributed field time series to HDF5/XDMF.
 * @details The class orchestrates collective filesystem access so each rank contributes its portion of every
 *          field, while rank 0 keeps the timeline and companion XDMF description consistent with the files.
 *          It reuses the shared helpers implemented by HDF5TimeSeriesExporterBase and only overrides the
 *          transport-specific details (collective dataset creation, MPI barriers, and rank ownership).
 */
class HDF5TimeSeriesExporterMPI final : public HDF5TimeSeriesExporterBase
{
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

    /** @copydoc HDF5TimeSeriesExporterBase::add_field */
    void add_field(
                        std::string             field_name,
                        hsize_t                 comps_np
                    ) override;

    /** @copydoc HDF5TimeSeriesExporterBase::write_mesh */
    void write_mesh( void ) override;

    /** @copydoc HDF5TimeSeriesExporterBase::initialize_time_series */
    void initialize_time_series() override;

    /** @copydoc HDF5TimeSeriesExporterBase::append_time */
    void append_time( cusfloat append_value ) override;

    /** @copydoc HDF5TimeSeriesExporterBase::append_step */
    void append_step(
                        std::string                 field_name,
                        cut::CusTensor<cusfloat>*   field_data
                    ) override;

    /** @copydoc HDF5TimeSeriesExporterBase::write_xdmf */
    void write_xdmf( void ) override;

private:
    /**
     * @brief Ensure the requested group exists, creating it when missing.
     * @param file Open HDF5 file where the group lives.
     * @param path Absolute group path inside the file.
     */
    static void ensureGroup(hid_t file, const std::string& path);

private:
    /** @brief Pointer to the MPI configuration describing communicator layout. */
    MpiConfig*  _mpi_config = nullptr;
    /** @brief Total number of mesh nodes across all ranks. */
    std::size_t _global_np  = 0;
    /** @brief Number of nodes handled by this rank when writing hyperslabs. */
    std::size_t _local_np   = 0;
    /** @brief Global index of the first node owned by this rank. */
    std::size_t _node_offset = 0;
};
