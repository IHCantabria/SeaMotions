
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
#include <cstddef>
#include "hdf5.h"
#include <string>
#include <unordered_map>
#include <vector>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/math/custensor/custensor.hpp"
#include "../../src/mesh/mesh.hpp"


/**
 * @brief Abstract base class that centralizes the shared logic for HDF5/XDMF time-series exporters.
 * @details The class keeps the bookkeeping required to emit SeaMotions meshes and transient fields, while
 *          derived implementations provide the concrete I/O backend (serial ZLIB or MPI-IO). The helper
 *          methods offer reusable building blocks—mesh serialization, dataset registration, timeline tracking
 *          and XDMF emission—so specialised exporters only implement the transport-specific details.
 */
class HDF5TimeSeriesExporterBase
{
public:
    /** @brief Convenience alias that stores the dimensions (time, nodes, components) for every field. */
    using FieldDimensions      = std::vector<hsize_t>;
    /** @brief Map that associates each registered field name with its dataset dimensions. */
    using FieldDimensionsMap   = std::unordered_map<std::string, FieldDimensions>;

    virtual ~HDF5TimeSeriesExporterBase() = default;

    /** @brief Register a field dataset to be grown along time. */
    virtual void add_field(std::string field_name, hsize_t comps_np) = 0;
    /** @brief Append the latest timestep data for a specific field. */
    virtual void append_step(std::string field_name, cut::CusTensor<cusfloat>* field_data) = 0;
    /** @brief Record the simulation time corresponding to the next append_step invocation. */
    virtual void append_time(cusfloat append_value) = 0;
    /** @brief Create or refresh the extendable datasets that host the transient fields. */
    virtual void initialize_time_series() = 0;
    /** @brief Serialize the static mesh geometry/topology artefacts required by the XDMF description. */
    virtual void write_mesh() = 0;
    /** @brief Emit the XDMF metadata file referencing the mesh and time-dependent datasets. */
    virtual void write_xdmf() = 0;

protected:
    /**
     * @brief Build the exporter base by binding a mesh and output directory.
     * @param folder_path Destination directory that will host mesh, field, and XDMF files.
     * @param mesh Pointer to the SeaMotions mesh whose nodes and topology are exported.
     */
    HDF5TimeSeriesExporterBase(std::string folder_path, Mesh* mesh);

    /**
     * @brief Store metadata for a newly declared field and return its dataset dimensions.
     * @param field_name Dataset identifier (stored under /Fields in the HDF5 file).
     * @param comps_np Number of components per node (1 for scalars, 3 for vectors, ...).
     * @return Vector describing the hyperslab being appended at each timestep (time, nodes, components).
     */
    FieldDimensions register_field_metadata(const std::string& field_name, hsize_t comps_np);

    /**
     * @brief Append a simulation time entry and advance the internal step counter.
     * @param append_value Time value (seconds) associated with the next append_step() call.
     */
    void record_time_value(cusfloat append_value);

    /**
     * @brief Serialize the mesh geometry and Mixed topology into the configured mesh file.
     * @details The function writes /Mesh/Points and /Mesh/MixedCells datasets, updating triangle/quadrilateral
     *          counters so the XDMF metadata can build accurate topology descriptions.
     */
    void write_mesh_payload();

    /**
     * @brief Create (or truncate) the HDF5 container that stores extendable field datasets.
     * @details The function ensures the target file exists and that a /Fields group is available before
     *          derived classes attempt to append datasets.
     */
    void create_fields_container();

    /**
     * @brief Emit the XDMF document that references the mesh/field files and recorded timesteps.
     */
    void write_xdmf_payload() const;

    /** @brief Expose the bound mesh pointer to derived exporters. */
    Mesh* mesh() const noexcept { return this->_mesh; }
    /** @brief Return the total number of nodes in the bound mesh. */
    std::size_t nodes_count() const noexcept;

    /** @brief Provide mutable access to the field-dimensions map for advanced backends. */
    FieldDimensionsMap& field_dimensions() noexcept { return this->_field_dims; }
    /** @brief Provide read-only access to the field-dimensions map. */
    const FieldDimensionsMap& field_dimensions() const noexcept { return this->_field_dims; }

    /** @brief Return the full path to the mesh HDF5 file. */
    const std::string& mesh_file() const noexcept { return this->_mesh_file; }
    /** @brief Return the full path to the field HDF5 file. */
    const std::string& field_file() const noexcept { return this->_field_file; }
    /** @brief Return the full path to the XDMF metadata file. */
    const std::string& xdmf_file() const noexcept { return this->_xdmf_file; }

protected:
    /** @brief Pointer to the mesh owning the geometry and topology description. */
    Mesh*                           _mesh       = nullptr;
    /** @brief Destination directory shared by the mesh, field, and XDMF artefacts. */
    std::string                     _folder_path;
    /** @brief Absolute path to the HDF5 file that stores transient fields. */
    std::string                     _field_file;
    /** @brief List of registered field names for XDMF attribute emission. */
    std::vector<std::string>        _field_names;
    /** @brief Mapping between field names and dataset dimensions (time, nodes, components). */
    FieldDimensionsMap              _field_dims;
    /** @brief Absolute path to the static mesh HDF5 file. */
    std::string                     _mesh_file;
    /** @brief Number of quadrilateral elements encountered while writing the mesh. */
    std::size_t                     _quads_np   = 0;
    /** @brief Number of timesteps that have been appended so far. */
    std::size_t                     _steps      = 0;
    /** @brief Timeline values associated with each appended timestep. */
    std::vector<cusfloat>           _time;
    /** @brief Number of triangular elements encountered while writing the mesh. */
    std::size_t                     _tri_np     = 0;
    /** @brief Absolute path to the XDMF metadata file. */
    std::string                     _xdmf_file;
};