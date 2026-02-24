
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

#include "hdf5_time_series_exporter_interface.hpp"

// Include local modules
#include "../../src/inout/hdf5_wrappers.hpp"
#include "../../src/math/custensor/custensor.hpp"
#include "../../src/math/math_tools.hpp"
#include "../../src/mesh/mesh.hpp"


/**
 * @brief Convenience exporter that serializes SeaMotions meshes and time-dependent fields to HDF5/XDMF.
 * @details The class owns the bookkeeping needed to create static mesh geometry, extendable field datasets,
 *          and synchronised XDMF descriptions so the generated files can be explored directly in ParaView.
 *          Typical usage creates the exporter with a mesh, registers each field via add_field(), then appends
 *          time steps with append_time() and append_step(). Finally write_xdmf() writes the companion XML file.
 *          The class leverages HDF5TimeSeriesExporterBase to reuse the mesh/XDMF bookkeeping shared with the
 *          MPI-capable exporter counterpart.
 */
class HDF5TimeSeriesExporter final : public HDF5TimeSeriesExporterBase
{
public:
    /*** Declare class constructors  ***/

    /**
     * @brief Build an exporter bound to a target folder and SeaMotions mesh.
     * @param folder_path Destination directory that will hold mesh, fields, and XDMF files.
     * @param mesh Pointer to the mesh whose geometry/topology are written once at construction.
     */
    HDF5TimeSeriesExporter(
                                const   std::string&    folder_path,
                                        Mesh*           mesh
                            );


    /***  Declare class methods  ***/

    /**
     * @brief Declare a new field dataset to be tracked along time.
     * @param field_name HDF5 dataset name stored under /Fields.
     * @param samples_multiplier Number of samples stored per mesh node (1 for classical node fields).
     * @param components_np Number of components in each sample (e.g., 1 for scalar, 3 for vector).
     * @details The first dimension is promoted to unlimited so append_step() can extend it at each timestep.
     */
    void add_field( 
                        std::string             field_name,
                        hsize_t                 samples_multiplier,
                        hsize_t                 components_np
                    ) override;

    /**
     * @brief Append a tensor snapshot for a previously registered field.
     * @param field_name Identifier passed to add_field().
     * @param field_data Tensor with data laid out in row-major order matching the registered dimensions.
     * @details The dataset is extended along the time axis and the new block is written via append_to_hdf5_dataset().
     */
    void append_step(
                        std::string                 field_name,
                        cut::CusTensor<cusfloat>*   field_data
                    ) override;

    /**
     * @brief Record the physical time corresponding to the next append_step() call.
     * @param append_value Time value in seconds (or consistent units) to store in the XDMF timeline.
     */
    void append_time( cusfloat append_value ) override;

    /**
     * @brief Create the HDF5 container that will hold extendable field datasets.
     * @details Called during construction to ensure /Fields exists before any field registrations occur.
     */
    void initialize_time_series(  ) override;

    /**
     * @brief Serialize the static mesh geometry and mixed topology once.
     * @details Writes /Mesh/Points and /Mesh/MixedCells datasets so the XDMF file can reference them for all timesteps.
     */
    void write_mesh( void ) override;

    /**
     * @brief Emit the XDMF companion file describing the mesh and field time series.
     * @details References the HDF5 datasets via HyperSlab selections so ParaView can stream the results without duplication.
     */
    void write_xdmf( void ) override;

};

/**
 * @example tests/tools/test_hdf5_paraview_class_time.cpp
 * @brief Demonstrates how to export SeaMotions panel meshes and synthetic pressure/velocity fields to ParaView.
 * @details The example builds a mesh from disk, registers pressure and velocity fields as extendable datasets,
 *          and advances a small harmonic solution in time before writing the final XDMF companion file. Running
 *          this translation unit produces the files mesh.h5, fields.h5, and data.xdmf2 inside the target folder.
 */
void example_hdf5_time_series_exporter( void );