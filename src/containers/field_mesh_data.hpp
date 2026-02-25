
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

// Include local modules
#include "../containers/rad_diff_data.hpp"
#include "../containers/mpi_config.hpp"
#include "../inout/hdf5_time_series_exporter_interface.hpp"
#include "../mesh/mesh.hpp"


struct FieldMeshDataConfig
{
    int                     compression_level       = 0;        // Compression level for output files (0-9, where 0 is no compression and 9 is maximum compression). Not used if MPI is enabled.
    std::size_t             dofs_np                 = 0;        // Number of degrees of freedom for the field data (e.g., 1 for scalar fields, 3 for vector fields, etc.)
    std::size_t             freqs_np                = 0;        // Number of frequencies for the field data
    std::vector<cusfloat*>  headings                = {};       // Array of heading values for the field data
    std::size_t             headings_np             = 0;        // Number of headings for the field data
    std::string             mesh_file_path          ;           // Path to the mesh file containing the field points
    std::string             body_name               ;           // Name of the body to extract from the mesh file
    bool                    out_components          = false;    // Flag to indicate if field components data should be outputted
    bool                    out_parallel            = false;    // Flag to indicate if output files should be written in parallel using MPI. If true, compression_level is ignored and no compression is applied.
    bool                    out_potential           = false;    // Flag to indicate if potential field data should be outputted
    bool                    out_pressure            = false;    // Flag to indicate if pressure field data should be outputted
    bool                    out_velocity            = false;    // Flag to indicate if velocity field data should be outputted
};


template<typename T, typename ModeComp>
struct FieldMeshData
{
private:
    /*** Declare private class attributes ***/
    int                             _compresion_level   = 0;        // Compression level for output files (0-9, where 0 is no compression and 9 is maximum compression). Not used if MPI is enabled.
    cusfloat*                       _data_scalar_r      = nullptr;  // Pointer to the scalar data storaged in the root process
    cusfloat*                       _data_vector        = nullptr;  // Pointer to the vector data
    cusfloat*                       _data_vector_r      = nullptr;  // Pointer to the vector data storaged in the root process
    HDF5TimeSeriesExporterBase*     _exporter           = nullptr;  // Pointer to the HDF5 exporter interface for writing mesh and field data
    bool                            _is_data_heap       = false;    // Flag to indicate if memory is allocated on heap for the root data
    bool                            _is_heap_exporter   = false;    // Flag to indicate if memory is allocated on heap for the exporter interface
    bool                            _is_heap_mesh       = false;    // Flag to indicate if memory is allocated on heap               
    bool                            _is_heap_rdd        = false;    // Flag to indicate if memory is allocated on heap               
    Mesh*                           _mesh               = nullptr;  // Pointer to the mesh containing the field points
    MpiConfig*                      _mpi_config         = nullptr;  // Pointer to MPI configuration
    RadDiffData<T, ModeComp>*       _rdd                = nullptr;  // Pointer to the radiation and diffraction data
    bool                            _out_components     = false;    // Flag to indicate if field components data should be outputted
    bool                            _out_parallel       = false;    // Flag to indicate if output files should be written in parallel using MPI. If true, compression_level is ignored and no compression is applied.
    bool                            _out_potential      = false;    // Flag to indicate if potential field data should be outputted
    bool                            _out_pressure       = false;    // Flag to indicate if pressure field data should be outputted
    bool                            _out_velocity       = false;    // Flag to indicate if velocity field data should be outputted
    std::size_t                     _step               = 0;        // Current time/freq step index for appending field data

    /*** Declare private methods ***/
    void    _configure_exporter( void );

public:
    /*** Declare class constructor ***/
    FieldMeshData( ) = default;

    FieldMeshData( FieldMeshDataConfig& config, MpiConfig* mpi_config );

    ~FieldMeshData( );

    /*** Declare public methods ***/
    void add_step( cusfloat freq );

    template<FieldTypeE field_type, FieldComponentE field_comp>
    const cusfloat* get_field_data_scalar( std::size_t heading_index ) const;

    template<FieldTypeE field_type, FieldComponentE field_comp>
    const cusfloat* get_field_data_vector( std::size_t heading_index ) const;

};

// Include template implementation file
#include "field_mesh_data.txx"