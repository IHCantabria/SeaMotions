
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
#include "../containers/field_mesh_data_config.hpp"
#include "../inout/hdf5_time_series_exporter_interface.hpp"
#include "../mesh/mesh.hpp"
#include "../math/custensor/custensor.hpp"
#include "../tools/memory_report.hpp"


template<typename T, typename ModeComp>
struct FieldMeshData
{
private:
    /*** Declare private class attributes ***/
    int                             _compresion_level   = 0;        // Compression level for output files (0-9, where 0 is no compression and 9 is maximum compression). Not used if MPI is enabled.
    cut::CusTensor<cusfloat>        _data_scalar_r;                 // Scalar data stored in the root process
    cut::CusTensor<cusfloat>        _data_vector;                   // Vector data buffer for local packing
    cut::CusTensor<cusfloat>        _data_vector_r;                 // Vector data stored in the root process
    HDF5TimeSeriesExporterBase*     _exporter           = nullptr;  // Pointer to the HDF5 exporter interface for writing mesh and field data
    FieldPointsDef*                 _field_points_def   = nullptr;  // Pointer to the field points definition
    Input*                          _input              = nullptr;  // Pointer to the input data
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

    FieldMeshData( 
                        Input*          input,
                        FieldPointsDef* field_points_def,
                        std::string     out_folder_path,
                        MpiConfig*      mpi_config 
                    );

    ~FieldMeshData( );

    /*** Declare public methods ***/
    void add_step( cusfloat freq );

    template<FieldTypeE field_type, FieldComponentE field_comp, ComplexDataTypeE complex_data_type>
    cut::CusTensor<cusfloat>* get_field_data_scalar( std::size_t heading_index );

    template<FieldTypeE field_type, FieldComponentE field_comp, ComplexDataTypeE complex_data_type>
    cut::CusTensor<cusfloat>* get_field_data_vector( std::size_t heading_index );

    RadDiffData<T, ModeComp>* get_rdd( void ) const;

    std::size_t memory_bytes( void ) const;

    void append_memory_report(
                                    std::vector<MemoryReportEntry>& entries,
                                    const std::string& prefix
                                ) const;

};

// Include template implementation file
#include "field_mesh_data.txx"