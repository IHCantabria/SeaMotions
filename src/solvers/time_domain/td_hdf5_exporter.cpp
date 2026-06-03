/*
 * Copyright (c) 2025 Sergio Fernández Ruano / IHCantabria
 *
 * This file is part of SeaMotionsTimeDev.
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

#include "td_hdf5_exporter.hpp"
#include "../../inout/hdf5_wrappers.hpp"

#include <filesystem>
#include <iostream>
#include <sstream>

#include "../../config.hpp"


/*****************************************************************************
 * Destructor
 *****************************************************************************/

TimeDomainHDF5Exporter::~TimeDomainHDF5Exporter( )
{
    this->close( );
}


/*****************************************************************************
 * initialize
 *****************************************************************************/

void TimeDomainHDF5Exporter::initialize(
    const std::string&  folder_path,
    int                 n_panels,
    int                 n_bodies
)
{
    namespace fs = std::filesystem;

    _n_panels = n_panels;
    _n_bodies = n_bodies;

    // Build path: <case>/1_results/results.td.h5
    const fs::path results_dir = fs::path( folder_path ) / fs::path( RESULTS_FOLDER_NAME );
    const fs::path file_path   = results_dir / "results.td.h5";

    _file_path = file_path.string( );

    _file = H5Fcreate(
        _file_path.c_str( ),
        H5F_ACC_TRUNC,
        H5P_DEFAULT,
        H5P_DEFAULT
    );

    if ( _file == H5I_INVALID_HID )
    {
        std::cerr << "[ERROR] TimeDomainHDF5Exporter: could not create file: "
                  << file_path << std::endl;
        return;
    }

    // ---- Create top-level groups ----
    {
        hid_t g = H5Gcreate( _file, "/panels", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
        H5Gclose( g );
    }
    {
        hid_t g = H5Gcreate( _file, "/bodies", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
        H5Gclose( g );
    }

    // ---- /time  (1-D, unlimited) ----
    create_hdf5_extendable_dataset(
        _file, "/time", cusfloat_h5,
        { 0 },
        { H5S_UNLIMITED }
    );

    // ---- Panel pressure datasets  (2-D: n_steps × n_panels) ----
    if ( n_panels > 0 )
    {
        const hsize_t np = static_cast<hsize_t>( n_panels );
        create_hdf5_extendable_dataset(
            _file, "/panels/pressure",
            cusfloat_h5, { 0, np }, { H5S_UNLIMITED, np } );
        create_hdf5_extendable_dataset(
            _file, "/panels/phi_dt_term",
            cusfloat_h5, { 0, np }, { H5S_UNLIMITED, np } );
        create_hdf5_extendable_dataset(
            _file, "/panels/kinetic_term",
            cusfloat_h5, { 0, np }, { H5S_UNLIMITED, np } );
        create_hdf5_extendable_dataset(
            _file, "/panels/hydrostatic_term",
            cusfloat_h5, { 0, np }, { H5S_UNLIMITED, np } );
    }

    // ---- Body kinematics + force datasets  (2-D: n_steps × 6) ----
    for ( int ib = 0; ib < n_bodies; ib++ )
    {
        std::ostringstream grp;
        grp << "/bodies/body_" << ib;
        {
            hid_t g = H5Gcreate( _file, grp.str( ).c_str( ),
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
            H5Gclose( g );
        }

        create_hdf5_extendable_dataset(
            _file, grp.str( ) + "/position",
            cusfloat_h5, { 0, 6 }, { H5S_UNLIMITED, 6 } );
        create_hdf5_extendable_dataset(
            _file, grp.str( ) + "/velocity",
            cusfloat_h5, { 0, 6 }, { H5S_UNLIMITED, 6 } );
        create_hdf5_extendable_dataset(
            _file, grp.str( ) + "/acceleration",
            cusfloat_h5, { 0, 6 }, { H5S_UNLIMITED, 6 } );

        // Forces: total + three Bernoulli components
        const std::string fgrp = grp.str( ) + "/forces";
        {
            hid_t g = H5Gcreate( _file, fgrp.c_str( ),
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
            H5Gclose( g );
        }
        create_hdf5_extendable_dataset(
            _file, fgrp + "/total",
            cusfloat_h5, { 0, 6 }, { H5S_UNLIMITED, 6 } );
        create_hdf5_extendable_dataset(
            _file, fgrp + "/phi_dt_term",
            cusfloat_h5, { 0, 6 }, { H5S_UNLIMITED, 6 } );
        create_hdf5_extendable_dataset(
            _file, fgrp + "/radiation",
            cusfloat_h5, { 0, 6 }, { H5S_UNLIMITED, 6 } );
        create_hdf5_extendable_dataset(
            _file, fgrp + "/wave_excitation",
            cusfloat_h5, { 0, 6 }, { H5S_UNLIMITED, 6 } );
        create_hdf5_extendable_dataset(
            _file, fgrp + "/kinetic_term",
            cusfloat_h5, { 0, 6 }, { H5S_UNLIMITED, 6 } );
        create_hdf5_extendable_dataset(
            _file, fgrp + "/hydrostatic_term",
            cusfloat_h5, { 0, 6 }, { H5S_UNLIMITED, 6 } );
    }

    // Flush the superblock + dataset headers immediately so the file is a
    // valid HDF5 container even if the solver dies before the first
    // append_step / destructor call.
    H5Fflush( _file, H5F_SCOPE_GLOBAL );

    std::cout << " -> HDF5 time-series: " << file_path << std::endl;
}


/*****************************************************************************
 * append_step
 *****************************************************************************/

void TimeDomainHDF5Exporter::append_step(
    cusfloat                                        t,
    const cusfloat*                                 pressure,
    const cusfloat*                                 phi_dt_comp,
    const cusfloat*                                 kinetic_comp,
    const cusfloat*                                 hydrostatic_comp,
    int                                             n_panels,
    const std::vector<std::array<cusfloat, 6>>&     body_pos,
    const std::vector<std::array<cusfloat, 6>>&     body_vel,
    const std::vector<std::array<cusfloat, 6>>&     body_acc,
    const std::vector<std::array<cusfloat, 6>>&     body_force_total,
    const std::vector<std::array<cusfloat, 6>>&     body_force_phi_dt,
    const std::vector<std::array<cusfloat, 6>>&     body_force_radiation,
    const std::vector<std::array<cusfloat, 6>>&     body_force_wave,
    const std::vector<std::array<cusfloat, 6>>&     body_force_kinetic,
    const std::vector<std::array<cusfloat, 6>>&     body_force_hydrostatic
)
{
    if ( _file == H5I_INVALID_HID ) { return; }

    // ---- time ----
    append_to_hdf5_dataset( _file, "/time", cusfloat_h5, { 1 }, &t );

    // ---- panel fields ----
    if ( n_panels > 0 && n_panels == _n_panels )
    {
        const hsize_t np = static_cast<hsize_t>( n_panels );
        append_to_hdf5_dataset( _file, "/panels/pressure",
                                cusfloat_h5, { 1, np }, pressure );
        append_to_hdf5_dataset( _file, "/panels/phi_dt_term",
                                cusfloat_h5, { 1, np }, phi_dt_comp );
        append_to_hdf5_dataset( _file, "/panels/kinetic_term",
                                cusfloat_h5, { 1, np }, kinetic_comp );
        append_to_hdf5_dataset( _file, "/panels/hydrostatic_term",
                                cusfloat_h5, { 1, np }, hydrostatic_comp );
    }

    // ---- body kinematics + forces ----
    const int nb = static_cast<int>( body_pos.size( ) );
    for ( int ib = 0; ib < nb && ib < _n_bodies; ib++ )
    {
        std::ostringstream grp;
        grp << "/bodies/body_" << ib;

        append_to_hdf5_dataset( _file, grp.str( ) + "/position",
                                cusfloat_h5, { 1, 6 }, body_pos[ib].data( ) );
        append_to_hdf5_dataset( _file, grp.str( ) + "/velocity",
                                cusfloat_h5, { 1, 6 }, body_vel[ib].data( ) );
        append_to_hdf5_dataset( _file, grp.str( ) + "/acceleration",
                                cusfloat_h5, { 1, 6 }, body_acc[ib].data( ) );

        const std::string fgrp = grp.str( ) + "/forces";
        append_to_hdf5_dataset( _file, fgrp + "/total",
                                cusfloat_h5, { 1, 6 }, body_force_total[ib].data( ) );
        append_to_hdf5_dataset( _file, fgrp + "/phi_dt_term",
                                cusfloat_h5, { 1, 6 }, body_force_phi_dt[ib].data( ) );
        append_to_hdf5_dataset( _file, fgrp + "/radiation",
                                cusfloat_h5, { 1, 6 }, body_force_radiation[ib].data( ) );
        append_to_hdf5_dataset( _file, fgrp + "/wave_excitation",
                                cusfloat_h5, { 1, 6 }, body_force_wave[ib].data( ) );
        append_to_hdf5_dataset( _file, fgrp + "/kinetic_term",
                                cusfloat_h5, { 1, 6 }, body_force_kinetic[ib].data( ) );
        append_to_hdf5_dataset( _file, fgrp + "/hydrostatic_term",
                                cusfloat_h5, { 1, 6 }, body_force_hydrostatic[ib].data( ) );
    }

    // Make partial results survive an abnormal solver termination
    // (segfault, Ctrl-C, OOM, etc.). HDF5 chunk caches the per-step writes
    // and would otherwise only commit them on H5Fclose.
    H5Fflush( _file, H5F_SCOPE_GLOBAL );
}


/*****************************************************************************
 * close
 *****************************************************************************/

void TimeDomainHDF5Exporter::close( )
{
    if ( _file != H5I_INVALID_HID )
    {
        H5Fclose( _file );
        _file = H5I_INVALID_HID;
    }
}


/*****************************************************************************
 * reopen
 *****************************************************************************/

void TimeDomainHDF5Exporter::reopen( )
{
    if ( _file == H5I_INVALID_HID || _file_path.empty( ) ) { return; }

    // Closing the file releases all HDF5 caches and triggers the underlying
    // VFD's fsync/_commit, so every byte buffered by HDF5 reaches the OS and
    // (with the default sec2 driver) the disk.
    H5Fclose( _file );
    _file = H5I_INVALID_HID;

    // Reopen in read/write mode so subsequent append_step() calls keep working.
    _file = H5Fopen( _file_path.c_str( ), H5F_ACC_RDWR, H5P_DEFAULT );

    if ( _file == H5I_INVALID_HID )
    {
        std::cerr << "[ERROR] TimeDomainHDF5Exporter::reopen: H5Fopen failed for "
                  << _file_path << " — subsequent steps will be skipped." << std::endl;
    }
}
