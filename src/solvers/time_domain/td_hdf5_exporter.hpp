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

#pragma once

#include <array>
#include <string>
#include <vector>

#include "hdf5.h"
#include "../../config.hpp"

/**
 * @brief HDF5 time-series exporter for the time-domain BEM solver.
 *
 * Creates a single file  <folder_path>/1_results/results.td.h5  and appends
 * one record per simulated time step.  The file layout is:
 *
 * @code
 *   /time                          (n_steps,)          simulation time [s]
 *   /panels/
 *       pressure                   (n_steps, n_panels)  total Bernoulli pressure [Pa]
 *       phi_dt_term                (n_steps, n_panels)  -rho * dphi/dt  [Pa]
 *       kinetic_term               (n_steps, n_panels)  -rho * 0.5 * |grad phi|^2  [Pa]
 *       hydrostatic_term           (n_steps, n_panels)  -rho * g * z_center  [Pa]
 *   /bodies/
 *       body_0/
 *           position               (n_steps, 6)  CoG [x,y,z,rx,ry,rz]  [m, rad]
 *           velocity               (n_steps, 6)  [m/s, rad/s]
 *           acceleration           (n_steps, 6)  [m/s², rad/s²]
 *           forces/
 *               total              (n_steps, 6)  [Fx,Fy,Fz,Mx,My,Mz]  [N, N·m]
 *               phi_dt_term        (n_steps, 6)  -rho * dphi/dt  contribution
 *               kinetic_term       (n_steps, 6)  -rho * 0.5 * |grad phi|^2 contribution
 *               hydrostatic_term   (n_steps, 6)  -rho * g * z  contribution
 *       body_1/ ...
 * @endcode
 *
 * All datasets are extendable along the first (time) dimension and compressed
 * with deflate-4.  The float precision matches @c cusfloat (double by default,
 * float with -DSIMPLE_PREC).
 *
 * Usage:
 * @code
 *   TimeDomainHDF5Exporter exp;
 *   exp.initialize(case_folder, n_panels, n_bodies);
 *   // ... inside time loop:
 *   exp.append_step(t, pressure, phi_dt_comp, kinetic_comp, hydrostatic_comp,
 *                   n_panels, body_pos, body_vel, body_acc);
 *   // ... after loop:
 *   exp.close();
 * @endcode
 */
class TimeDomainHDF5Exporter
{
public:
    TimeDomainHDF5Exporter( ) = default;
    ~TimeDomainHDF5Exporter( );

    /**
     * @brief Open the HDF5 file and create all extendable datasets.
     *
     * @param folder_path  Case root directory (the file is created in
     *                     @c folder_path/1_results/results.td.h5 ).
     * @param n_panels     Number of BEM panels (fixed for the whole run).
     * @param n_bodies     Number of rigid bodies.
     */
    void initialize(
        const std::string&  folder_path,
        int                 n_panels,
        int                 n_bodies
    );

    /**
     * @brief Append one time step to all datasets.
     *
     * @param t                Simulation time [s].
     * @param pressure         Total Bernoulli pressure per panel [n_panels].
     * @param phi_dt_comp      -rho * dphi/dt per panel [n_panels].
     * @param kinetic_comp     -rho * 0.5 * |grad phi|^2 per panel [n_panels].
     * @param hydrostatic_comp -rho * g * z_center per panel [n_panels].
     * @param n_panels         Number of panels (must match the value given to initialize()).
     * @param body_pos         Body CoG positions  [n_bodies][6].
     * @param body_vel         Body CoG velocities [n_bodies][6].
     * @param body_acc         Body CoG accelerations [n_bodies][6].
     * @param body_force_total       Integrated hydrodynamic forces [n_bodies][6].
     * @param body_force_phi_dt      -rho*dphi/dt  force contribution [n_bodies][6].
     * @param body_force_kinetic     -rho*0.5*|∇φ|² force contribution [n_bodies][6].
     * @param body_force_hydrostatic -rho*g*z  force contribution [n_bodies][6].
     */
    void append_step(
        cusfloat                                            t,
        const cusfloat*                                     pressure,
        const cusfloat*                                     phi_dt_comp,
        const cusfloat*                                     kinetic_comp,
        const cusfloat*                                     hydrostatic_comp,
        int                                                 n_panels,
        const std::vector<std::array<cusfloat, 6>>&         body_pos,
        const std::vector<std::array<cusfloat, 6>>&         body_vel,
        const std::vector<std::array<cusfloat, 6>>&         body_acc,
        const std::vector<std::array<cusfloat, 6>>&         body_force_total,
        const std::vector<std::array<cusfloat, 6>>&         body_force_phi_dt,
        const std::vector<std::array<cusfloat, 6>>&         body_force_kinetic,
        const std::vector<std::array<cusfloat, 6>>&         body_force_hydrostatic
    );

    /** @brief Flush and close the HDF5 file.  Safe to call multiple times. */
    void close( );

    /** @return @c true when the file is open and ready to receive data. */
    bool is_open( ) const { return _file != H5I_INVALID_HID; }

private:
    hid_t   _file       = H5I_INVALID_HID;
    int     _n_panels   = 0;
    int     _n_bodies   = 0;
};
