
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
#include <string>
#include <vector>

// Include local modules
#include "../../config.hpp"
#include "../../containers/body_def.hpp"


/**
 * @brief Input structure for the time-domain solver.
 *
 * Reads the following case.input.dat format:
 *
 * @code{.text}
 *   -------- SeaMotions input file --------
 *   1.1.0           version
 *   ------- Solver Controls -------
 *   100.0           SimTime        - Total simulation time [s]
 *   0.01            Dt             - Time step
 *   true            UseDuhamel     - Enable Duhamel (radiation memory) integral
 *   false           UseGkDuhamel   - Use GK + sigma interpolation (true) or trapezoidal sum (false)
 *   10              GkNSeg         - Number of composite GK segments to divide HistTime into (>=1)
 *   10.0            RampTime       - Force ramp-up duration [s] (0 = disabled)
 *   0.0             HistTime       - Duhamel history window [s] (0 = full sim_time)
 *   100             OutFlushNIter  - Close+reopen the time-series HDF5 file every N output steps to force OS-level commit (0 = disabled)
 *   ------- Body Definition -------
 *   s11.bd.dat      BodyFN         - body file name(s)
 *   ------- Field Points ----------
 *   false           UseFP
 *   free_surface.fp.dat  FieldPN
 *   ------- Output Channels -------
 *   true            OutFK
 *   true            OutFH
 *   true            OutFT
 *   false           OutPot
 *   false           OutPress
 *   true            OutMesh
 *   true            OutSources
 *   true            OutStMass
 *   ------- Site Conditions -------
 *   1025.0          RhoW
 *   9.81            GravAcc
 *   50.0            WaterDepth
 *   1.0             WaveAmp        - Wave amplitude [m] (0 = no incident wave)
 *   3.0             WaveP          - Wave period [s] (0 = no incident wave)
 *   0               WaveH          - Wave heading [deg]
 * @endcode
 */
struct InputT
{
public:
    // Solver controls
    cusfloat                    sim_time            = 100.0;
    cusfloat                    dt                  = 0.01;
    bool                        use_duhamel         = true;    // enable radiation memory (Duhamel) integral
    bool                        use_gk_duhamel      = false;   // use GK + sigma interpolation instead of trapezoidal sum
    int                         gk_n_seg            = 10;      // number of composite GK macro-intervals over the full Duhamel history
    cusfloat                    wave_heading        = 0.0;     // heading in degrees

    // Body definitions
    std::vector<std::string>    bodies_finame       ;
    int                         bodies_np           = 0;
    BodyDef**                   bodies              = nullptr;

    // Field points
    bool                        use_field_points    = false;
    std::string                 field_points_finame = "";

    // Output channels
    bool                        out_fk              = true;
    bool                        out_fh              = true;
    bool                        out_ft              = true;
    bool                        out_potential       = false;
    bool                        out_pressure        = false;
    bool                        out_mesh            = true;
    bool                        out_sources         = true;
    bool                        out_struct_mass     = true;

    // Site conditions
    cusfloat                    water_density       = 1025.0;
    cusfloat                    grav_acc            = 9.81;
    cusfloat                    water_depth         = 50.0;
    cusfloat                    wave_amp            = 0.0;     // wave amplitude [m]  (0 = no incident wave)
    cusfloat                    wave_period         = 0.0;     // wave period [s]     (0 = no incident wave)
    cusfloat                    ang_freq            = 0.0;     // angular frequency [rad/s] derived from wave_period
    cusfloat                    ramp_time           = 0.0;     // force ramp-up duration [s] (0 = disabled)
    cusfloat                    duhamel_hist_time   = 0.0;     // Duhamel history window [s] (0 = full sim_time)
    int                         output_flush_niter  = 100;     // Close+reopen the time-series HDF5 file every N steps (0 = disabled)
    int                         added_mass_update_niter = 0;   // Recompute the infinite-frequency added mass every N steps (0 = frozen at init)

    // Paths
    std::string                 case_fopath         = "";
    std::string                 folder_path         = "";

    // DOFs
    const int                   dofs_np             = 6;

    // Constructor and destructor
    InputT( ) = default;
    ~InputT( );

    // Methods
    void    load( const std::string& folder_path );
    void    initialize( );

    static  void    read_body(
                                    const std::string& folder_path,
                                    const std::string& target_file,
                                    BodyDef* body
                              );

private:
    void    read_case( const std::string& folder_path );
    void    read_bodies( const std::string& folder_path );
};
