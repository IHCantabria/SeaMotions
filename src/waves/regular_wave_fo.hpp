
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
#include "../config.hpp"
#include "../math/math_constants.hpp"
#include "waves_common.hpp"
#include "wave_dispersion_base_fo.hpp"


/**
 * @brief First-order regular incident wave — single source of truth.
 *
 * Describes the linear incident wave used by the time-domain solver and owns
 * everything derived from it: the dispersion wavenumber `k`, the heading `mu`
 * (radians) and the soft-start amplitude ramp.  Both consumers ask this one
 * object instead of re-deriving the wave locally:
 *
 *   * the BEM kernel — for the Froude-Krylov potential time/space derivatives
 *     (@ref phi_dt / @ref phi_dx / @ref phi_dy / @ref phi_dz), and
 *   * the body diffraction boundary condition in TimeSolver — for the incident
 *     normal velocity / acceleration at each panel
 *     (@ref normal_velocity / @ref normal_acceleration).
 *
 * The soft-start ramp is applied to the amplitude, so — because the underlying
 * first-order wave kernels are linear in amplitude — the Froude-Krylov pressure
 * AND the diffraction boundary condition (hence the scattered part of σ and its
 * Duhamel memory) ramp from the *same* definition.  This keeps the whole wave
 * excitation soft-starting consistently and removes the previous duplication of
 * the wave parameters and the ramp across two files (single-responsibility).
 *
 * Default-constructed instances are inactive (no wave), so a body run without
 * waves simply gets zero contributions.
 */
class RegularWaveFO
{
public:
    /// Inactive wave (no incident excitation).
    RegularWaveFO( ) = default;

    /**
     * @param amp          Wave amplitude [m].
     * @param ang_freq     Angular frequency w [rad/s].
     * @param heading_deg  Wave heading [deg] (converted to radians internally).
     * @param water_depth  Water depth h [m].
     * @param grav_acc     Gravitational acceleration g [m/s^2].
     * @param ramp_time    Soft-start duration [s]; <= 0 disables the ramp.
     */
    RegularWaveFO(
                    cusfloat amp,
                    cusfloat ang_freq,
                    cusfloat heading_deg,
                    cusfloat water_depth,
                    cusfloat grav_acc,
                    cusfloat ramp_time
                 )
        : _aw( amp )
        , _w( ang_freq )
        , _h( water_depth )
        , _g( grav_acc )
        , _ramp_time( ramp_time )
    {
        this->_mu     = heading_deg * PI / static_cast<cusfloat>( 180.0 );
        this->_active = ( _aw > static_cast<cusfloat>( 0.0 )
                          && _w > static_cast<cusfloat>( 0.0 ) );
        this->_k      = this->_active ? w2k( _w, _h, _g )
                                      : static_cast<cusfloat>( 0.0 );
    }

    /// True when a non-trivial incident wave is present.
    bool active( ) const { return this->_active; }

    /// Linear soft-start factor in [0, 1]; 1 when the ramp is disabled.
    cusfloat ramp( cusfloat t ) const
    {
        return ( this->_ramp_time > static_cast<cusfloat>( 0.0 ) && t < this->_ramp_time )
               ? ( t / this->_ramp_time )
               : static_cast<cusfloat>( 1.0 );
    }

    /// Ramped wave amplitude at time t.
    cusfloat amplitude( cusfloat t ) const { return this->_aw * this->ramp( t ); }

    // -- Accessors (raw, un-ramped wave parameters) -----------------------
    cusfloat ang_freq( )    const { return this->_w; }
    cusfloat wavenumber( )  const { return this->_k; }
    cusfloat heading_rad( ) const { return this->_mu; }

    // ---------------------------------------------------------------------
    // Froude-Krylov incident potential derivatives at a field point (ramped).
    // ---------------------------------------------------------------------
    cusfloat phi_dt( cusfloat x, cusfloat y, cusfloat z, cusfloat t ) const
    {
        if ( !this->_active ) { return static_cast<cusfloat>( 0.0 ); }
        return wave_potential_fo_time_dt( this->amplitude( t ), _w, _k, _h, _g, x, y, z, _mu, t );
    }
    cusfloat phi_dx( cusfloat x, cusfloat y, cusfloat z, cusfloat t ) const
    {
        if ( !this->_active ) { return static_cast<cusfloat>( 0.0 ); }
        return wave_potential_fo_time_dx( this->amplitude( t ), _w, _k, _h, _g, x, y, z, _mu, t );
    }
    cusfloat phi_dy( cusfloat x, cusfloat y, cusfloat z, cusfloat t ) const
    {
        if ( !this->_active ) { return static_cast<cusfloat>( 0.0 ); }
        return wave_potential_fo_time_dy( this->amplitude( t ), _w, _k, _h, _g, x, y, z, _mu, t );
    }
    cusfloat phi_dz( cusfloat x, cusfloat y, cusfloat z, cusfloat t ) const
    {
        if ( !this->_active ) { return static_cast<cusfloat>( 0.0 ); }
        return wave_potential_fo_time_dz( this->amplitude( t ), _w, _k, _h, _g, x, y, z, _mu, t );
    }

    // ---------------------------------------------------------------------
    // Diffraction boundary-condition helpers at a panel (ramped).
    //   n is the 3-component physical normal (normal_vec[0..2]).
    // ---------------------------------------------------------------------

    /// Incident-wave normal fluid velocity  v·n  (the diffraction BC magnitude).
    cusfloat normal_velocity(
                                cusfloat        x,
                                cusfloat        y,
                                cusfloat        z,
                                const cusfloat* n,
                                cusfloat        t
                            ) const
    {
        if ( !this->_active ) { return static_cast<cusfloat>( 0.0 ); }
        const cusfloat a  = this->amplitude( t );
        const cusfloat vx = wave_potential_fo_time_dx( a, _w, _k, _h, _g, x, y, z, _mu, t );
        const cusfloat vy = wave_potential_fo_time_dy( a, _w, _k, _h, _g, x, y, z, _mu, t );
        const cusfloat vz = wave_potential_fo_time_dz( a, _w, _k, _h, _g, x, y, z, _mu, t );
        return vx * n[0] + vy * n[1] + vz * n[2];
    }

    /// Incident-wave normal fluid acceleration  a·n  (∂t of the normal velocity).
    cusfloat normal_acceleration(
                                    cusfloat        x,
                                    cusfloat        y,
                                    cusfloat        z,
                                    const cusfloat* n,
                                    cusfloat        t
                                ) const
    {
        if ( !this->_active ) { return static_cast<cusfloat>( 0.0 ); }
        const cusfloat a   = this->amplitude( t );
        const cusfloat ax  = wave_potential_fo_time_dtdx( a, _w, _k, _h, _g, x, y, z, _mu, t );
        const cusfloat ay  = wave_potential_fo_time_dtdy( a, _w, _k, _h, _g, x, y, z, _mu, t );
        const cusfloat az  = wave_potential_fo_time_dtdz( a, _w, _k, _h, _g, x, y, z, _mu, t );
        return ax * n[0] + ay * n[1] + az * n[2];
    }

private:
    cusfloat _aw        = static_cast<cusfloat>( 0.0 );   ///< amplitude [m]
    cusfloat _w         = static_cast<cusfloat>( 0.0 );   ///< angular frequency [rad/s]
    cusfloat _k         = static_cast<cusfloat>( 0.0 );   ///< wavenumber [1/m]
    cusfloat _h         = static_cast<cusfloat>( 0.0 );   ///< water depth [m]
    cusfloat _g         = static_cast<cusfloat>( 0.0 );   ///< gravity [m/s^2]
    cusfloat _mu        = static_cast<cusfloat>( 0.0 );   ///< heading [rad]
    cusfloat _ramp_time = static_cast<cusfloat>( 0.0 );   ///< soft-start [s]
    bool     _active    = false;
};
