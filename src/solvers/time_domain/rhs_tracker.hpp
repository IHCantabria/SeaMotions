/*
 * Copyright (c) 2025 Sergio Fernández Ruano / IHCantabria
 *
 * This file is part of SeaMotions Software.
 *
 * SeaMotions is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

#pragma once

// ---------------------------------------------------------------------------
// rhs_tracker.hpp
//
// Diagnostic CSV logger for the Duhamel convolution.  Three output files are
// written per MPI rank:
//
//   rhs_track_nodes_rank<R>.csv      one row per (step, obs, src, seg, kt)
//   rhs_track_segments_rank<R>.csv   one row per (step, obs, src, seg)
//   rhs_track_gp_rank<R>.csv         one row per (step, obs, src, seg, kt, gp)
//
// The "nodes" file breaks down the spatial-integrated kernel values at every
// time-quadrature node so you can see exactly which (lag, src panel,
// quadrature node) is responsible for a sudden change in _rhs_duhamel between
// consecutive steps.
//
// The "gp" file adds the underlying per-Gauss-point geometry (R, dX, dY, dZ,
// dZp, beta, mu) for every spatial Gauss point on the source panel at every
// tracked time-quadrature node.  This is useful when a step boundary is
// suspected of pushing a GP through a kernel singularity (e.g. a waterline
// crossing turning μ ≈ 0 or R ≈ 0).
//
// CONFIGURATION (edit the static sets in tracked_steps()/tracked_obs() and
// recompile).  When the current step or observer panel is not in the set,
// the tracking branch is bypassed and there is no run-time cost.
// ---------------------------------------------------------------------------

#include <cstdio>
#include <fstream>
#include <mutex>
#include <set>
#include <string>

#include "mpi.h"

#include "../../config.hpp"


namespace tdtrack {

class RhsTracker
{
public:
    static RhsTracker& instance()
    {
        static RhsTracker t;
        return t;
    }

    // ----- CONFIG ---------------------------------------------------------
    // Time steps to log.  step = round(t_current / dt).
    static const std::set<int>& tracked_steps()
    {
        static const std::set<int> s = { 157, 158, 159, 160 };
        return s;
    }
    // Observer-panel indices (the "obs" loop variable in the accumulator).
    static const std::set<int>& tracked_obs()
    {
        static const std::set<int> s = { 345, 346, 347, 348, 349, 350 };
        return s;
    }
    // ----------------------------------------------------------------------

    bool is_tracking_step( int step ) const
    {
        return tracked_steps().count( step ) > 0;
    }
    bool is_tracking_obs( int obs ) const
    {
        return tracked_obs().count( obs ) > 0;
    }

    void log_node(
                    int     step,
                    double  t_current,
                    int     obs,
                    int     src,
                    int     seg,
                    int     kt,
                    double  lag,
                    double  sigma_q,
                    double  tmp_dtn,
                    double  tmp_dtt,
                    double  tmp_dtx,
                    double  tmp_dty,
                    double  tmp_dtz,
                    double  tmp_dttn,
                    double  scale,
                    double  node_contrib_dtn,
                    double  mu_min,
                    double  mu_max,
                    double  beta_max,
                    double  R_min
                 )
    {
        open_files_if_needed();
        std::lock_guard<std::mutex> lk( _mu );
        _nodes_ofs
            << step          << ',' << t_current  << ','
            << obs           << ',' << src        << ',' << seg << ',' << kt << ','
            << lag           << ',' << sigma_q    << ','
            << tmp_dtn       << ',' << tmp_dtt    << ','
            << tmp_dtx       << ',' << tmp_dty    << ',' << tmp_dtz << ','
            << tmp_dttn      << ','
            << scale         << ',' << node_contrib_dtn << ','
            << mu_min        << ',' << mu_max     << ','
            << beta_max      << ',' << R_min
            << '\n';
    }

    void log_segment(
                        int     step,
                        double  t_current,
                        int     obs,
                        int     src,
                        int     seg,
                        double  dtn_val,
                        double  dttn_val,
                        double  contrib_rhs
                    )
    {
        open_files_if_needed();
        std::lock_guard<std::mutex> lk( _mu );
        _segs_ofs
            << step      << ',' << t_current << ','
            << obs       << ',' << src       << ',' << seg << ','
            << dtn_val   << ',' << dttn_val  << ',' << contrib_rhs
            << '\n';
    }

    // One row per spatial Gauss point at a tracked (step, obs, src, seg, kt)
    // time-quadrature node.  All quantities are in SI / kernel-native units.
    void log_gp(
                    int     step,
                    double  t_current,
                    int     obs,
                    int     src,
                    int     seg,
                    int     kt,
                    int     gp,
                    double  lag,
                    double  R,
                    double  dX,
                    double  dY,
                    double  dZ,
                    double  dZp,
                    double  beta,
                    double  mu
               )
    {
        open_files_if_needed();
        std::lock_guard<std::mutex> lk( _mu );
        _gp_ofs
            << step          << ',' << t_current  << ','
            << obs           << ',' << src        << ',' << seg << ',' << kt << ',' << gp << ','
            << lag           << ','
            << R             << ','
            << dX            << ',' << dY         << ',' << dZ  << ',' << dZp << ','
            << beta          << ',' << mu
            << '\n';
    }

    ~RhsTracker()
    {
        if ( _nodes_ofs.is_open() ) _nodes_ofs.close();
        if ( _segs_ofs .is_open() ) _segs_ofs .close();
        if ( _gp_ofs   .is_open() ) _gp_ofs   .close();
    }

private:
    RhsTracker() = default;

    void open_files_if_needed()
    {
        if ( _files_open ) return;
        std::lock_guard<std::mutex> lk( _mu );
        if ( _files_open ) return;

        int rank      = 0;
        int mpi_init  = 0;
        MPI_Initialized( &mpi_init );
        if ( mpi_init ) MPI_Comm_rank( MPI_COMM_WORLD, &rank );

        char nodes_name[256], segs_name[256], gp_name[256];
        std::snprintf( nodes_name, sizeof( nodes_name ),
                       "rhs_track_nodes_rank%d.csv",    rank );
        std::snprintf( segs_name,  sizeof( segs_name  ),
                       "rhs_track_segments_rank%d.csv", rank );
        std::snprintf( gp_name,    sizeof( gp_name    ),
                       "rhs_track_gp_rank%d.csv",       rank );

        _nodes_ofs.open( nodes_name, std::ios::out );
        _segs_ofs .open( segs_name,  std::ios::out );
        _gp_ofs   .open( gp_name,    std::ios::out );

        _nodes_ofs.setf( std::ios::scientific );
        _nodes_ofs.precision( 10 );
        _segs_ofs .setf( std::ios::scientific );
        _segs_ofs .precision( 10 );
        _gp_ofs   .setf( std::ios::scientific );
        _gp_ofs   .precision( 10 );

        _nodes_ofs
            << "step,t_current,obs,src,seg,kt,lag,sigma,"
               "tmp_dtn,tmp_dtt,tmp_dtx,tmp_dty,tmp_dtz,tmp_dttn,"
               "scale,node_contrib_dtn,mu_min,mu_max,beta_max,R_min\n";
        _segs_ofs
            << "step,t_current,obs,src,seg,dtn_val,dttn_val,contrib_rhs\n";
        _gp_ofs
            << "step,t_current,obs,src,seg,kt,gp,lag,R,dX,dY,dZ,dZp,beta,mu\n";

        _files_open = true;
    }

    std::ofstream _nodes_ofs;
    std::ofstream _segs_ofs;
    std::ofstream _gp_ofs;
    std::mutex    _mu;
    bool          _files_open = false;
};

}   // namespace tdtrack
