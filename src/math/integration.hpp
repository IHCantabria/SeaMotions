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
#include <functional>

// Include local modules
#include "../config.hpp"
#include "../containers/panel_geom.hpp"
#include "../math/gauss.hpp"
#include "../math/gauss_kronrod_t.hpp"

// Generate alias for long type definitions in oder to clarify
// the function declaration
typedef std::function <cuscomplex(cuscomplex)> fcc_type;



#define GAUSS_2D_LOOP( storage, fcn )                               \
    storage += (                                                    \
                            GaussPointsT<2,NGP>::weights_x[i]       \
                            *                                       \
                            GaussPointsT<2,NGP>::weights_y[i]       \
                            *                                       \
                            fcn[i]                                  \
                            *                                       \
                            panel->jac_det_gauss_points[i]          \
                        );                                          \


/********************************************************/
/*************** Function Declaration *******************/
/********************************************************/
template<
            typename T, 
            typename U
        >                   inline  cuscomplex  adaptive_quadrature_panel(
                                                                                        T*              panel,
                                                                                        U               fcn,
                                                                                        cusfloat        abs_tol,
                                                                                        cusfloat        rel_tol,
                                                                                        bool            block_adaption=false,
                                                                                        bool            even_order=true,
                                                                                        int             go_fixed=1
                                                                                    );

                                            cuscomplex  complex_integration(
                                                                                        fcc_type        f_def,
                                                                                        cuscomplex      a,
                                                                                        cuscomplex      b,
                                                                                        cusfloat        tol
                                                                            );    
                                            
                                            cuscomplex  complex_integration_path(     
                                                                                        fcc_type        f_def,
                                                                                        int             num_segments,
                                                                                        cuscomplex*     waypoints,
                                                                                        cusfloat        tol,
                                                                                        bool            close_path,
                                                                                        bool            verbose
                                                                                );

template<
            int NGP, 
            typename T, 
            typename F, 
            typename L
        >                           inline void         gauss1d_loop( 
                                                                                        T&              storage, 
                                                                                        const F*        f, 
                                                                                        const L&        len 
                                                                    );

template<
            int NGP, 
            typename T, 
            typename F, 
            typename P
        >
                                    inline  void        gauss2d_loop(
                                                                                        T&              storage, 
                                                                                        F&&             f, 
                                                                                        const P*        panel
                                                                    );

template<
            int NGP, 
            typename T, 
            typename F, 
            typename P
        >
                                    inline  void        gauss2d_loop(
                                                                                        T&              storage, 
                                                                                        const F*        f, 
                                                                                        const P*        panel
                                                                    );

template<
            typename T, 
            typename U
        >                                   cuscomplex  quadrature_panel(
                                                                                        T*              panel,
                                                                                        U               target_fcn,
                                                                                        int             gp_order
                                                                        );

template<
            typename T, 
            typename U
        >                                   cuscomplex  quadrature_panel(
                                                                                        T*              panel,
                                                                                        U               target_fcn,
                                                                                        GaussPoints*    gp
                                                                    );

template<
            typename T, 
            typename U,
            auto Kernel,
            int NGP
        >                                   void        quadrature_panel_t(
                                                                                        T*                  panel,
                                                                                        U                   target_fcn,
                                                                                        cuscomplex&         result_G,
                                                                                        cuscomplex&         result_G_dn_sf,
                                                                                        cuscomplex&         result_G_dn_pf,
                                                                                        cuscomplex&         result_G_dx,
                                                                                        cuscomplex&         result_G_dy,
                                                                                        cuscomplex&         result_G_dz,
                                                                                        bool                verbose=false
                                                                        );

template<
            typename T, 
            typename U,
            int NGP,
            int NGPT
        >                                   void        quadrature_panel_time_t(
                                                                                        T*                  panel,
                                                                                        U&                  target_fcn,
                                                                                        cusfloat            t0,
                                                                                        cusfloat            t1,
                                                                                        cusfloat&           result_G_dtn,
                                                                                        cusfloat&           result_G_dtx,
                                                                                        cusfloat&           result_G_dty,
                                                                                        cusfloat&           result_G_dtz,
                                                                                        cusfloat&           result_G_dtt,
                                                                                        cusfloat&           result_G_dttn,
                                                                                        cusfloat&           result_G_dttx,
                                                                                        cusfloat&           result_G_dtty,
                                                                                        cusfloat&           result_G_dttz,
                                                                                        bool                verbose=false
                                                                        );

/**
 * @brief Like quadrature_panel_time_t but weights the spatial integral at each
 *        time quadrature point by sigma_at(t_q).
 *
 * The returned results already include the sigma(t) weighting, i.e.:
 *   result_G_dtn = ∫_{t0}^{t1} sigma(t) · [∫_panel dG/dn dA] dt
 *
 * This enables higher-order Duhamel integration when sigma varies within a
 * time sub-interval: pass a linear interpolator of sigma_hist as sigma_at.
 *
 * @tparam T        Panel geometry type.
 * @tparam U        Green's function interface type.
 * @tparam NGP      Number of spatial Gauss points per edge.
 * @tparam NGPT     Number of temporal Gauss points.
 * @tparam SigmaFcn Callable with signature cusfloat(cusfloat t_abs).
 */
template<
            typename T,
            typename U,
            int NGP,
            int NGPT,
            typename SigmaFcn
        >                                   void        quadrature_panel_time_t_sigma(
                                                                                        T*                  panel,
                                                                                        U&                  target_fcn,
                                                                                        cusfloat            t0,
                                                                                        cusfloat            t1,
                                                                                        SigmaFcn            sigma_at,
                                                                                        cusfloat&           result_G_dtn,
                                                                                        cusfloat&           result_G_dtx,
                                                                                        cusfloat&           result_G_dty,
                                                                                        cusfloat&           result_G_dtz,
                                                                                        cusfloat&           result_G_dtt,
                                                                                        cusfloat&           result_G_dttn,
                                                                                        cusfloat&           result_G_dttx,
                                                                                        cusfloat&           result_G_dtty,
                                                                                        cusfloat&           result_G_dttz,
                                                                                        bool                verbose=false
                                                                        );

template<
            typename Functor
        >                                   cusfloat    romberg_quadrature(   
                                                                                        Functor         f, 
                                                                                        cusfloat        a, 
                                                                                        cusfloat        b, 
                                                                                        cusfloat        precision
                                                                            );

/**
 * @brief Like quadrature_panel_time_t but uses a Gauss-Kronrod rule for the
 *        time integration instead of standard Gauss-Legendre.
 *
 * @tparam T     Panel geometry type.
 * @tparam U     Green's function interface type.
 * @tparam NGP   Number of spatial Gauss points per edge.
 * @tparam NGKPT Total number of Kronrod points (7, 15, or 21).
 */
template<
            typename T,
            typename U,
            int NGP,
            int NGKPT
        >                                   void        quadrature_panel_time_t_gk(
                                                                                        T*                  panel,
                                                                                        U&                  target_fcn,
                                                                                        cusfloat            t0,
                                                                                        cusfloat            t1,
                                                                                        cusfloat&           result_G_dtn,
                                                                                        cusfloat&           result_G_dtx,
                                                                                        cusfloat&           result_G_dty,
                                                                                        cusfloat&           result_G_dtz,
                                                                                        cusfloat&           result_G_dtt,
                                                                                        cusfloat&           result_G_dttn,
                                                                                        cusfloat&           result_G_dttx,
                                                                                        cusfloat&           result_G_dtty,
                                                                                        cusfloat&           result_G_dttz,
                                                                                        bool                verbose=false
                                                                        );

/**
 * @brief Like quadrature_panel_time_t_sigma but uses a Gauss-Kronrod rule for
 *        the time integration and weights each Kronrod point by sigma_at(t_k).
 *
 * The returned results already include the sigma(t) weighting:
 *   result_G_dtn = ∫_{t0}^{t1} sigma(t) · [∫_panel dG/dn dA] dt
 *
 * @tparam T        Panel geometry type.
 * @tparam U        Green's function interface type.
 * @tparam NGP      Number of spatial Gauss points per edge.
 * @tparam NGKPT    Total number of Kronrod points (7, 15, or 21).
 * @tparam SigmaFcn Callable with signature cusfloat(cusfloat t_abs).
 */
template<
            typename T,
            typename U,
            int NGP,
            int NGKPT,
            typename SigmaFcn
        >                                   void        quadrature_panel_time_t_sigma_gk(
                                                                                        T*                  panel,
                                                                                        U&                  target_fcn,
                                                                                        cusfloat            t0,
                                                                                        cusfloat            t1,
                                                                                        SigmaFcn            sigma_at,
                                                                                        cusfloat&           result_G_dtn,
                                                                                        cusfloat&           result_G_dtx,
                                                                                        cusfloat&           result_G_dty,
                                                                                        cusfloat&           result_G_dtz,
                                                                                        cusfloat&           result_G_dtt,
                                                                                        cusfloat&           result_G_dttn,
                                                                                        cusfloat&           result_G_dttx,
                                                                                        cusfloat&           result_G_dtty,
                                                                                        cusfloat&           result_G_dttz,
                                                                                        bool                verbose=false
                                                                        );

// Include templates definition
#include "integration.txx"
