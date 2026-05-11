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
#include <cmath>
#include <iostream>
#include <limits>

// Include local modules
#include "math_tools.hpp"
#include "./sparse/sparse_math.hpp"
#include "generalized_alpha.hpp"

//--------------------------------------------------------------------
//-- Define NonLinear Generalized-Alpha solver
//--------------------------------------------------------------------

/**
 * @brief Nonlinear Generalized-Alpha time solver with Newton-Raphson iteration.
 *
 * Couples the Generalized-Alpha time-integration scheme with a full Newton-Raphson
 * iteration loop at each time step to handle nonlinear systems of the form:
 *
 *   M * a + f_int(u, v) = F_ext(t, u, v)
 *
 * where f_int may depend nonlinearly on displacement u and velocity v (e.g. geometric
 * nonlinearity, nonlinear damping, nonlinear restoring forces, Morison drag, mooring).
 *
 * Algorithm per time step:
 *   1. Predict kinematics at t_{n+1} using the Generalized-Alpha predictors.
 *   2. Newton-Raphson loop until ||R|| < tol:
 *      a. Evaluate residual:
 *           R = (1-alpha_m)*M*a_{n+1} + alpha_m*M*a_n
 *             + (1-alpha_f)*f_int(u_{n+1}, v_{n+1})
 *             + alpha_f*f_int(u_n, v_n)
 *             - F_ext(t_{n+1-alpha_f}, u_{n+1}, v_{n+1})
 *      b. Assemble effective tangent stiffness:
 *           K_eff = (1-alpha_f)*beta*dt^2 * K_T
 *                 + (1-alpha_f)*gamma*dt   * C_T
 *                 + (1-alpha_m)            * M
 *      c. Solve linear system:  K_eff * delta_u = -R
 *      d. Update kinematics (Newmark correctors):
 *           u_{n+1} += delta_u
 *           v_{n+1} += gamma/(beta*dt) * delta_u
 *           a_{n+1} += 1/(beta*dt^2)  * delta_u
 *   3. Apply restrictions.
 *   4. Advance time.
 *
 * The class is templated on two user-supplied functor types:
 *
 *   TResidual — computes the internal-force / residual contribution.
 *               Interface:
 *                 void operator()(
 *                     cusfloat    t,           // Current time t_{n+1-alpha_f}
 *                     cusfloat    dt,          // Time step
 *                     cusfloat*   u_pos,       // Current displacement trial u_{n+1}
 *                     cusfloat*   u_vel,       // Current velocity trial v_{n+1}
 *                     cusfloat*   u_acc,       // Current acceleration trial a_{n+1}
 *                     cusfloat*   u_pos_old,   // Displacement at t_n
 *                     cusfloat*   u_vel_old,   // Velocity at t_n
 *                     cusfloat*   u_acc_old,   // Acceleration at t_n
 *                     cusfloat*   residual     // Output: residual vector R (size rows_np)
 *                 );
 *
 *   TTangent — assembles the tangent stiffness K_T and tangent damping C_T.
 *              Interface:
 *                 void operator()(
 *                     cusfloat    t,           // Current time t_{n+1-alpha_f}
 *                     cusfloat    dt,          // Time step
 *                     cusfloat*   u_pos,       // Current displacement trial u_{n+1}
 *                     cusfloat*   u_vel,       // Current velocity trial v_{n+1}
 *                     cusfloat*   u_acc,       // Current acceleration trial a_{n+1}
 *                     CSRMatrix*& K_T,         // Output: tangent stiffness  (caller-owned)
 *                     CSRMatrix*& C_T          // Output: tangent damping    (caller-owned)
 *                 );
 *              The functor must allocate new CSRMatrix objects for K_T and C_T.
 *              GeneralizedAlphaNL will delete them after building K_eff.
 *
 * @tparam TResidual  Functor type for residual assembly (see above)
 * @tparam TTangent   Functor type for tangent assembly  (see above)
 * @tparam TLinear    Functor type forwarded to the inner GeneralizedAlpha for the
 *                    linear predictor step (can be a zero-force functor).
 */
template<typename TResidual, typename TTangent, typename TLinear>
class GeneralizedAlphaNL
{
private:
    /* Declare private class attributes */
    bool    _is_restrictions    = false;    // Switch to check if restrictions were internally allocated
    int*    _restrictions       = nullptr;  // Restriction flags (0 = free, 1 = fixed)

    /* Declare private class methods */

    /**
     * @brief Apply DOF restrictions to a vector (zero out fixed DOFs).
     * @param _vec  Vector to restrict (size rows_np)
     */
    void    _apply_restrictions( cusfloat* _vec );

    /**
     * @brief Delegate constructor — builds all internal state.
     */
    void    _build(
                        TResidual*      residual_fcn_in,
                        TTangent*       tangent_fcn_in,
                        TLinear*        linear_fcn_in,
                        CSRMatrix*      mass_mat_in,
                        cusfloat        time_step_in,
                        cusfloat        t0_in,
                        cusfloat        rho_inf_in,
                        int             max_iter_in,
                        cusfloat        abs_tol_in,
                        cusfloat        rel_tol_in,
                        cusfloat*       y0_pos_in,
                        cusfloat*       y0_vel_in,
                        cusfloat*       y0_acc_in,
                        int*            restrictions
                   );

public:
    // ----------------------------------------------------------------
    // Public attributes
    // ----------------------------------------------------------------
    cusfloat            alpha_f         = 0.0;      // Generalized-Alpha alpha_f parameter
    cusfloat            alpha_m         = 0.0;      // Generalized-Alpha alpha_m parameter
    cusfloat            beta            = 0.0;      // Newmark beta parameter
    int                 count           = 0;        // Step counter
    TResidual*          residual_fcn    = nullptr;  // Residual functor pointer
    TTangent*           tangent_fcn     = nullptr;  // Tangent functor pointer
    cusfloat            gamma           = 0.0;      // Newmark gamma parameter
    sparse_matrix_t*    mass_mat        = nullptr;  // Mass matrix (MKL handle)
    int                 rows_np         = 0;        // Number of DOFs
    cusfloat            rho_inf         = 0.0;      // Spectral radius at infinite frequency
    cusfloat            time            = 0.0;      // Current time
    cusfloat            time_init       = 0.0;      // Initial time
    cusfloat            time_old        = 0.0;      // Previous step time
    cusfloat            time_step       = 0.0;      // Time step dt
    cusfloat*           y0_acc          = nullptr;  // Initial acceleration (user-provided)
    cusfloat*           y0_pos          = nullptr;  // Initial position (user-provided)
    cusfloat*           y0_vel          = nullptr;  // Initial velocity (user-provided)
    cusfloat*           y_acc           = nullptr;  // Current acceleration a_{n+1}
    cusfloat*           y_pos           = nullptr;  // Current displacement u_{n+1}
    cusfloat*           y_vel           = nullptr;  // Current velocity v_{n+1}
    cusfloat*           y_acc_old       = nullptr;  // Previous acceleration a_n
    cusfloat*           y_pos_old       = nullptr;  // Previous displacement u_n
    cusfloat*           y_vel_old       = nullptr;  // Previous velocity v_n

    // Newton-Raphson convergence settings
    int                 max_iter        = 50;       // Maximum Newton-Raphson iterations
    cusfloat            abs_tol         = 1e-10;    // Absolute residual tolerance
    cusfloat            rel_tol         = 1e-8;     // Relative residual tolerance (||R||/||R0||)

    // Newton-Raphson diagnostics (last step)
    int                 nr_iter         = 0;        // Number of NR iterations in last step
    bool                nr_converged    = false;    // Convergence flag for last step
    cusfloat            nr_res_abs      = 0.0;      // Final absolute residual norm
    cusfloat            nr_res_rel      = 0.0;      // Final relative residual norm

    // ----------------------------------------------------------------
    // Constructors and destructor
    // ----------------------------------------------------------------

    /**
     * @brief Constructor without DOF restrictions.
     *
     * @param residual_fcn_in   Pointer to the residual functor
     * @param tangent_fcn_in    Pointer to the tangent functor
     * @param linear_fcn_in     Pointer to the linear (predictor) external force functor
     * @param mass_mat_in       Mass matrix in CSRMatrix format
     * @param time_step_in      Time step dt
     * @param t0_in             Start time
     * @param rho_inf_in        Spectral radius at infinite frequency ∈ [0, 1]
     * @param max_iter_in       Maximum Newton-Raphson iterations per step
     * @param abs_tol_in        Absolute residual convergence tolerance
     * @param rel_tol_in        Relative residual convergence tolerance
     * @param y0_pos_in         Initial displacement
     * @param y0_vel_in         Initial velocity
     * @param y0_acc_in         Initial acceleration (corrected internally)
     */
    GeneralizedAlphaNL(
                            TResidual*  residual_fcn_in,
                            TTangent*   tangent_fcn_in,
                            TLinear*    linear_fcn_in,
                            CSRMatrix*  mass_mat_in,
                            cusfloat    time_step_in,
                            cusfloat    t0_in,
                            cusfloat    rho_inf_in,
                            int         max_iter_in,
                            cusfloat    abs_tol_in,
                            cusfloat    rel_tol_in,
                            cusfloat*   y0_pos_in,
                            cusfloat*   y0_vel_in,
                            cusfloat*   y0_acc_in
                       );

    /**
     * @brief Constructor with DOF restrictions.
     *
     * @param residual_fcn_in   Pointer to the residual functor
     * @param tangent_fcn_in    Pointer to the tangent functor
     * @param linear_fcn_in     Pointer to the linear (predictor) external force functor
     * @param mass_mat_in       Mass matrix in CSRMatrix format
     * @param time_step_in      Time step dt
     * @param t0_in             Start time
     * @param rho_inf_in        Spectral radius at infinite frequency ∈ [0, 1]
     * @param max_iter_in       Maximum Newton-Raphson iterations per step
     * @param abs_tol_in        Absolute residual convergence tolerance
     * @param rel_tol_in        Relative residual convergence tolerance
     * @param y0_pos_in         Initial displacement
     * @param y0_vel_in         Initial velocity
     * @param y0_acc_in         Initial acceleration (corrected internally)
     * @param restrictions_in   Integer array: 0 = free DOF, 1 = fixed DOF
     */
    GeneralizedAlphaNL(
                            TResidual*  residual_fcn_in,
                            TTangent*   tangent_fcn_in,
                            TLinear*    linear_fcn_in,
                            CSRMatrix*  mass_mat_in,
                            cusfloat    time_step_in,
                            cusfloat    t0_in,
                            cusfloat    rho_inf_in,
                            int         max_iter_in,
                            cusfloat    abs_tol_in,
                            cusfloat    rel_tol_in,
                            cusfloat*   y0_pos_in,
                            cusfloat*   y0_vel_in,
                            cusfloat*   y0_acc_in,
                            int*        restrictions_in
                       );

    ~GeneralizedAlphaNL( void );

    // ----------------------------------------------------------------
    // Public methods
    // ----------------------------------------------------------------

    /**
     * @brief Retrieve linearly interpolated kinematics at an arbitrary time within the
     *        current step interval [time_old, time].
     *
     * @param time_at    Evaluation time (must satisfy time_old ≤ time_at ≤ time)
     * @param y_pos_at   Output: interpolated displacement
     * @param y_vel_at   Output: interpolated velocity
     * @param y_acc_at   Output: interpolated acceleration
     */
    void    get_values_at(
                                cusfloat    time_at,
                                cusfloat*   y_pos_at,
                                cusfloat*   y_vel_at,
                                cusfloat*   y_acc_at
                         );

    /**
     * @brief Initialize the solver: reset time counter, allocate working vectors,
     *        and correct initial acceleration via M*a0 = R(t0).
     */
    void    initialize( void );

    /**
     * @brief Linear interpolation of the solution within the last completed step.
     *
     * @param time       Evaluation time (must satisfy time_old ≤ time ≤ time)
     * @param y_pos_itp  Output: interpolated displacement
     * @param y_vel_itp  Output: interpolated velocity
     * @param y_acc_itp  Output: interpolated acceleration
     */
    void    interpolate_solution(
                                        cusfloat    time,
                                        cusfloat*   y_pos_itp,
                                        cusfloat*   y_vel_itp,
                                        cusfloat*   y_acc_itp
                                );

    /**
     * @brief Advance the solution by one time step using Generalized-Alpha + Newton-Raphson.
     *
     * On non-convergence the solver prints a warning but continues. Check
     * nr_converged after each step call if convergence monitoring is required.
     */
    void    step( void );

};


#include "generalized_alpha_nl.txx"
