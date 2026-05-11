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

// Include general usage libraries
#include <iomanip>

// Include local modules
#include "generalized_alpha_nl.hpp"


//--------------------------------------------------------------------
//-- Private methods
//--------------------------------------------------------------------

template<typename TResidual, typename TTangent, typename TLinear>
void GeneralizedAlphaNL<TResidual, TTangent, TLinear>::_apply_restrictions(
                                                                                cusfloat* _vec
                                                                            )
{
    for ( int i=0; i<this->rows_np; i++ )
    {
        if ( this->_restrictions[i] > 0 )
        {
            _vec[i] = 0.0;
        }
    }
}


template<typename TResidual, typename TTangent, typename TLinear>
void GeneralizedAlphaNL<TResidual, TTangent, TLinear>::_build(
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
                                                                    int*        restrictions
                                                                )
{
    // linear_fcn_in is accepted for interface consistency but not stored (NL solver
    // uses the residual functor for all force evaluations).
    (void)linear_fcn_in;

    // Save scalar inputs
    this->residual_fcn  = residual_fcn_in;
    this->tangent_fcn   = tangent_fcn_in;
    this->rows_np       = mass_mat_in->rows_np;
    this->time_step     = time_step_in;
    this->time_init     = t0_in;
    this->rho_inf       = rho_inf_in;
    this->max_iter      = max_iter_in;
    this->abs_tol       = abs_tol_in;
    this->rel_tol       = rel_tol_in;
    this->y0_pos        = y0_pos_in;
    this->y0_vel        = y0_vel_in;
    this->y0_acc        = y0_acc_in;

    // Compute Generalized-Alpha parameters from rho_inf
    this->alpha_m   = ( 2.0*rho_inf_in - 1.0 ) / ( rho_inf_in + 1.0 );
    this->alpha_f   = rho_inf_in / ( rho_inf_in + 1.0 );
    this->gamma     = 0.5 - this->alpha_m + this->alpha_f;
    this->beta      = pow2s( 1.0 - this->alpha_m + this->alpha_f ) / 4.0;

    // Convert mass matrix from CSRMatrix to MKL sparse handle
    this->mass_mat = new sparse_matrix_t;
    CALL_AND_CHECK_STATUS(
        mkl_sparse_d_create_csr(
                                    this->mass_mat,
                                    SPARSE_INDEX_BASE_ZERO,
                                    mass_mat_in->rows_np,
                                    mass_mat_in->rows_np,
                                    mass_mat_in->row_index_cum,
                                    mass_mat_in->row_index_cum+1,
                                    mass_mat_in->col_index,
                                    mass_mat_in->values
                                ),
        "Error after MKL_SPARSE_D_CREATE_CSR - mass_mat (GeneralizedAlphaNL) \n"
    );

    // Store restrictions pointer if not already set by constructor
    if ( !this->_is_restrictions )
    {
        this->_restrictions = restrictions;
    }

    // Initialize solver state
    this->initialize();
}


//--------------------------------------------------------------------
//-- Public methods
//--------------------------------------------------------------------

template<typename TResidual, typename TTangent, typename TLinear>
void GeneralizedAlphaNL<TResidual, TTangent, TLinear>::initialize( void )
{
    struct matrix_descr descrA;
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;

    //////////////////////////////////////////////////////////////////
    ///////////// CORRECT INITIAL ACCELERATION  //////////////////////
    //////////////////////////////////////////////////////////////////
    // Evaluate the residual at t0 with zero initial guess for a0
    // and solve M*a0 = -R(u0, v0, 0, t0) to get a consistent a0.
    // The residual functor is expected to return:
    //   R = f_int(u0, v0) - F_ext(t0, u0, v0)
    // so that M*a0 = F_ext - f_int  →  a0 = M^{-1}*(F_ext - f_int)

    CSRMatrix*      mass_csr    = convert_mkl_to_csrmatrix( this->mass_mat );
    SparseSolver    acc_crr_sps = SparseSolver( mass_csr, false );

    cusfloat* acc_crr_rhs = generate_empty_vector<cusfloat>( this->rows_np );

    // Evaluate residual: R = f_int(u0,v0) - F_ext(t0,u0,v0)
    ( *this->residual_fcn )(
                                this->time_init,
                                this->time_step,
                                this->y0_pos,
                                this->y0_vel,
                                this->y0_acc,
                                this->y0_pos,
                                this->y0_vel,
                                this->y0_acc,
                                acc_crr_rhs
                            );

    // Negate to get the rhs: M*a0 = -R  →  a0 = M^{-1}*(-R) = M^{-1}*(F_ext - f_int)
    for ( int i=0; i<this->rows_np; i++ )
    {
        acc_crr_rhs[i] = -acc_crr_rhs[i];
    }

    // Solve for corrected initial acceleration
    clear_vector( this->rows_np, this->y0_acc );
    acc_crr_sps.solve( acc_crr_rhs, this->y0_acc );

    //////////////////////////////////////////////////////////////////
    /////////////// RESET INTEGRATION VARIABLES //////////////////////
    //////////////////////////////////////////////////////////////////
    this->time      = this->time_init;
    this->time_old  = this->time_init;
    this->count     = 0;
    this->y_pos     = generate_empty_vector<cusfloat>( this->rows_np );
    this->y_vel     = generate_empty_vector<cusfloat>( this->rows_np );
    this->y_acc     = generate_empty_vector<cusfloat>( this->rows_np );
    this->y_pos_old = generate_empty_vector<cusfloat>( this->rows_np );
    this->y_vel_old = generate_empty_vector<cusfloat>( this->rows_np );
    this->y_acc_old = generate_empty_vector<cusfloat>( this->rows_np );

    copy_vector( this->rows_np, this->y0_pos, this->y_pos );
    copy_vector( this->rows_np, this->y0_vel, this->y_vel );
    copy_vector( this->rows_np, this->y0_acc, this->y_acc );

    mkl_free( acc_crr_rhs );
}


template<typename TResidual, typename TTangent, typename TLinear>
void GeneralizedAlphaNL<TResidual, TTangent, TLinear>::get_values_at(
                                                                            cusfloat    time_at,
                                                                            cusfloat*   y_pos_at,
                                                                            cusfloat*   y_vel_at,
                                                                            cusfloat*   y_acc_at
                                                                        )
{
    cusfloat alpha = ( time_at - this->time_old ) / ( this->time - this->time_old );

    for ( int i=0; i<this->rows_np; i++ )
    {
        y_pos_at[i] = ( 1.0 - alpha ) * this->y_pos_old[i] + alpha * this->y_pos[i];
        y_vel_at[i] = ( 1.0 - alpha ) * this->y_vel_old[i] + alpha * this->y_vel[i];
        y_acc_at[i] = ( 1.0 - alpha ) * this->y_acc_old[i] + alpha * this->y_acc[i];
    }
}


template<typename TResidual, typename TTangent, typename TLinear>
void GeneralizedAlphaNL<TResidual, TTangent, TLinear>::interpolate_solution(
                                                                                    cusfloat    time,
                                                                                    cusfloat*   y_pos_itp,
                                                                                    cusfloat*   y_vel_itp,
                                                                                    cusfloat*   y_acc_itp
                                                                                )
{
    constexpr cusfloat TIME_PRECISION = 1e-6;
    if ( time < this->time_old - TIME_PRECISION )
    {
        std::cerr << std::setprecision(16) << std::endl;
        std::cerr << "GeneralizedAlphaNL::interpolate_solution - ";
        std::cerr << "Interpolation time is below time_old value: ";
        std::cerr << "time_old [s] = " << this->time_old << " - time_itp = " << time << std::endl;
        exit(300);
    }
    if ( time > this->time + TIME_PRECISION )
    {
        std::cerr << std::setprecision(16) << std::endl;
        std::cerr << "GeneralizedAlphaNL::interpolate_solution - ";
        std::cerr << "Interpolation time is above time value: ";
        std::cerr << "time [s] = " << this->time << " - time_itp = " << time << std::endl;
        exit(300);
    }

    cusfloat alpha = ( time - this->time_old ) / this->time_step;

    clear_vector( this->rows_np, y_pos_itp );
    clear_vector( this->rows_np, y_vel_itp );
    clear_vector( this->rows_np, y_acc_itp );

    sv_add( this->rows_np, (1.0-alpha), this->y_pos_old, alpha, this->y_pos, y_pos_itp );
    sv_add( this->rows_np, (1.0-alpha), this->y_vel_old, alpha, this->y_vel, y_vel_itp );
    sv_add( this->rows_np, (1.0-alpha), this->y_acc_old, alpha, this->y_acc, y_acc_itp );
}


//--------------------------------------------------------------------
//-- Constructors and destructor
//--------------------------------------------------------------------

template<typename TResidual, typename TTangent, typename TLinear>
GeneralizedAlphaNL<TResidual, TTangent, TLinear>::GeneralizedAlphaNL(
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
                                                                        )
{
    this->_restrictions     = generate_empty_vector<int>( mass_mat_in->rows_np );
    this->_is_restrictions  = true;

    this->_build(
                    residual_fcn_in,
                    tangent_fcn_in,
                    linear_fcn_in,
                    mass_mat_in,
                    time_step_in,
                    t0_in,
                    rho_inf_in,
                    max_iter_in,
                    abs_tol_in,
                    rel_tol_in,
                    y0_pos_in,
                    y0_vel_in,
                    y0_acc_in,
                    this->_restrictions
                );
}


template<typename TResidual, typename TTangent, typename TLinear>
GeneralizedAlphaNL<TResidual, TTangent, TLinear>::GeneralizedAlphaNL(
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
                                                                        )
{
    this->_build(
                    residual_fcn_in,
                    tangent_fcn_in,
                    linear_fcn_in,
                    mass_mat_in,
                    time_step_in,
                    t0_in,
                    rho_inf_in,
                    max_iter_in,
                    abs_tol_in,
                    rel_tol_in,
                    y0_pos_in,
                    y0_vel_in,
                    y0_acc_in,
                    restrictions_in
                );
}


template<typename TResidual, typename TTangent, typename TLinear>
GeneralizedAlphaNL<TResidual, TTangent, TLinear>::~GeneralizedAlphaNL( void )
{
    // Delete internally allocated restrictions vector
    if ( this->_is_restrictions )
    {
        mkl_free( this->_restrictions );
    }

    // Destroy mass matrix MKL handle
    mkl_sparse_destroy( *( this->mass_mat ) );
    delete this->mass_mat;

    // Free kinematics vectors
    mkl_free( this->y_pos     );
    mkl_free( this->y_vel     );
    mkl_free( this->y_acc     );
    mkl_free( this->y_pos_old );
    mkl_free( this->y_vel_old );
    mkl_free( this->y_acc_old );
}


//--------------------------------------------------------------------
//-- Step: Generalized-Alpha + Newton-Raphson
//--------------------------------------------------------------------

template<typename TResidual, typename TTangent, typename TLinear>
void GeneralizedAlphaNL<TResidual, TTangent, TLinear>::step( void )
{
    struct matrix_descr descrA;
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;

    //////////////////////////////////////////////////////////////////
    ///// COPY PREVIOUS STEP DATA INTO THE OLD VECTOR HANDLE /////////
    //////////////////////////////////////////////////////////////////
    copy_vector( this->rows_np, this->y_pos, this->y_pos_old );
    copy_vector( this->rows_np, this->y_vel, this->y_vel_old );
    copy_vector( this->rows_np, this->y_acc, this->y_acc_old );

    //////////////////////////////////////////////////////////////////
    ////////////////// GENERALIZED-ALPHA PREDICTORS //////////////////
    //////////////////////////////////////////////////////////////////
    // Displacement predictor:  u_{n+1}^(0) = u_n + dt*v_n + dt^2*(0.5-beta)*a_n
    // Velocity predictor:      v_{n+1}^(0) = v_n + dt*(1-gamma)*a_n
    // Acceleration predictor:  a_{n+1}^(0) = 0
    cusfloat* u_pred = generate_empty_vector<cusfloat>( this->rows_np );
    cusfloat* v_pred = generate_empty_vector<cusfloat>( this->rows_np );
    cusfloat* a_pred = generate_empty_vector<cusfloat>( this->rows_np );

    // u_pred = u_n + dt*v_n
    sv_add(
                this->rows_np,
                1.0, this->y_pos_old,
                this->time_step, this->y_vel_old,
                u_pred
            );
    // u_pred += dt^2*(0.5-beta)*a_n
    sv_add(
                this->rows_np,
                1.0, u_pred,
                pow2s( this->time_step ) * ( 0.5 - this->beta ), this->y_acc_old,
                u_pred
            );

    // v_pred = v_n + dt*(1-gamma)*a_n
    sv_add(
                this->rows_np,
                1.0, this->y_vel_old,
                this->time_step * ( 1.0 - this->gamma ), this->y_acc_old,
                v_pred
            );

    // a_pred = 0  (already zero from generate_empty_vector)

    // Copy predictors into the trial kinematics that will be updated by Newton-Raphson
    copy_vector( this->rows_np, u_pred, this->y_pos );
    copy_vector( this->rows_np, v_pred, this->y_vel );
    copy_vector( this->rows_np, a_pred, this->y_acc );

    //////////////////////////////////////////////////////////////////
    /////////////////// NEWTON-RAPHSON LOOP //////////////////////////
    //////////////////////////////////////////////////////////////////
    // Alpha-f interpolated time: t_{n+1-alpha_f} = t_n + (1-alpha_f)*dt
    cusfloat t_alpha_f = this->time + ( 1.0 - this->alpha_f ) * this->time_step;

    cusfloat* residual      = generate_empty_vector<cusfloat>( this->rows_np );
    cusfloat* delta_u       = generate_empty_vector<cusfloat>( this->rows_np );

    cusfloat res_norm_0     = std::numeric_limits<cusfloat>::max();
    cusfloat res_norm       = std::numeric_limits<cusfloat>::max();

    this->nr_converged      = false;
    this->nr_iter           = 0;

    for ( int iter=0; iter<this->max_iter; iter++ )
    {
        //------------------------------------------------------------
        // 1. Evaluate residual R at the alpha-f level:
        //    R = f_int(u_{n+1}, v_{n+1}, a_{n+1}) - F_ext(t_{n+1-alpha_f}, ...)
        //    (internally the functor must handle alpha blending if needed,
        //     or the residual can be the full nonlinear equation evaluated at
        //     the trial state — both interpretations are supported)
        //------------------------------------------------------------
        clear_vector( this->rows_np, residual );

        ( *this->residual_fcn )(
                                    t_alpha_f,
                                    this->time_step,
                                    this->y_pos,        // u_{n+1} trial
                                    this->y_vel,        // v_{n+1} trial
                                    this->y_acc,        // a_{n+1} trial
                                    this->y_pos_old,    // u_n
                                    this->y_vel_old,    // v_n
                                    this->y_acc_old,    // a_n
                                    residual
                                );

        // Add mass inertia contribution to the residual:
        // R += (1-alpha_m)*M*a_{n+1} + alpha_m*M*a_n
        // -> R += (1-alpha_m)*M*a_{n+1}
        CALL_AND_CHECK_STATUS(
            mkl_sparse_d_mv(
                                SPARSE_OPERATION_NON_TRANSPOSE,
                                ( 1.0 - this->alpha_m ),
                                *(this->mass_mat),
                                descrA,
                                this->y_acc,
                                1.0,
                                residual
                            ),
            "Error after MKL_SPARSE_D_MV - R += (1-alpha_m)*M*a_{n+1} \n"
        );
        // -> R += alpha_m*M*a_n
        CALL_AND_CHECK_STATUS(
            mkl_sparse_d_mv(
                                SPARSE_OPERATION_NON_TRANSPOSE,
                                this->alpha_m,
                                *(this->mass_mat),
                                descrA,
                                this->y_acc_old,
                                1.0,
                                residual
                            ),
            "Error after MKL_SPARSE_D_MV - R += alpha_m*M*a_n \n"
        );

        // Apply restrictions to residual
        this->_apply_restrictions( residual );

        //------------------------------------------------------------
        // 2. Check convergence
        //------------------------------------------------------------
        cusfloat _res_sq = 0.0;
        for ( int i=0; i<this->rows_np; i++ ) { _res_sq += residual[i] * residual[i]; }
        res_norm = std::sqrt( _res_sq );

        if ( iter == 0 )
        {
            res_norm_0 = ( res_norm > 0.0 ) ? res_norm : 1.0;
        }

        cusfloat res_rel = res_norm / res_norm_0;

        if ( res_norm < this->abs_tol || res_rel < this->rel_tol )
        {
            this->nr_converged  = true;
            this->nr_iter       = iter;
            this->nr_res_abs    = res_norm;
            this->nr_res_rel    = res_rel;
            break;
        }

        //------------------------------------------------------------
        // 3. Assemble effective tangent stiffness:
        //    K_eff = (1-alpha_f)*beta*dt^2 * K_T
        //          + (1-alpha_f)*gamma*dt   * C_T
        //          + (1-alpha_m)            * M
        //------------------------------------------------------------
        CSRMatrix* K_T  = nullptr;
        CSRMatrix* C_T  = nullptr;

        ( *this->tangent_fcn )(
                                    t_alpha_f,
                                    this->time_step,
                                    this->y_pos,
                                    this->y_vel,
                                    this->y_acc,
                                    K_T,
                                    C_T
                                );

        // Build K_eff from an empty sparse matrix, then accumulate each contribution
        // K_eff = (1-alpha_m)*M + (1-alpha_f)*gamma*dt*C_T + (1-alpha_f)*beta*dt^2*K_T
        sparse_matrix_t K_eff_mkl;
        sparse_matrix_t K_T_mkl;
        sparse_matrix_t C_T_mkl;

        CALL_AND_CHECK_STATUS(
            mkl_sparse_d_create_csr(
                                        &K_T_mkl,
                                        SPARSE_INDEX_BASE_ZERO,
                                        K_T->rows_np,
                                        K_T->rows_np,
                                        K_T->row_index_cum,
                                        K_T->row_index_cum+1,
                                        K_T->col_index,
                                        K_T->values
                                    ),
            "Error after MKL_SPARSE_D_CREATE_CSR - K_T_mkl \n"
        );

        CALL_AND_CHECK_STATUS(
            mkl_sparse_d_create_csr(
                                        &C_T_mkl,
                                        SPARSE_INDEX_BASE_ZERO,
                                        C_T->rows_np,
                                        C_T->rows_np,
                                        C_T->row_index_cum,
                                        C_T->row_index_cum+1,
                                        C_T->col_index,
                                        C_T->values
                                    ),
            "Error after MKL_SPARSE_D_CREATE_CSR - C_T_mkl \n"
        );

        // Create empty starting matrix for K_eff
        int*        keff_row_cum_dum            = generate_empty_vector<int>( this->rows_np+1 );
        int         keff_col_index_dum[1]       = {0};
        cusfloat    keff_values_dum[1]          = {0.0};
        keff_row_cum_dum[this->rows_np]         = 1;

        CALL_AND_CHECK_STATUS(
            mkl_sparse_d_create_csr(
                                        &K_eff_mkl,
                                        SPARSE_INDEX_BASE_ZERO,
                                        this->rows_np,
                                        this->rows_np,
                                        keff_row_cum_dum,
                                        keff_row_cum_dum+1,
                                        keff_col_index_dum,
                                        keff_values_dum
                                    ),
            "Error after MKL_SPARSE_D_CREATE_CSR - K_eff_mkl (empty) \n"
        );

        // K_eff += (1-alpha_m)*M
        // MATRIX CONTENTS -> K_eff = (1-alpha_m)*M
        CALL_AND_CHECK_STATUS(
            mkl_sparse_d_add(
                                SPARSE_OPERATION_NON_TRANSPOSE,
                                *(this->mass_mat),
                                ( 1.0 - this->alpha_m ),
                                K_eff_mkl,
                                &K_eff_mkl
                            ),
            "Error after MKL_SPARSE_D_ADD - K_eff += (1-alpha_m)*M \n"
        );

        // K_eff += (1-alpha_f)*gamma*dt * C_T
        // MATRIX CONTENTS -> K_eff = (1-alpha_m)*M + (1-alpha_f)*gamma*dt*C_T
        CALL_AND_CHECK_STATUS(
            mkl_sparse_d_add(
                                SPARSE_OPERATION_NON_TRANSPOSE,
                                C_T_mkl,
                                ( 1.0 - this->alpha_f ) * this->gamma * this->time_step,
                                K_eff_mkl,
                                &K_eff_mkl
                            ),
            "Error after MKL_SPARSE_D_ADD - K_eff += (1-alpha_f)*gamma*dt*C_T \n"
        );

        // K_eff += (1-alpha_f)*beta*dt^2 * K_T
        // MATRIX CONTENTS -> K_eff = (1-alpha_m)*M + (1-alpha_f)*gamma*dt*C_T + (1-alpha_f)*beta*dt^2*K_T
        CALL_AND_CHECK_STATUS(
            mkl_sparse_d_add(
                                SPARSE_OPERATION_NON_TRANSPOSE,
                                K_T_mkl,
                                ( 1.0 - this->alpha_f ) * this->beta * pow2s( this->time_step ),
                                K_eff_mkl,
                                &K_eff_mkl
                            ),
            "Error after MKL_SPARSE_D_ADD - K_eff += (1-alpha_f)*beta*dt^2*K_T \n"
        );

        mkl_free( keff_row_cum_dum );

        // Extract K_eff as a CSRMatrix and solve the linear system:
        //   K_eff * delta_a = -R
        // where delta_a is the acceleration increment.
        // K_eff is assembled w.r.t. delta_a (see formulation above).
        CSRMatrix*  K_eff_csr   = convert_mkl_to_csrmatrix( &K_eff_mkl );
        SparseSolver linsolver  = SparseSolver( K_eff_csr, false );

        clear_vector( this->rows_np, delta_u );

        // Negate residual to form rhs: K_eff * delta_a = -R
        cusfloat* neg_residual = generate_empty_vector<cusfloat>( this->rows_np );
        for ( int i=0; i<this->rows_np; i++ )
        {
            neg_residual[i] = -residual[i];
        }
        this->_apply_restrictions( neg_residual );

        // delta_u here stores delta_a (the acceleration increment)
        linsolver.solve( neg_residual, delta_u );

        //------------------------------------------------------------
        // 4. Update kinematics from the acceleration increment delta_a
        //    (Newmark correctors):
        //
        //    a_{n+1} += delta_a
        //    v_{n+1} += gamma*dt        * delta_a
        //    u_{n+1} += beta*dt^2       * delta_a
        //
        //    Derivation: Newmark relations give
        //      u = u_pred + beta*dt^2 * a  →  delta_u = beta*dt^2 * delta_a
        //      v = v_pred + gamma*dt  * a  →  delta_v = gamma*dt  * delta_a
        //------------------------------------------------------------
        cusfloat coeff_u = this->beta  * pow2s( this->time_step );  // beta*dt^2
        cusfloat coeff_v = this->gamma * this->time_step;           // gamma*dt

        sv_add( this->rows_np, 1.0, this->y_acc, 1.0,    delta_u, this->y_acc );  // a += delta_a
        sv_add( this->rows_np, 1.0, this->y_vel, coeff_v, delta_u, this->y_vel ); // v += gamma*dt*delta_a
        sv_add( this->rows_np, 1.0, this->y_pos, coeff_u, delta_u, this->y_pos ); // u += beta*dt^2*delta_a

        // Apply restrictions to updated kinematics
        this->_apply_restrictions( this->y_pos );
        this->_apply_restrictions( this->y_vel );
        this->_apply_restrictions( this->y_acc );

        // Cleanup tangent matrices
        delete K_T;
        delete C_T;
        delete K_eff_csr;
        mkl_free( neg_residual );
        mkl_sparse_destroy( K_T_mkl );
        mkl_sparse_destroy( C_T_mkl );
        mkl_sparse_destroy( K_eff_mkl );

    }   // end Newton-Raphson loop

    // If not yet recorded (max_iter reached without convergence)
    if ( !this->nr_converged )
    {
        this->nr_iter       = this->max_iter;
        this->nr_res_abs    = res_norm;
        this->nr_res_rel    = res_norm / res_norm_0;

        std::cerr << "[WARN] GeneralizedAlphaNL::step - Newton-Raphson did not converge at"
                  << " t = " << this->time + this->time_step
                  << " | ||R|| = " << this->nr_res_abs
                  << " | ||R||/||R0|| = " << this->nr_res_rel
                  << " | iter = " << this->max_iter << "\n";
    }

    // Apply final restrictions
    this->_apply_restrictions( this->y_vel );
    this->_apply_restrictions( this->y_acc );

    //////////////////////////////////////////////////////////////////
    /////////////// UPDATE INTEGRATION VARIABLES /////////////////////
    //////////////////////////////////////////////////////////////////
    this->count++;
    this->time_old  = this->time;
    this->time      += this->time_step;

    // Cleanup
    mkl_free( u_pred   );
    mkl_free( v_pred   );
    mkl_free( a_pred   );
    mkl_free( residual );
    mkl_free( delta_u  );
}
