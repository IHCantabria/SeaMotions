
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
#include <iostream>

// Include local modules
#include "newmark.hpp"


inline void get_components_matrix( sparse_matrix_t lhs_mkl )
{
    struct matrix_descr descrA;
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;

    cusfloat e0[3] = {1.0, 0.0, 0.0};
    cusfloat e1[3] = {0.0, 1.0, 0.0};
    cusfloat e2[3] = {0.0, 0.0, 1.0};

    cusfloat r0[3];
    cusfloat r1[3];
    cusfloat r2[3];

    mkl_sparse_d_mv(
                        SPARSE_OPERATION_NON_TRANSPOSE,
                        1.0,
                        lhs_mkl,
                        descrA,
                        e0,
                        0.0,
                        r0
                    );
    mkl_sparse_d_mv(
                        SPARSE_OPERATION_NON_TRANSPOSE,
                        1.0,
                        lhs_mkl,
                        descrA,
                        e1,
                        0.0,
                        r1
                    );
    mkl_sparse_d_mv(
                        SPARSE_OPERATION_NON_TRANSPOSE,
                        1.0,
                        lhs_mkl,
                        descrA,
                        e2,
                        0.0,
                        r2
                    );
    
    print_vector(3, r0, 1, 6);
    print_vector(3, r1, 1, 6);
    print_vector(3, r2, 1, 6);
    
}


template<typename T>
void    NewmarkBeta<T>::_apply_restrictions(
                                                cusfloat*   _vec
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


template<typename T>
void    NewmarkBeta<T>::_build( 
                                        T               fext_in,
                                        CSRMatrix*      mass_mat_in,
                                        CSRMatrix*      stiff_mat_in,
                                        CSRMatrix*      damp_mat_in,
                                        cusfloat        time_step_in,
                                        cusfloat        t0_in,
                                        cusfloat*       y0_pos_in,
                                        cusfloat*       y0_vel_in,
                                        cusfloat*       y0_acc_in,
                                        int*            restrictions
                                )
{
    // Save input variables
    this->fext_fcn      =   fext_in;
    this->rows_np       =   mass_mat_in->rows_np;
    this->time_step     =   time_step_in;
    this->time_init     =   t0_in;
    this->y0_pos        =   y0_pos_in;
    this->y0_vel        =   y0_vel_in;
    this->y0_acc        =   y0_acc_in;
    this->beta          =   0.25;
    this->gamma         =   0.5;

    // Convert sparse input arrays from CSRMatrix to sparse_matrix_t type
    // and stores it into memory
    this->damp_mat = new sparse_matrix_t;
    CALL_AND_CHECK_STATUS(
        mkl_sparse_d_create_csr ( 
                                    this->damp_mat,
                                    SPARSE_INDEX_BASE_ZERO,
                                    damp_mat_in->rows_np,
                                    damp_mat_in->rows_np,
                                    damp_mat_in->row_index_cum,
                                    damp_mat_in->row_index_cum+1,
                                    damp_mat_in->col_index,
                                    damp_mat_in->values
                                ),
        "Error after MKL_SPARSE_D_CREATE_CSR - damp_mat \n"
    );
    
    
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
        "Error after MKL_SPARSE_D_CREATE_CSR - mass_mat \n"
    );

    this->stiff_mat = new sparse_matrix_t;
    CALL_AND_CHECK_STATUS(
        mkl_sparse_d_create_csr ( 
                                    this->stiff_mat,
                                    SPARSE_INDEX_BASE_ZERO,
                                    stiff_mat_in->rows_np,
                                    stiff_mat_in->rows_np,
                                    stiff_mat_in->row_index_cum,
                                    stiff_mat_in->row_index_cum+1,
                                    stiff_mat_in->col_index,
                                    stiff_mat_in->values
                                ),
        "Error after MKL_SPARSE_D_CREATE_CSR - stiff_mat \n"
    );

    // Storage restrictions pointer if not done previously
    if ( !this->_is_restrictions )
    {
        this->_restrictions = restrictions;
    }

    // Check for restrictions in input kinematics+
    this->_check_init_kinematics_retrictions( );

    // Initialize solver
    this->initialize( );
}


template<typename T>
void    NewmarkBeta<T>::_check_init_kinematics_retrictions(
                                                                void
                                                            )
{
    // Build warn message if any
    for ( int i=0; i<this->rows_np; i++ )
    {
        if ( this->_restrictions[i] > 0 )
        {
            // Check velocity restrictions inconsistency
            if ( check_zero_eps( this->y0_vel[i], EPS_PRECISION ) )
            {
                _print_inconsistency_msg( i, this->y0_vel[i], "Velocity" );
            }

            // Check acceleration restrictions inconsistency
            if ( check_zero_eps( this->y0_acc[i], EPS_PRECISION ) )
            {
                _print_inconsistency_msg( i, this->y0_acc[i], "Acceleration" );
            }
        }
    }
}


template<typename T>
void NewmarkBeta<T>::get_values_at(
                                        cusfloat  time_at,
                                        cusfloat* y_pos_at,
                                        cusfloat* y_vel_at,
                                        cusfloat* y_acc_at
                                    )
{
    // Get interpolation coefficient
    cusfloat alpha = ( time_at - this->time_old ) / ( this->time - this->time_old );

    // Perform interpolation
    for ( int i=0; i<this->rows_np; i++ )
    {
        y_pos_at[i] = ( 1 - alpha ) * this->y_pos_old[i] + alpha * this->y_pos[i];
        y_vel_at[i] = ( 1 - alpha ) * this->y_vel_old[i] + alpha * this->y_vel[i];
        y_acc_at[i] = ( 1 - alpha ) * this->y_acc_old[i] + alpha * this->y_acc[i];
    }
}


template<typename T>
void NewmarkBeta<T>::initialize( void )
{
    // Generate a matrix descriptor to use it in the 
    // mkl operation routines
    struct matrix_descr descrA;
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;

    //////////////////////////////////////////////////////////////////
    /////////////////////// CREATE LHS MATRIX/////////////////////////
    //////////////////////////////////////////////////////////////////
    // Create mkl sparse matrix handle to calculate the 
    // left hand side of the Newmark-beta method
    sparse_matrix_t lhs_mkl;

    // Add mass matrix to the left hand side
    // OPERATION -> LHS = M
    // MATRIX CONTENTS -> LHS = M
    CALL_AND_CHECK_STATUS(
        mkl_sparse_copy(
                            *(this->mass_mat),
                            descrA,
                            &lhs_mkl
                        ),
        "Error after MKL_SPARSE_COPY - lhs_mkl \n"
    );

    // Add damping matrix to the left hand side
    // OPERATION -> LHS += gamma*dt*C
    // MATRIX CONTENTS -> LHS = M + gamma*dt*C
    CALL_AND_CHECK_STATUS(
        mkl_sparse_d_add(
                            SPARSE_OPERATION_NON_TRANSPOSE,
                            *(this->damp_mat),
                            this->gamma*this->time_step,
                            lhs_mkl,
                            &lhs_mkl
                        ),
        "Error after MKL_SPARSE_D_ADD - lhs_mkl \n"
    );

    // Add stiffness matrix to the left hand side
    // OPERATION -> LHS += beta*dt^2*K
    // MATRIX CONTENTS -> LHS = M + gamma*dt*C + beta*dt^2*K
    CALL_AND_CHECK_STATUS(
        mkl_sparse_d_add(
                            SPARSE_OPERATION_NON_TRANSPOSE,
                            *(this->stiff_mat),
                            this->beta*pow2s(this->time_step),
                            lhs_mkl,
                            &lhs_mkl
                        ),
        "Error after MKL_SPARSE_D_ADD - lhs_mkl \n"
    );

    //////////////////////////////////////////////////////////////////
    //////////////// CREATE STIFF_DAMP_VEL MATRIX ////////////////////
    //////////////////////////////////////////////////////////////////
    // Create mkl sparse matrix handle to calculate the 
    // stiffness-damping-velocity matrix for the Newmark-beta method
    this->stiff_damp_vel            =   new sparse_matrix_t;
    int*    sdv_row_cum_dum         =   generate_empty_vector<int>(this->rows_np+1);
    int     sdv_col_index[1]        =   {0};
    cusfloat  sdv_values[1]         =   {0.0};
    sdv_row_cum_dum[this->rows_np]  =   1;

    CALL_AND_CHECK_STATUS(
        mkl_sparse_d_create_csr(
                                    this->stiff_damp_vel, 
                                    SPARSE_INDEX_BASE_ZERO, 
                                    this->rows_np, 
                                    this->rows_np, 
                                    sdv_row_cum_dum, 
                                    sdv_row_cum_dum+1, 
                                    sdv_col_index, 
                                    sdv_values
                                ),
        "Error after MKL_SPARSE_D_CREATE_CSR - stiff_damp_vel \n"
    );

    // Add damp matrix to the stiff_damp_vel matrix
    // OPERATION -> SDV = C
    // MATRIX CONTENTS -> SDV = C
    CALL_AND_CHECK_STATUS(
        mkl_sparse_d_add(
                            SPARSE_OPERATION_NON_TRANSPOSE,
                            *(this->damp_mat),
                            1.0,
                            *(this->stiff_damp_vel),
                            this->stiff_damp_vel
                        ),
        "Error after MKL_SPARSE_D_ADD - stiff_damp_vel \n"
    );

    // Add stiffness matrix to the stiff_damp_vel matrix
    // OPERATION -> SDV += dt*K
    // MATRIX CONTENTS -> SDV = C + dt*K
    CALL_AND_CHECK_STATUS(
        mkl_sparse_d_add(
                            SPARSE_OPERATION_NON_TRANSPOSE,
                            *(this->stiff_mat),
                            this->time_step,
                            *(this->stiff_damp_vel),
                            this->stiff_damp_vel
                        ),
        "Error after MKL_SPARSE_D_ADD - stiff_damp_vel \n"
    );

    //////////////////////////////////////////////////////////////////
    //////////////// CREATE STIFF_DAMP_ACC MATRIX ////////////////////
    //////////////////////////////////////////////////////////////////
    // Create mkl sparse matrix handle to calculate the 
    // stiffness-damping-acceleration matrix for the Newmark-beta method
    this->stiff_damp_acc            =   new sparse_matrix_t;
    int*    sda_row_cum_dum         =   generate_empty_vector<int>(this->rows_np+1);
    int     sda_col_index[1]        =   {0};
    cusfloat  sda_values[1]         =   {0.0};
    sda_row_cum_dum[this->rows_np]  =   1;


    CALL_AND_CHECK_STATUS(
        mkl_sparse_d_create_csr(
                                    this->stiff_damp_acc, 
                                    SPARSE_INDEX_BASE_ZERO, 
                                    this->rows_np, 
                                    this->rows_np, 
                                    sda_row_cum_dum, 
                                    sda_row_cum_dum+1, 
                                    sda_col_index, 
                                    sda_values
                                ),
        "Error after MKL_SPARSE_D_CREATE_CSR - stiff_damp_acc \n"
    );

    // Add damp matrix to the stiff_damp_acc matrix
    // OPERATION -> SDA = dt*(1-gamma)*C
    // MATRIX CONTENTS -> SDA = dt*(1-gamma)*C
    CALL_AND_CHECK_STATUS(
        mkl_sparse_d_add(
                            SPARSE_OPERATION_NON_TRANSPOSE,
                            *(this->damp_mat),
                            this->time_step*(1-this->gamma),
                            *(this->stiff_damp_acc),
                            this->stiff_damp_acc
                        ),
        "Error after MKL_SPARSE_D_ADD - stiff_damp_acc \n"
    );

    // Add stiff matrix to the stiff_damp_acc matrix
    // OPERATION -> SDA += 0.5*dt**2.0*(1-2*beta)*K
    // MATRIX CONTENTS -> SDA = dt*(1-gamma)*C + 0.5*dt**2.0*(1-2*beta)*K
    CALL_AND_CHECK_STATUS(
        mkl_sparse_d_add(
                            SPARSE_OPERATION_NON_TRANSPOSE,
                            *(this->stiff_mat),
                            0.5*pow2s(this->time_step)*(1-2*this->beta),
                            *(this->stiff_damp_acc),
                            this->stiff_damp_acc
                        ),
        "Error after MKL_SPARSE_D_ADD - stiff_damp_acc \n"
    );

    //////////////////////////////////////////////////////////////////
    //////////// EXTRACT sparse data to CSRMatrix ////////////////////
    //////////////////////////////////////////////////////////////////
    this->lhs_csr               = convert_mkl_to_csrmatrix(&lhs_mkl);

    //////////////////////////////////////////////////////////////////
    ///////////////// INITIALIZE SPARSE SOLVER ///////////////////////
    //////////////////////////////////////////////////////////////////
    this->linsolve              = new SparseSolver(lhs_csr, false);

    //////////////////////////////////////////////////////////////////
    ///////////// CORRECT INITIAL ACCELERATION  //////////////////////
    //////////////////////////////////////////////////////////////////
    /* In this section it will be corrected the acceleration in order 
    to force the initial conditions to satisfy the the Newton 
    second law equation: M·a + C·v + K·d = F */

    // Start an instance of the SparseSolver with the mass matrix as a 
    // system matrix
    CSRMatrix*      mass_csr    = convert_mkl_to_csrmatrix(this->mass_mat);
    SparseSolver    acc_crr_sps = SparseSolver(mass_csr, false);

    // Create rhs vector to storage the result of the sum of the
    // damping, stiffness and external forces
    cusfloat* acc_crr_rhs       = generate_empty_vector<cusfloat>(this->rows_np);

    // Add external forces contribution
    // OPERATION -> ACC_CRR_RHS = F_ext
    // VECTOR CONTENTS -> ACC_CRR_RSH = F_ext
    (*this->fext_fcn)(
                            this->time_old,
                            this->time_step,
                            this->y0_pos,
                            this->y0_vel,
                            this->y0_acc,
                            acc_crr_rhs
                        );

    // Add stiffnes forces for the rhs vector
    // OPERATION -> ACC_CRR_RHS += -K*d0
    // VECTOR CONTENTS -> ACC_CRR_RSH = F_ext - K*d0
    CALL_AND_CHECK_STATUS(
        mkl_sparse_d_mv(
                            SPARSE_OPERATION_NON_TRANSPOSE,
                            -1.0,
                            *(this->stiff_mat),
                            descrA,
                            this->y0_pos,
                            1.0,
                            acc_crr_rhs
                        ),
        "Error after MKL_SPARSE_D_MV - ACC_CRR_RHS = -K*d0 \n"
    );

    // Add stiffnes forces for the rhs vector
    // OPERATION -> ACC_CRR_RHS += -C*v0
    // VECTOR CONTENTS -> ACC_CRR_RSH = F_ext - K*d0 - C*v0
    CALL_AND_CHECK_STATUS(
        mkl_sparse_d_mv(
                            SPARSE_OPERATION_NON_TRANSPOSE,
                            -1.0,
                            *(this->damp_mat),
                            descrA,
                            this->y0_vel,
                            1.0,
                            acc_crr_rhs
                        ),
        "Error after MKL_SPARSE_D_MV - ACC_CRR_RHS += -C*v0 \n"
    );

    // Solve the system of equations: M·a0 = F_ext - K·d0 - C·v0
    // in order to get the initial acceleration
    clear_vector( this->rows_np, this->y0_acc );
    acc_crr_sps.solve( acc_crr_rhs, this->y0_acc );

    //////////////////////////////////////////////////////////////////
    /////////////// RESET INTEGRATION VARIABLES //////////////////////
    //////////////////////////////////////////////////////////////////
    // Reset integration variables
    this->time      = this->time_init;
    this->time_old  = this->time_init;
    this->count     = 0;
    this->y_pos     = generate_empty_vector<cusfloat>( this->rows_np );
    this->y_vel     = generate_empty_vector<cusfloat>( this->rows_np );
    this->y_acc     = generate_empty_vector<cusfloat>( this->rows_np );
    this->y_pos_old = generate_empty_vector<cusfloat>( this->rows_np );
    this->y_vel_old = generate_empty_vector<cusfloat>( this->rows_np );
    this->y_acc_old = generate_empty_vector<cusfloat>( this->rows_np );

    // Copy initial values in the current step array
    copy_vector( this->rows_np, this->y0_pos, this->y_pos );
    copy_vector( this->rows_np, this->y0_vel, this->y_vel );
    copy_vector( this->rows_np, this->y0_acc, this->y_acc );

    // Initialize right hand side vector on heap memory in 
    // order to avoid of regenerate this memory during the 
    // time steps
    this->rhs       = generate_empty_vector<cusfloat>( this->rows_np );

    // Allocate heap memory for u_gamma and u_beta auxiliar 
    // vectors to storage information about kinematics
    this->_u_beta   = generate_empty_vector<cusfloat>( this->rows_np );
    this->_u_gamma  = generate_empty_vector<cusfloat>( this->rows_np );

    // Delete heap memory of variables to be destroyed along 
    // this function
    mkl_free( acc_crr_rhs       );
    mkl_free( sdv_row_cum_dum   );
    mkl_free( sda_row_cum_dum   );

}


template<typename T>
void NewmarkBeta<T>::interpolate_solution(
                                                cusfloat time, 
                                                cusfloat* y_pos_itp,
                                                cusfloat* y_vel_itp,
                                                cusfloat* y_acc_itp
                                        )
{
    // Check for time bounds
    cusfloat TIME_PRECISION = 1e-6;
    if ( time < this->time_old-TIME_PRECISION )
    {
        std::cerr << std::setprecision(16) << std::endl;
        std::cerr << "NewmarkBeta::interpolate_solution - ";
        std::cerr << "Interpolation time is below time_old value: ";
        std::cerr << "time_old [s] = " << this->time_old << " - time_itp = ";
        std::cerr << time << std::endl;
        exit(300);
    }
    if ( time > this->time+TIME_PRECISION )
    {
        std::cerr << std::setprecision(16) << std::endl;
        std::cerr << "NewmarkBeta::interpolate_solution - ";
        std::cerr << "Interpolation time is above time value: ";
        std::cerr << "time_old [s] = " << this->time_old << " - time_itp = ";
        std::cerr << time << std::endl;
        exit(300);
    }

    // Calculate alpha value
    cusfloat alpha = ( time - this->time_old ) / this->time_step;

    // Clear target vectors in order to avoid spurious data
    // to sum
    clear_vector( this->rows_np, y_pos_itp );
    clear_vector( this->rows_np, y_vel_itp );
    clear_vector( this->rows_np, y_acc_itp );

    // Perform interpolation
    sv_add(
                this->rows_np,
                (1-alpha),
                this->y_pos_old,
                alpha,
                this->y_pos,
                y_pos_itp
            );
    
    sv_add(
                this->rows_np,
                (1-alpha),
                this->y_vel_old,
                alpha,
                this->y_vel,
                y_vel_itp
            );
    
    sv_add(
                this->rows_np,
                (1-alpha),
                this->y_acc_old,
                alpha,
                this->y_acc,
                y_acc_itp
            );
    
}


template<typename T>
NewmarkBeta<T>::NewmarkBeta(
                                    T               fext_in,
                                    CSRMatrix*      mass_mat_in,
                                    CSRMatrix*      stiff_mat_in,
                                    CSRMatrix*      damp_mat_in,
                                    cusfloat        time_step_in,
                                    cusfloat        t0_in,
                                    cusfloat*       y0_pos_in,
                                    cusfloat*       y0_vel_in,
                                    cusfloat*       y0_acc_in
                        )
{
    // Buid restrictions vector assuming no restrictions are imposed
    this->_restrictions     = generate_empty_vector<int>( mass_mat_in.rows_np );
    this->_is_restrictions  = true;

    // Call class builder
    this->_build( 
                    fext_in,
                    mass_mat_in,
                    stiff_mat_in,
                    damp_mat_in,
                    time_step_in,
                    t0_in,
                    y0_pos_in,
                    y0_vel_in,
                    y0_acc_in,
                    this->_restrictions
                );
    
}


template<typename T>
NewmarkBeta<T>::NewmarkBeta(
                                    T               fext_in,
                                    CSRMatrix*      mass_mat_in,
                                    CSRMatrix*      stiff_mat_in,
                                    CSRMatrix*      damp_mat_in,
                                    cusfloat        time_step_in,
                                    cusfloat        t0_in,
                                    cusfloat*       y0_pos_in,
                                    cusfloat*       y0_vel_in,
                                    cusfloat*       y0_acc_in,
                                    int*            restrictions_in
                        )
{
    // Call class builder
    this->_build( 
                    fext_in,
                    mass_mat_in,
                    stiff_mat_in,
                    damp_mat_in,
                    time_step_in,
                    t0_in,
                    y0_pos_in,
                    y0_vel_in,
                    y0_acc_in,
                    restrictions_in
                );
    
}


template<typename T>
NewmarkBeta<T>::~NewmarkBeta( void )
{
    // Delete restrictions matrix if any
    if ( this->_is_restrictions )
    {
        mkl_free( this->_restrictions );
    }

    // Delete sparse solve
    delete this->linsolve;

    // Delete heap memory associated to system sparse matrixes
    mkl_sparse_destroy( *( this->damp_mat ) );
    delete this->lhs_csr;
    mkl_sparse_destroy( *( this->mass_mat       ) );
    mkl_sparse_destroy( *( this->stiff_damp_acc ) );
    mkl_sparse_destroy( *( this->stiff_damp_vel ) );
    mkl_sparse_destroy( *( this->stiff_mat      ) );

    // Delete heap memory associated with the right hand side
    // and the auxiliar vectors
    mkl_free( this->rhs             );
    mkl_free( this->y_acc           );
    mkl_free( this->y_pos           );
    mkl_free( this->y_vel           );
    mkl_free( this->y_acc_old       );
    mkl_free( this->y_pos_old       );
    mkl_free( this->y_vel_old       );
    mkl_free( this->_u_beta         );
    mkl_free( this->_u_gamma        );
}


template<typename T>
void NewmarkBeta<T>::step( void )
{
    // Generate matrix descriptor for the mkl operation 
    // routines
    struct matrix_descr descrA;
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;

    //////////////////////////////////////////////////////////////////
    ///// COPY PREVIOUS STEP DATA INTO THE OLD VECTOR HANDLE /////////
    //////////////////////////////////////////////////////////////////
    copy_vector( this->rows_np, this->y_pos, this->y_pos_old );
    copy_vector( this->rows_np, this->y_vel, this->y_vel_old );
    copy_vector( this->rows_np, this->y_acc, this->y_acc_old );

    clear_vector( this->rows_np, this->y_pos );
    clear_vector( this->rows_np, this->y_vel );
    clear_vector( this->rows_np, this->y_acc );

    //////////////////////////////////////////////////////////////////
    //////////// CALCULATE RIGHT HAND SIDE VECTOR ////////////////////
    //////////////////////////////////////////////////////////////////
    // Clear rhs vector in order to use it along the
    // sumation process withouh taking into account
    // spurious data from previous steps or the initalization
    clear_vector( this->rows_np, this->rhs );

    // Add external forces contribution
    ( *this->fext_fcn )(
                            this->time,
                            this->time_step,
                            this->y_pos_old,
                            this->y_vel_old,
                            this->y_acc_old,
                            this->rhs
                        );

    // Apply restrictions to external forces
    this->_apply_restrictions( this->rhs );

    // Add first summand with the stiffness contibution
    // OPEARTION -> RHS = -K*u_pos
    // MATRIX CONTENTS -> RHS = -K*u_pos
    CALL_AND_CHECK_STATUS(
        mkl_sparse_d_mv(
                            SPARSE_OPERATION_NON_TRANSPOSE,
                            -1.0,
                            *(this->stiff_mat),
                            descrA,
                            this->y_pos_old,
                            1.0,
                            this->rhs
                        ),
        "Error after MKL_SPARSE_D_MV - RHS = -K*u_pos \n"
    );

    // Add first summand with the stiffness contibution
    // OPEARTION -> RHS += -SDV*u_vel
    // MATRIX CONTENTS -> RHS = -K*u_pos - SDV*u_vel
    CALL_AND_CHECK_STATUS(
        mkl_sparse_d_mv(
                            SPARSE_OPERATION_NON_TRANSPOSE,
                            -1.0,
                            *(this->stiff_damp_vel),
                            descrA,
                            this->y_vel_old,
                            1.0,
                            this->rhs
                        ),
        "Error after MKL_SPARSE_D_MV - RHS += -SDV*u_vel \n"
    );

    // Add first summand with the stiffness contibution
    // OPEARTION -> RHS += -SDA*u_acc
    // MATRIX CONTENTS -> RHS = -K*u_pos - SDV*u_vel - SDA*u_acc
    CALL_AND_CHECK_STATUS(
        mkl_sparse_d_mv(
                            SPARSE_OPERATION_NON_TRANSPOSE,
                            -1.0,
                            *(this->stiff_damp_acc),
                            descrA,
                            this->y_acc_old,
                            1.0,
                            this->rhs
                        ),
        "Error after MKL_SPARSE_D_MV - RHS += -SDA*u_acc \n"
    );

    //////////////////////////////////////////////////////////////////
    ///////////////// CALCULATE NEW KINEMATICS ///////////////////////
    //////////////////////////////////////////////////////////////////
    // Calculate new acceleration
    this->linsolve->solve( this->rhs, this->y_acc );

    // Calculate u_gamma parameter has the combination of the
    // the gamma parameter and the accelerations
    // OPERATION -> u_gamma = (1-gamma)*u_acc_old + gamma*u_acc
    sv_add(
                this->rows_np,
                (1-this->gamma),
                this->y_acc_old,
                this->gamma,
                this->y_acc,
                this->_u_gamma
            );

    // Calculate new velocity as a function of the old one and 
    // the time step
    // OPERATION -> u_vel = u_vel_old + dt*u_gamma
    sv_add(
                this->rows_np,
                1.0,
                this->y_vel_old,
                this->time_step,
                this->_u_gamma,
                this->y_vel
            );

    // Calculate u_beta parameter as a function of the 
    // old and the new acceleration values
    // OPERATION -> u_beta = (1-2*beta)*u_acc_old + 2*beta*u_acc
    sv_add(
                this->rows_np,
                (1-2*this->beta),
                this->y_acc_old,
                2*this->beta,
                this->y_acc,
                this->_u_beta
            );
    
    // Calculate new position as a function of the old position,
    // velocity and acceleration
    // OPERATION -> u_pos = u_pos_old + dt*u_vel_old
    // VECTOR CONTENTS -> u_pos = u_pos_old + dt*u_vel_old
    sv_add(
                this->rows_np,
                1.0,
                this->y_pos_old,
                this->time_step,
                this->y_vel_old,
                this->y_pos
            );

    // Calculate new position as a function of the old position,
    // velocity and acceleration
    // OPERATION -> u_pos += 0.5*dt**2.0*u_beta
    // VECTOR CONTENTS -> u_pos = u_pos_old + dt*u_vel_old + 0.5*dt**2.0*u_beta
    sv_add(
                this->rows_np,
                1.0,
                this->y_pos,
                0.5*pow2s(this->time_step),
                this->_u_beta,
                this->y_pos
            );

    // Apply restrictions to predicted kinematics
    this->_apply_restrictions( this->y_vel );
    this->_apply_restrictions( this->y_acc );

    //////////////////////////////////////////////////////////////////
    /////////////// UPDATE INTEGRATION VARIABLES /////////////////////
    //////////////////////////////////////////////////////////////////
    this->count++;
    this->time_old  = this->time;
    this->time      += this->time_step;

}