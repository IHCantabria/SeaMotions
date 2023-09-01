
#ifndef __nonlinear_solvers_hpp
#define __nonlinear_solvers_hpp

// Include general usage libraries
#include <functional>

// Include general usage scientific libraries
#include "mkl.h"

// Include local modules
#include "math_tools.hpp"


// Deinfe macro to simply the sparse operation status
#define CHECK_NLINSOL_STATUS( function, is_debug, error_message ) do {  \
        if(function > 0)                                                \
        {                                                               \
            if ( is_debug )                                             \
            {                                                           \
                std::cout << "FILE: " << __FILE__ << " - ";             \
                std::cout << "LINE: " << __LINE__ << " - ";             \
            }                                                           \
        std::cout << error_message << "\n";                             \
        }                                                               \
} while( 0 )

#define SHOW_ERROR_MESAGE( error_message ) do {             \
          std::cerr << "FILE: " << __FILE__ << " - ";       \
          std::cerr << "LINE: " << __LINE__ << " - ";       \
          std::cerr << error_message << "\n";               \
} while( 0 )


// Define the type for the external function pointer in the case of 
// the Trust Region algorithm
typedef void (*fdef_type_tr)(int*, int*, double*, double*);


////////////////////////////////////////////////////////////////
/******************** Declare functions ***********************/
////////////////////////////////////////////////////////////////
template<typename T>
inline void     calculate_jacobi_mat(
                                    T       fdef,       // in:     external objective function
                                    int     nfi,       // in:     number of function variables
                                    int     nfo,       // in:     dimension of function value
                                    int     ordering,   // in:      Row major or column major
                                    double* fjac,       // out:    jacobi matrix
                                    double* x,          // in:     solution vector
                                    double  jac_eps    // in:     jacobi calculation precision
                                    )
{
    // Generate auxiliar input and output vectors
    double* fm = generate_empty_vector<double>( nfo );
    double* fp = generate_empty_vector<double>( nfo );
    double* xm = generate_empty_vector<double>( nfi );
    double* xp = generate_empty_vector<double>( nfi );

    // Loop over input and output number of variables to 
    // complete the jacobian matrix calculation
    for ( int i=0; i<nfo; i++ )
    {
        for ( int j=0; j<nfi; j++ )
        {
            // Calculate input vectors
            for ( int k=0; k<nfi; k++ )
            {
                xm[k] = x[k];
                xp[k] = x[k];
            }
            xm[j] -= jac_eps;
            xp[j] += jac_eps;

            // Calculate function value
            fdef( &nfo, &nfi, xm, fm );
            fdef( &nfo, &nfi, xp, fp );
            
            // Calculate jacobian
            if ( ordering == 0) // Row major
            {
                fjac[i*nfi+j] = (fp[i]-fm[i])/2.0/jac_eps;
            }
            else
            {
                fjac[j*nfo+i] = (fp[i]-fm[i])/2.0/jac_eps;
            }
        }
    }

    // Delete auxiliar vectors
    mkl_free( fm );
    mkl_free( fp );
    mkl_free( xm );
    mkl_free( xp );
}


template<typename T>
inline  void    bisection(
                            T       f_def, 
                            double  a, 
                            double  b, 
                            double  fabs_tol,
                            double  xrel_tol, 
                            int     max_iter, 
                            bool    verbose,
                            double  &sol, 
                            int     &info
                            )
{
    /**
     * @brief Bisection method to solve non-linear monoparametric equations.
     * 
     * The bisection method solves non-linear equations dependent on one parameter. It is an 
     * iterative solver based on the interval sub-division. More information can be found 
     * at: Numerical Recipes The Art of Scientific Computing. W.Press S.Teukolsky
     * 
     * \param f_def_raw Target function
     * \param a Lower bound of the interval in which is assumed to have the zero.
     * \param b Upper bound of the interval in which is assumed to have the zero.
     * \param fabs_tol Absolute value tolerance of target function value to stop iterations: f_def(x) < fabs_tol
     * \param xrel_tol Relative tolerance of the independent variable to stop the iterations: dx < x*xrel_tol
     * \param max_iter Maximum iterations allowed to find the zero of the target function.
     * \param verbose Bool flag to print out to screen the iterative path to find the solution.
     * \param sol Scalar parameter in which store the solution of the iterative process.
     * \param info This parameter describes the finish status of the iterative process.
     *              - info == 0 -> Solution found.
     *              - info == 1 -> Solution not found or convergence problems.
     */

    // Start the iterative method
    double fa = f_def(a);
    double c, cb=(a+b)/2.0+2*xrel_tol+1.0;
    double fc;

    // Start iterative loop
    double abs_err = 0.0;
    int count_iter = 0;
    info = 0;
    double rel_err = 0.0;
    while (true)
    {
        // Calculate new value
        c = (a+b)/2.0;
        fc = f_def(c);
        
        // Calculate errors
        abs_err = abs(fc);
        rel_err = abs(c-cb);

        // Print iterative process
        if (verbose)
        {
            std::cout << "Iter: " << count_iter << " - x: " << c << " - f(x): " << fc << std::endl;
            std::cout << "  -> a: " << a << " - f(a): " << f_def(a) << std::endl;
            std::cout << "  -> b: " << b << " - f(b): " << f_def(b) << std::endl;
            std::cout << "      - Abs.Error: " << abs_err << " - Rel.Error: " << rel_err << std::endl;
        }

        // Check for convergence
        if ((abs(fc)<=fabs_tol) || (abs(c-cb)<abs(c*xrel_tol+1e-14)))
        {
            break;
        }

        // Check for maximum iterations limit
        if (count_iter > max_iter)
        {
            std::cerr << "WARNING: Bisection method could not find the solution ";
            std::cerr << "with the accurancy requested. Residual Value: " << fc << std::endl;
            info = 1;
            break;
        }
        count_iter++;

        // Update interval values for the next iteration
        cb = c;
        if (fa*fc>0)
        {
            a = c;
            fa = fc;
        }
        else
        {
            b = c;
        }
    }

    // Store solution in the output variable
    sol = c;
}


template<typename T>
inline  void    newton_raphson(
                                T       f_def,
                                T       f_der_def,
                                double  x0, 
                                double  fabs_tol,
                                double  xrel_tol,
                                int     max_iter,
                                bool    verbose,
                                double  &sol, 
                                int     &info
                                )
{
    /**
     * @brief Newthon-Rapshon method to solve non-linear equations of the type x=T(x).
     * 
     * This function implements the Newthon-Rapshon method to solve non-linear equations. The
     * solver is prepared to solve monoparametric equations (1D). The solver is based on the 
     * clasical algorithm that can be found at: 
     *  - Numerical Recipes The Art of Scientific Computing. W.Press S.Teukolsky
     * 
     * \param f_def Target function definition
     * \param f_der_der Target function derivative definition
     * \param x0 First pont to start iterations
     * \param fabs_tol Absolute value tolerance of target function value to stop iterations: f_def(x) < fabs_tol
     * \param xrel_tol Relative tolerance of the independent variable to stop the iterations: dx < x*xrel_tol
     * \param max_iter Maximum iterations allowed
     * \param verbose Flag to activate the iterative process print out to screen.
     * \param sol Solution to the equation provided.
     * \param info This parameter describes the finish status of the iterative process.
     *              - info == 0 -> Solution found.
     *              - info == 1 -> Solution not found or convergence problems.
     * 
     */
    // Set info variable to sucess. It will be set to failed
    // if the convergence is poor or it fails.
    info = 0;

    // Perform iterative loop to look for the solution
    double x = x0, dx = 0.0;
    int count_iter = 0;
    while (true)
    {
        // Calculate new increment
        dx = f_def(x)/f_der_def(x);

        // Calculate new solution point
        x -= dx;

        // Print-out to screen the iterative process status
        if (verbose)
        {
            std::cout << "Iter: " << count_iter << " - fabs(x): " << f_def(x);
            std::cout << " - x: " << x << " - dx: " << dx << std::endl; 
        }

        // Check for convergence
        if ((std::abs(f_def(x))<fabs_tol || (std::abs(dx)<=std::abs(xrel_tol*x+1e-14))))
        {
            break;
        }

        // Update iterations and check for limits
        if (count_iter > max_iter)
        {
            std::cerr << "WARNING: Newton-Raphson method could not find the solution ";
            std::cerr << "with the accurancy requested. Residual Value: " << f_def(x) << std::endl;
            info = 1;
            break;
        }
        count_iter++;
    }

    // Storage result in output channel
    sol = x;
}


template<typename T>
inline  void newton_raphson_nd(
                                T       f_def,
                                double* x0,
                                int     nfi,
                                int     nfo,
                                double* fabs_tol,
                                double* frel_tol,
                                double* xrel_tol,
                                int     max_iter,
                                bool    verbose,
                                double* sol, 
                                int     &info
                                )
{
    // Check if there are the input arguments as the output ones
    if ( nfi != nfo )
    {
        std::cerr << "NetonRapshon: the number of function inputs and outputs must be the same." << std::endl;
        info = 5;
        return;
    }

    // Configure solver settings
    int     count       = 0;
    double  jac_eps     = 1e-6;

    // Allocate internal vectors to storage the function value and 
    // the jacobian matrix
    double* f_jac = generate_empty_vector<double>( nfi*nfo );
    double* f_vec = generate_empty_vector<double>( nfo );

    // Configure linear equations solver
    int*    ipiv        = generate_empty_vector<int>( nfi );
    int     info_ls     = 0;
    int     rows_np     = nfi;
    int     rhs_np      = 1;

    // Calculate first solution vector
    f_def( &nfo, &nfi, x0, f_vec);

    // Declare some auxiliar variables to work during the 
    // iterative loop
    bool    fabs_pass   = false;
    bool    frel_pass   = false;
    bool    xrel_pass   = false;
    int     index       = 0;
    double  err_i       = 0.0;
    double* f_last      = generate_empty_vector<double>( nfi );
    double* x1          = generate_empty_vector<double>( nfi );
    double* dx          = generate_empty_vector<double>( nfi );

    // Loop in order to find the solution
        // Calculate new jacobian
    while ( true )
    {
        // Clear solution vector
        clear_vector( nfo, x1 );

        // Calculate first jacobian matrix
        calculate_jacobi_mat( f_def, nfi, nfo, 1, f_jac, x0, jac_eps );

        // Calculate new step
        svs_mult( nfo, f_vec, -1.0, f_vec);
        dgesv( &rows_np, &rhs_np, f_jac, &rows_np, ipiv, f_vec, &rows_np, &info_ls );
        sv_add( rows_np, x0, f_vec, x1 );

        // Copy last step increment into a new variable to reuse it in
        // the error check section
        copy_vector( nfi, f_vec, dx );

        // Evaluate right hand side
        clear_vector( nfo, f_vec );
        f_def( &nfo, &nfi, x1, f_vec);

        // Calculate errors
        index = cblas_idamax( rows_np, f_vec, 1 );
        err_i = f_vec[index];
        if ( verbose )
        {
            std::cout << "Iteration: " << count << " - F AbsErr: " << err_i << std::endl;
        }

        // Check tolerances
        fabs_pass = true;
        frel_pass = true;
        xrel_pass = true;

        for ( int i=0; i<nfi; i++ )
        {
            // Check absolute value of the function
            if ( std::abs( f_vec[i] ) > fabs_tol[i] )
            {
                fabs_pass = false;
            }

            // Check relative value of the function
            if ( std::abs( ( f_vec[i] - f_last[i] )/f_last[i] ) > frel_tol[i] )
            {
                frel_pass = false;
            }

            // Check relative value of the input variables
            if ( std::abs( dx[i] ) > xrel_tol[i] )
            {
                xrel_pass = false;
            }
        }

        if ( fabs_pass || frel_pass || xrel_pass )
        {
            break;
        }

        // Copy new solution into the old one to perform the next step
        copy_vector( rows_np, x1, x0 );
        
        // Advance iteration and check for limits
        count += 1;
        if ( count > max_iter )
        {
            // std::cout << std::abs( err_i ) << " - " << std::abs( ( err_i - err_last )/err_i ) << std::endl;
            // std::cout << x1[0] << " - " << x1[1] << std::endl;
            // std::cout << dx[0] << " - " << dx[1] << std::endl;
            std::cerr << "Maximum iterations reached" << std::endl;
            exit(200);
        }

        copy_vector( nfi, f_vec, f_last );
    }

    // Copy data from x1 to sol vector
    copy_vector( nfi, x1, sol);

    // Deallocate heap memory matrixes
    mkl_free( dx );
    mkl_free( f_jac );
    mkl_free( f_last );
    mkl_free( f_vec );
    mkl_free( ipiv );
    mkl_free( x1 );
}


template<typename T>
inline  void trust_region(
                    T               fdef,
                    int             nfi,
                    int             nfo,
                    double*         x0,
                    double*         x,
                    int             &tr_status,
                    bool            verbose,
                    double          fabs_tol=1e-5,
                    double          frel_tol=1e-5,
                    double          xrel_tol=1e-5
                )
{
    /**
     * @brief Trust region algorithm to solve non linear equations.
     * 
     * \param fext      Function definition to find the roots
     * \param nfi       Number of function input variables.
     * \param nfo       Number of function ouput values.
     * \param x0        Initial guess for the solution.
     * \param x         Solution of the system of equations
     * \param tr_status State variable containing the end status of the TR method
     * \param verbose   Switch to activate the verbose mode to get more insight on the program status
     */
    ////////////////////////////////////////////////////////////////
    /************ Set solver configuration settings ***************/
    ////////////////////////////////////////////////////////////////
    const double    eps[6]      = {xrel_tol, fabs_tol, 1e-6, 1e-6, frel_tol, 1e-6}; // set precisions for stop-criteria
    const int       iter1       = 1000;                                             // iter1 - maximum number of iterations
    const int       iter2       = 100;                                              // iter2 - maximum number of iterations of calculation of trial-step
    double          jac_eps     = 1e-6;                                             // set precision of the Jacobian matrix calculation
    int             mem_error   = 0;                                                // memory allocated correctly
    double          rs          = 0.0;                                              // initial step bound
    

    ////////////////////////////////////////////////////////////////
    /** Initialize solver (allocate memory, set initial values) ***/
    ////////////////////////////////////////////////////////////////
    // Initialize local variables to use along the function in order
    // to manage the execution of the solver
    int             info[6];            // results of input parameter checking
    int             iter        = 0;    // number of iterations (out statistics)
    int             error       = 0;    // variable to manage the error code in case of failure
    double          r1          = 0.0;  // Initial residual (out statistics)
    double          r2          = 0.0;  // Final residual (out statistics)
    int             RCI_Request = 0;    // reverse communication interface variable
    int             res         = -1;   // variable to manage execution status
    int             st_cr       = 0;    // number of stop criterion (out statistics)
    int             successful  = 0;    // Management parameter for the rci cycle stage
    _TRNSP_HANDLE_t handle;             // TR solver handle

    // Copy initial data into the solution vector
    for ( int i=0; i<nfi; i++ )
    {
        x[i] = x0[i];
    }

    // Allocate memory for function vector values and the jacobian
    // matrix
    double* fjac = generate_empty_vector<double>( nfi*nfo );
    double* fvec = generate_empty_vector<double>( nfo );
    
    res = dtrnlsp_init (
                        &handle,    // in|out:  TR solver handle
                        &nfi,       // in:      number of function variables
                        &nfo,       // in:      dimension of function value
                        x,          // in:      solution vector. contains values x for f(x)
                        eps,        // in:      precisions for stop-criteria
                        &iter1,     // in:      maximum number of iterations
                        &iter2,     // in:      maximum number of iterations of calculation of trial-step
                        &rs         // in:      initial step bound
                        );
    if ( res != TR_SUCCESS )
    {
        printf ("| \n");
        SHOW_ERROR_MESAGE( "Error in dtrnlsp_init" );
        error = 1;
        goto end;
    }

    // Checks the correctness of handle and arrays containing Jacobian matrix,
    // objective function, lower and upper bounds, and stopping criteria
    res = dtrnlsp_check (
                        &handle,    // in|out:  TR solver handle
                        &nfi,       // in:      number of function variables
                        &nfo,       // in:      dimension of function value
                        fjac,       // in:      jacobian of the target function
                        fvec,       // in:      values of the target function
                        eps,        // in:      precisions for stop-criteria
                        info        // out:     results of input parameter checking
                        );
    if ( res != TR_SUCCESS )
    {
        SHOW_ERROR_MESAGE( "Error in dtrnlspbc_init" );
        error = 1;
        goto end;
    }
    else
    {
        if (info[0] != 0 || // The handle is not valid.
            info[1] != 0 || // The fjac array is not valid.
            info[2] != 0 || // The fvec array is not valid.
            info[3] != 0    // The eps array is not valid.
           )
        {
            if ( info[0] != 0 )
            {
                std::cerr << "The handle is not valid." << std::endl;
            }
            else if ( info[1] != 0 )
            {
                std::cerr << "The fjac array is not valid." << std::endl;
            }
            else if ( info[2] != 0 )
            {
                std::cerr << "The fvec array is not valid." << std::endl;
            }
            else if ( info[3] != 0)
            {
                std::cerr << "The eps array is not valid." << std::endl;
            }
            SHOW_ERROR_MESAGE( "Input parameters for dtrnlsp_solve are not valid" );
            error = 1;
            goto end;
        }
    }


    ////////////////////////////////////////////////////////////////
    /**************** Solve system of equations *******************/
    ////////////////////////////////////////////////////////////////
    /* set initial rci cycle variables */
    

    // Perform rci cycle
    while ( successful == 0 )
    {
        // Call TR solver
        res = dtrnlsp_solve(
                            &handle,        // in|out:  TR solver handle
                            fvec,           // in:      values of the target function
                            fjac,           // in:      jacobian of the target function
                            &RCI_Request    // in|out:  return number which denote next step for performing
                            );
        if ( res != TR_SUCCESS )
        {
            SHOW_ERROR_MESAGE( "Error in dtrnlsp_solve" );
            error = 1;
            goto end;
        }
        // According with rci_request value we do next step
        if (
            RCI_Request == -1 ||
            RCI_Request == -2 ||
            RCI_Request == -3 ||
            RCI_Request == -4 ||
            RCI_Request == -5 ||
            RCI_Request == -6
            )
        {
            // Exit rci cycle
            successful = 1;
        }
        if (RCI_Request == 1)
        {
            // Recalculate function value
            fdef(
                &nfo,       // in:     dimension of function value
                &nfi,       // in:     number of function variables
                x,          // in:     solution vector
                fvec        // out:    function value f(x)
                );
        }
        if (RCI_Request == 2)
        {
            // Compute jacobi matrix
            calculate_jacobi_mat(
                                    fdef,       // in:     external objective function
                                    nfi,        // in:     number of function variables
                                    nfo,        // in:     dimension of function value
                                    1,          // in:       matrix ordering
                                    fjac,       // out:    jacobi matrix
                                    x,          // in:     solution vector
                                    jac_eps     // in:     jacobi calculation precision
                                    );
        }
    }

    ////////////////////////////////////////////////////////////////
    /******************* Get solution status **********************/
    ////////////////////////////////////////////////////////////////
    res = dtrnlsp_get(
                    &handle,    // in:        TR solver handle
                    &iter,      // out:       number of iterations
                    &st_cr,     // out:       number of stop criterion
                    &r1,        // out:       initial residuals
                    &r2         // out:       final residuals
                    );
    if ( res != TR_SUCCESS )
    {
        SHOW_ERROR_MESAGE("error in dtrnlsp_get");
        error = 1;
        goto end;
    }


    ////////////////////////////////////////////////////////////////
    /*********************** Free memory **************************/
    ////////////////////////////////////////////////////////////////
    if (dtrnlsp_delete (&handle) != TR_SUCCESS)
    {
        /* if function does not complete successfully then print error message */
        SHOW_ERROR_MESAGE("Error in dtrnlsp_delete");
        error = 1;
        goto end;
    }

end:
    // Delete heap memory associated to the function evaluation
    // and the jacobian matrix
    mkl_free (fjac);
    mkl_free (fvec);

    /* Release internal memory that might be used for computations.          */
    /* NOTE: It is important to call the routine below to avoid memory leaks */
    /* unless you disable  Intel oneMKL Memory Manager                       */
    MKL_Free_Buffers ();
    if (error != 0)
    {
        exit(500);
    }
    if (mem_error == 1) 
    {
        SHOW_ERROR_MESAGE("Insufficient memory");
        exit(501);
    }

    ////////////////////////////////////////////////////////////////
    /********* Check residuals to validate solution ***************/
    ////////////////////////////////////////////////////////////////
    /* if final residual less then required precision then print pass */
    switch ( st_cr )
    {
        case 1:
            std::cerr << "Maximum number of iterations exceed! - STATUS: FAILED" << std::endl;
            tr_status = 1;
            break;

        case 2:
            if ( verbose )
            {
                std::cout << "Δ < eps(1) - STATUS: PASS" << std::endl;
            }
            tr_status = 0;
            break;
        
        case 3:
            if ( verbose )
            {
                std::cout << "||F(x)||2 < eps(2) - STATUS: PASS" << std::endl;
            }
            tr_status = 0;
            break;
        
        case 4:
            std::cerr << "||J(x)(1:m,j)||2 < eps(3) - STATUS: FAILED" << std::endl;
            tr_status = 1;
            break;
        
        case 5:
            if ( verbose )
            {
                std::cout << "||s||2 < eps(4) - STATUS: PASS" << std::endl;
            }
            tr_status = 0;
            break;
        
        case 6:
            if ( verbose )
            {
                std::cout << "||F(x)||2 - ||F(x) - J(x)s||2 < eps(5) - STATUS: PASS" << std::endl;
            }
            tr_status = 0;
            break;
    }
}

#endif