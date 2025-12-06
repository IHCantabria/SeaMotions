
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
#include "math_tools.hpp"
#include "./sparse/sparse_math.hpp"

//--------------------------------------------------------------------
//-- Define class auxiliar functions
//--------------------------------------------------------------------
/**
 * @brief   Function used to unify kinematics inconsistency message format
 * 
 * @param   Index of the degree of freedom with the inconsistency
 * @param   Value for the inconsistency
 * @para    Field name where the inconsistency was detected
 */
template<typename T>
inline  void    _print_inconsistency_msg( 
                                                int                 index,
                                                const T&            value,
                                                const std::string   quantity_name
                                            )
{
    std::cerr << "[WARN]: Newmark-Beta solver detected inconsistency between "
                 "restricted degrees of freedom and initial boundary conditions kinematics.\n";
    std::cerr << "\t\t-> " << quantity_name << "[ " << index << " ]: " << value << "\n";
}          

//--------------------------------------------------------------------
//-- Define Newmark-Beta solver
//--------------------------------------------------------------------

/**
 * @brief Newmark-Beta time solver implementation with constant time step.
 *
 * @tparam Functor matching the expected interface that calculates the external force vector
 */
template<typename T>
class NewmarkBeta
{
private:
    /* Declare private class attributes */
    bool    _is_restrictions    = false;        // Switch to check if restrictions vector was internally allocated
    int*    _restrictions       = nullptr;      // Vector composed of 0 or 1 to apply restrictions to the selected degrees of freedom

    /* Declare private class methods */
    
    /**
     * @brief   Apply restrictions to a given vector.
     * 
     * This methods pretends to unify the imposition of restrictions for kinematics and 
     * external forces through a unique interface
     * 
     * @param   Vector to restrict
     * 
     */
    void    _apply_restrictions( cusfloat* _vec );

    /**
     * @brief It is the delegate of the class constructor
     * 
     * This class is used to build the object so the constructors 
     * calls are used as interfaces so it is possible to input 
     * with different parameters and to have some default states
     * 
     * @tparam  Functor matching the expected interface that calculates the external force vector
     * @param   Mass matrix in CSRMatrix format
     * @param   Stiffness matrix in CSRMatrix format
     * @param   Linear damping matrix in CSRMatrix format
     * @param   Time step to be used for time integration
     * @param   Start time for the simulation
     * @param   Initial position for the degrees of freedom described by CSRMatrixes
     * @param   Initial velocity for the degrees of freedom described by CSRMatrixes
     * @param   Initial acceleration for the degrees of freedom described by CSRMatrixes
     * @param   Restrictions for the degrees of freedom described by CSRMatrixes
     */
    void    _build( 
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
                    );

    /**
     * @brief Check input kinematic arguments to be aligned with the imposed restrictions
     */
    void    _check_init_kinematics_retrictions(
                                                void
                                            );

public:
    // Define class attributes
    cusfloat            beta            = 0.0;
    int                 count           = 0;
    sparse_matrix_t*    damp_mat        = nullptr;
    T                   fext_fcn          ;
    cusfloat            gamma           = 0.0;
    CSRMatrix*          lhs_csr         = nullptr;
    SparseSolver*       linsolve        = nullptr;
    sparse_matrix_t*    mass_mat        = nullptr;
    int                 mat_size        = 0;
    sparse_matrix_t*    stiff_damp_acc  = nullptr;
    sparse_matrix_t*    stiff_damp_vel  = nullptr;
    sparse_matrix_t*    stiff_mat       = nullptr;
    cusfloat*           rhs             = nullptr;      
    int                 rows_np         = 0;
    cusfloat            time            = 0.0;
    cusfloat            time_init       = 0.0;
    cusfloat            time_old        = 0.0;
    cusfloat            time_step       = 0.0;
    cusfloat*           _u_beta         = nullptr;
    cusfloat*           _u_gamma        = nullptr;
    cusfloat*           y0_acc          = nullptr;
    cusfloat*           y0_pos          = nullptr;
    cusfloat*           y0_vel          = nullptr;
    cusfloat*           y_acc           = nullptr;
    cusfloat*           y_pos           = nullptr;
    cusfloat*           y_vel           = nullptr;
    cusfloat*           y_acc_old       = nullptr;
    cusfloat*           y_pos_old       = nullptr;
    cusfloat*           y_vel_old       = nullptr;

    /* Declare class constructors and destructor */

    /**
     * @brief   Class constructor without imposing restrictions
     * 
     * @tparam  Functor matching the expected interface that calculates the external force vector
     * @param   Mass matrix in CSRMatrix format
     * @param   Stiffness matrix in CSRMatrix format
     * @param   Linear damping matrix in CSRMatrix format
     * @param   Time step to be used for time integration
     * @param   Start time for the simulation
     * @param   Initial position for the degrees of freedom described by CSRMatrixes
     * @param   Initial velocity for the degrees of freedom described by CSRMatrixes
     * @param   Initial acceleration for the degrees of freedom described by CSRMatrixes
     */
    NewmarkBeta(
                                        T               fext_in,
                                        CSRMatrix*      mass_mat_in,
                                        CSRMatrix*      stiff_mat_in,
                                        CSRMatrix*      damp_mat_in,
                                        cusfloat        time_step_in,
                                        cusfloat        t0_in,
                                        cusfloat*       y0_pos_in,
                                        cusfloat*       y0_vel_in,
                                        cusfloat*       y0_acc_in
                );

    /**
     * @brief   Class constructor imposing restrictions to the required DOFs
     * 
     * @tparam  Functor matching the expected interface that calculates the external force vector
     * @param   Mass matrix in CSRMatrix format
     * @param   Stiffness matrix in CSRMatrix format
     * @param   Linear damping matrix in CSRMatrix format
     * @param   Time step to be used for time integration
     * @param   Start time for the simulation
     * @param   Initial position for the degrees of freedom described by CSRMatrixes
     * @param   Initial velocity for the degrees of freedom described by CSRMatrixes
     * @param   Initial acceleration for the degrees of freedom described by CSRMatrixes
     * @param   Restrictions for the degrees of freedom described by CSRMatrixes
     */
    NewmarkBeta(
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
                );
    
    ~NewmarkBeta( void );

    // Define class methods
    void    get_values_at(
                                        cusfloat        time,
                                        cusfloat*       y_pos_at,
                                        cusfloat*       y_vel_at,
                                        cusfloat*       y_acc_at
                        );
    
    void    initialize(
                                        void
                        );
    
    void    interpolate_solution(
                                        cusfloat        time, 
                                        cusfloat*       y_pos_itp,
                                        cusfloat*       y_vel_itp,
                                        cusfloat*       y_acc_itp
                                );

    void step(
                                        void
            );

};


#include "newmark.txx"
