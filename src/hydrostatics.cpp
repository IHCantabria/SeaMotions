
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

// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "hydrostatics.hpp"
#include "./math/integration.hpp"
#include "./math/math_tools.hpp"
#include "tools.hpp"


void Hydrostatics::_calculate( 
                                Mesh*       mesh,
                                MpiConfig*  mpi_config
                            )
{
    /*********************************************/
    /**** Calculate volume and area integrals ****/
    /*********************************************/
    
    // Loop over panels to calculate total hydrostatic force
    // over the floating object
    cusfloat _volume        = 0.0;
    cusfloat _volume_mx     = 0.0;
    cusfloat _volume_my     = 0.0;
    cusfloat _volume_mz     = 0.0;
    cusfloat _area_wl       = 0.0;
    cusfloat _area_wl_mx    = 0.0;
    cusfloat _area_wl_my    = 0.0;
    cusfloat _area_wl_ixx   = 0.0;
    cusfloat _area_wl_ixy   = 0.0;
    cusfloat _area_wl_iyy   = 0.0;


    #ifndef MPI_BUILD
    int     start_elem      = mpi_config->proc_rank * 0;
    int     last_elem       = mesh->elems_np;
    #else
    int     elems_per_proc  = static_cast<int>( std::ceil( static_cast<cusfloat>( mesh->elems_np ) / static_cast<cusfloat>( mpi_config->procs_total ) ) );
    int     start_elem      = elems_per_proc * mpi_config->proc_rank;
    int     last_elem       = elems_per_proc * ( mpi_config->proc_rank + 1 );
    last_elem               = ( last_elem > mesh->elems_np ) ? mesh->elems_np: last_elem;
    #endif

    PanelGeom panel_proj;
    for ( int i=start_elem; i<last_elem; i++ )
    {
        // Get current panel
        PanelGeom* paneli   = mesh->panels[i];

        if ( paneli->type == 0 )
        {
            // Calculate water plane area
            auto    area_wl_fcn     =   [&]( int ){ return paneli->normal_vec[2]; };

            gauss2d_loop<NUM_GP>(
                                    _area_wl,
                                    area_wl_fcn,
                                    paneli
                                );

            // Calculate water plane area first order moments
            auto    area_wl_mx_fcn  =   [&]( int i )
                                        { 
                                            return (
                                                        paneli->gauss_points_global_y[i]
                                                        *
                                                        paneli->normal_vec[2]
                                                    ); 
                                        };

            gauss2d_loop<NUM_GP>(
                                    _area_wl_mx,
                                    area_wl_mx_fcn,
                                    paneli
                                );

            auto    area_wl_my_fcn  =   [&]( int i )
                                        { 
                                            return (
                                                        paneli->gauss_points_global_x[i]
                                                        *
                                                        paneli->normal_vec[2]
                                                    ); 
                                        };

            gauss2d_loop<NUM_GP>(
                                    _area_wl_my,
                                    area_wl_my_fcn,
                                    paneli
                                );

            // Calculate water plane area second order moments
            auto    area_wl_ixx_fcn  =   [&]( int i )
                                        { 
                                            return (
                                                        pow2s( paneli->gauss_points_global_y[i] )
                                                        *
                                                        paneli->normal_vec[2]
                                                    ); 
                                        };

            gauss2d_loop<NUM_GP>(
                                    _area_wl_ixx,
                                    area_wl_ixx_fcn,
                                    paneli
                                );

            auto    area_wl_iyy_fcn  =   [&]( int i )
                                        { 
                                            return (
                                                        pow2s( paneli->gauss_points_global_x[i] )
                                                        *
                                                        paneli->normal_vec[2]
                                                    ); 
                                        };

            gauss2d_loop<NUM_GP>(
                                    _area_wl_iyy,
                                    area_wl_iyy_fcn,
                                    paneli
                                );

            // Calculate volume
            auto    volume_fcn      =   [&](int i)
                                        { 
                                            return ( 
                                                        -paneli->gauss_points_global_z[i] 
                                                        * 
                                                        paneli->normal_vec[2] 
                                                    );
                                        };
            
            gauss2d_loop<NUM_GP>(
                                        _volume,
                                        volume_fcn,
                                        paneli
                                    );

            // Calculate volume moments
            auto    volume_mx_fcn   =   [&](int i)
                                        { 
                                            return (
                                                        -paneli->gauss_points_global_z[i] 
                                                        *  
                                                        paneli->gauss_points_global_x[i]
                                                        * 
                                                        paneli->normal_vec[2] 
                                                    ); 
                                        };
            
            gauss2d_loop<NUM_GP>(
                                        _volume_mx,
                                        volume_mx_fcn,
                                        paneli
                                    );

            auto    volume_my_fcn   =   [&](int i)
                                        { 
                                            return (
                                                        -paneli->gauss_points_global_z[i] 
                                                        *  
                                                        paneli->gauss_points_global_y[i]
                                                        * 
                                                        paneli->normal_vec[2] 
                                                    ); 
                                        };
            
            gauss2d_loop<NUM_GP>(
                                        _volume_my,
                                        volume_my_fcn,
                                        paneli
                                    );

            auto    volume_mz_fcn   =   [&](int i)
                                        { 
                                            return (
                                                        -paneli->gauss_points_global_z[i] 
                                                        *  
                                                        paneli->gauss_points_global_z[i] / 2.0
                                                        * 
                                                        paneli->normal_vec[2] 
                                                    ); 
                                        };
            
            gauss2d_loop<NUM_GP>(
                                        _volume_mz,
                                        volume_mz_fcn,
                                        paneli
                                    );
        }
    }

    #ifdef MPI_BUILD
    // Sum the volume in all processors
    MPI_Allreduce(
                    MPI_IN_PLACE,
                    &_volume,
                    1,
                    mpi_cusfloat,
                    MPI_SUM,
                    MPI_COMM_WORLD
                );

    // Sum all the volume X moment
    MPI_Allreduce(
                    MPI_IN_PLACE,
                    &_volume_mx,
                    1,
                    mpi_cusfloat,
                    MPI_SUM,
                    MPI_COMM_WORLD
                );

    // Sum all the volume Y moment
    MPI_Allreduce(
                    MPI_IN_PLACE,
                    &_volume_my,
                    1,
                    mpi_cusfloat,
                    MPI_SUM,
                    MPI_COMM_WORLD
                );

    // Sum all the volume Z moment
    MPI_Allreduce(
                    MPI_IN_PLACE,
                    &_volume_mz,
                    1,
                    mpi_cusfloat,
                    MPI_SUM,
                    MPI_COMM_WORLD
                );

    // Sum all the area
    MPI_Allreduce(
                    MPI_IN_PLACE,
                    &_area_wl,
                    1,
                    mpi_cusfloat,
                    MPI_SUM,
                    MPI_COMM_WORLD
                );

    // Sum all the area X moments
    MPI_Allreduce(
                    MPI_IN_PLACE,
                    &_area_wl_mx,
                    1,
                    mpi_cusfloat,
                    MPI_SUM,
                    MPI_COMM_WORLD
                );

    // Sum all the area Y moments
    MPI_Allreduce(
                    MPI_IN_PLACE,
                    &_area_wl_my,
                    1,
                    mpi_cusfloat,
                    MPI_SUM,
                    MPI_COMM_WORLD
                );

    // Sum all the area interia around X axis
    MPI_Allreduce(
                    MPI_IN_PLACE,
                    &_area_wl_ixx,
                    1,
                    mpi_cusfloat,
                    MPI_SUM,
                    MPI_COMM_WORLD
                );

    // Sum all the area interia in XY plane
    MPI_Allreduce(
                    MPI_IN_PLACE,
                    &_area_wl_ixy,
                    1,
                    mpi_cusfloat,
                    MPI_SUM,
                    MPI_COMM_WORLD
                );

    // Sum all the area interia around Y axis
    MPI_Allreduce(
                    MPI_IN_PLACE,
                    &_area_wl_iyy,
                    1,
                    mpi_cusfloat,
                    MPI_SUM,
                    MPI_COMM_WORLD
                );

    #endif

    // Calculate volume and displacement
    this->volume        = _volume;
    this->displacement  = _volume * rho_water;

    // Calculate KB
    this->cob[0]        = _volume_mx / volume;
    this->cob[1]        = _volume_my / volume;
    this->cob[2]        = _volume_mz / volume;
    this->kb            = this->cob[2] - mesh->z_min;

    // Storage area properties
    this->wl_area       = _area_wl;
    this->wl_area_mx    = _area_wl_mx;
    this->wl_area_my    = _area_wl_my;
    this->wl_area_ixx   = _area_wl_ixx;
    this->wl_area_ixy   = _area_wl_ixy;
    this->wl_area_iyy   = _area_wl_iyy;

    // Calculate water line area centre of gravity
    this->wl_area_cog[0] = this->wl_area_mx / this->wl_area;
    this->wl_area_cog[1] = this->wl_area_my / this->wl_area;
    this->wl_area_cog[2] = 0.0;

    /*********************************************/
    /***** Calculate hydrostatic properties ******/
    /*********************************************/

    // Calculate KG
    this->kg  = this->cog[2] - mesh->z_min;

    // Calculate metacentric radius
    this->bmx = this->wl_area_ixx / this->volume;
    this->bmy = this->wl_area_iyy / this->volume;

    // Calculate metacentric point
    this->kmx = this->kb  + this->bmx;
    this->kmy = this->kb  + this->bmy;

    // Calculate metacentric height
    this->gmx = this->kmx - this->kg;
    this->gmy = this->kmy - this->kg;

    /*********************************************/
    /******* Calculate hydrostatic matrix ********/
    /*********************************************/

    // Clean hydrostatic stiffness matrix
    clear_vector( 36, this->hydstiffmat );

    // K33 - Heave
    this->hydstiffmat[14] = this->grav_acc * this->rho_water * wl_area;

    // K34 - K43 - Heave/Roll
    this->hydstiffmat[15] = this->grav_acc * this->rho_water * this->wl_area_my;
    this->hydstiffmat[20] = this->grav_acc * this->rho_water * this->wl_area_my;

    // K35 - K53 - Heave/Pitch
    this->hydstiffmat[16] = - this->grav_acc * this->rho_water * this->wl_area_mx;
    this->hydstiffmat[26] = - this->grav_acc * this->rho_water * this->wl_area_mx;

    // K44 - Roll
    this->hydstiffmat[21] = this->grav_acc * this->rho_water * ( this->wl_area_ixx + this->volume * ( this->cob[2] - this->cog[2] ) );

    // K45 - K54 - Roll/Pitch
    this->hydstiffmat[22] = this->grav_acc * this->rho_water * this->wl_area_ixy;
    this->hydstiffmat[27] = this->grav_acc * this->rho_water * this->wl_area_ixy;

    // K55 - Pitch
    this->hydstiffmat[28] = this->grav_acc * this->rho_water * ( this->wl_area_iyy + this->volume * ( this->cob[2] - this->cog[2] ) );

    // K46 - K56 - Yaw
    this->hydstiffmat[23] = - this->grav_acc * this->rho_water * ( this->cob[0] - this->cog[0] ) * this->volume;
    this->hydstiffmat[29] = - this->grav_acc * this->rho_water * ( this->cob[1] - this->cog[1] ) * this->volume;

}


Hydrostatics::Hydrostatics( 
                                Mesh*       mesh,
                                cusfloat    rhow_in,
                                cusfloat    grav_acc_in,
                                cusfloat    mass_in,
                                cusfloat*   cog_in,
                                cusfloat*   rad_inertia_in,
                                MpiConfig*  mpi_config
                            )
{
    // Storage required input data
    this->grav_acc  = grav_acc_in;
    this->rho_water = rhow_in;
    this->mass      = mass_in;
    
    copy_vector( 3, cog_in, this->cog );
    copy_vector( 3, rad_inertia_in, this->rad_inertia );

    // Calculate hydrostatics
    this->_calculate( 
                        mesh,
                        mpi_config
                    );
}


void Hydrostatics::print( void )
{
    std::cout << std::fixed << std::setprecision( 6 );
    std::cout << std::endl;
    std::cout << "HYDROSTATICS: " << std::endl;
    std::cout << " - Volume properties: "             << std::endl;
    std::cout << "--------------------------------"   << std::endl;
    std::cout << "      + Displaced Volume [m3]:    " << this->volume << std::endl;
    std::cout << "      + Displacement [kg]:        " << this->displacement << std::endl;
    std::cout << "      + COB [m]:                  " << this->cob[0] << ", " << this->cob[1] << ", " << this->cob[2] << std::endl;
    std::cout << " - Area properties: " << std::endl;
    std::cout << "--------------------------------"   << std::endl;
    std::cout << "      + WL Area [m2]:             " << this->wl_area << std::endl;
    std::cout << "      + WL Area Mx [m3]:          " << this->wl_area_mx << std::endl;
    std::cout << "      + WL Area My [m3]:          " << this->wl_area_my << std::endl;
    std::cout << "      + WL area Ixx [m4]:         " << this->wl_area_ixx << std::endl;
    std::cout << "      + WL area Iyy [m4]:         " << this->wl_area_iyy << std::endl;
    std::cout << " - Mass Properties: "               << std::endl;
    std::cout << "--------------------------------"   << std::endl;
    std::cout << "      + Structural Mass [kg]:     " << this->mass << std::endl;
    std::cout << "      + COG [m]:                  " << this->cog[0] << ", " << this->cog[1] << ", " << this->cog[2] << std::endl;
    std::cout << " - Stability Properties: "          << std::endl;
    std::cout << "--------------------------------"   << std::endl;
    std::cout << "      + KB [m]:                   " << this->kb << std::endl;
    std::cout << "      + BMX [m]:                  " << this->bmx << std::endl;
    std::cout << "      + BMY [m]:                  " << this->bmy << std::endl;
    std::cout << "      + KMX [m]:                  " << this->kmx << std::endl;
    std::cout << "      + KMY [m]:                  " << this->kmy << std::endl;
    std::cout << "      + KG [m]:                   " << this->kg << std::endl;
    std::cout << "      + GMX [m]:                  " << this->gmx << std::endl;
    std::cout << "      + GMY [m]:                  " << this->gmy << std::endl;
    std::cout << " - Hydrostatic Stiffness Matrix: "  << std::endl;
    std::cout << "--------------------------------"   << std::endl;
    print_matrix( 6, 6, this->hydstiffmat, 6, 0, 1 );
}