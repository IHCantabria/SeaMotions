
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
    
    // Define lambda function to calculate the 
    // volume and the volume centroid
    auto get_z_coord     =  [ ]
                            ( PanelGeom* panel, cusfloat x, cusfloat y ) -> cusfloat
                            {
                                // Get local coordinates
                                cusfloat xi=0.0, eta=0.0;
                                panel->local_coords_from_z_proj( x, y, xi, eta );

                                // Calculate Z position over the mesh panel
                                cusfloat global_pos[3] = { 0.0, 0.0, 0.0 };
                                panel->local_to_global( xi, eta, global_pos );

                                return global_pos[2];
                            };
    
    // Define lambda function to calculate wl area,
    // area centroid and area inertias
    auto wl_area_fcn        =   [ ]
                                ( cusfloat, cusfloat, cusfloat , cusfloat , cusfloat ) -> cuscomplex
                                {
                                    return 1.0;
                                };

    auto wl_area_mom_x_fcn  =   [ this ]
                                ( cusfloat, cusfloat, cusfloat x, cusfloat , cusfloat ) -> cuscomplex
                                {
                                    return x - this->cog[0];
                                };

    auto wl_area_mom_y_fcn  =   [ this ]
                                ( cusfloat, cusfloat, cusfloat , cusfloat y, cusfloat ) -> cuscomplex
                                {
                                    return y - this->cog[1];
                                };

    auto wl_area_ixx_fcn    =   [ this ]
                                ( cusfloat, cusfloat, cusfloat , cusfloat y, cusfloat ) -> cuscomplex
                                {
                                    return pow2s( y - this->cog[1] );
                                };

    auto wl_area_ixy_fcn    =   [ this ]
                                ( cusfloat, cusfloat, cusfloat x, cusfloat y, cusfloat ) -> cuscomplex
                                {
                                    return ( x - this->cog[0] ) * ( y - this->cog[1] );
                                };
    
    auto wl_area_iyy_fcn    =   [ this ]
                                ( cusfloat, cusfloat, cusfloat x, cusfloat , cusfloat ) -> cuscomplex
                                {
                                    return pow2s( x - this->cog[0] );
                                };

    // Loop over panels to calculate total hydrostatic force
    // over the floating object
    cusfloat _ai            = 0.0;
    cusfloat _amxi          = 0.0;
    cusfloat _amyi          = 0.0;
    cusfloat _aixi          = 0.0;
    cusfloat _aixyi         = 0.0;
    cusfloat _aiyi          = 0.0;
    cusfloat _area_abs_eps  = 0.1;
    cusfloat _area_rel_eps  = 0.0;
    cusfloat normal_sign    = 0.0;
    cusfloat vi             = 0.0;
    cusfloat vim_x          = 0.0;
    cusfloat vim_y          = 0.0;
    cusfloat vim_z          = 0.0;
    cusfloat _vol_abs_eps   = 0.1;
    cusfloat _vol_rel_eps   = 0.0;
    cusfloat _volume        = 0.0;
    cusfloat _volume_mom_x  = 0.0;
    cusfloat _volume_mom_y  = 0.0;
    cusfloat _volume_mom_z  = 0.0;
    cusfloat _wl_area       = 0.0;
    cusfloat _wl_area_mx    = 0.0;
    cusfloat _wl_area_my    = 0.0;
    cusfloat _wl_area_ixx   = 0.0;
    cusfloat _wl_area_ixy   = 0.0;
    cusfloat _wl_area_iyy   = 0.0;


    #ifndef MPI_BUILD
    int     start_elem      = mpi_config->proc_rank * 0;
    int     last_elem       = mesh->elems_np;
    #else
    int     elems_per_proc  = static_cast<int>( std::ceil( static_cast<cusfloat>( mesh->elems_np ) / static_cast<cusfloat>( mpi_config->procs_total ) ) );
    int     start_elem      = elems_per_proc * mpi_config->proc_rank;
    int     last_elem       = elems_per_proc * ( mpi_config->proc_rank + 1 );
    last_elem               = ( last_elem > mesh->elems_np ) ? mesh->elems_np: last_elem;
    #endif

    for ( int i=start_elem; i<last_elem; i++ )
    {
        // Get current panel
        PanelGeom* panel    = mesh->panels[i];

        if ( panel->type == 0 )
        {
            // Get projection panel
            PanelGeom* panel_proj = new PanelGeom;
            panel->get_panel_xy_proj( panel_proj );

            // Get normal sign
            normal_sign = (cusfloat)sign( mesh->panels[i]->normal_vec[2] );

            // std::cout << "Panel: " << i << " - nz: " << mesh->panels[i]->normal_vec[2] << std::endl;

            if ( std::abs( mesh->panels[i]->normal_vec[2] ) > 1e-2 )
            {
                // Define volume functions
                auto volume_fcn         =   [ panel, get_z_coord ]
                                            ( cusfloat, cusfloat, cusfloat x, cusfloat y, cusfloat ) -> cuscomplex
                                            {
                                                // Calculate Z position over the mesh panel
                                                cusfloat zg = get_z_coord( panel, x, y );

                                                return cuscomplex( zg, 0.0 );
                                            };

                auto volume_mom_x_fcn   =   [ panel, get_z_coord ]
                                            ( cusfloat, cusfloat, cusfloat x, cusfloat y, cusfloat ) -> cuscomplex
                                            {
                                                // Calculate Z position over the mesh panel
                                                cusfloat zg = get_z_coord( panel, x, y );

                                                return cuscomplex( x*zg, 0.0 );
                                            };

                auto volume_mom_y_fcn   =   [ panel, get_z_coord ]
                                            ( cusfloat, cusfloat, cusfloat x, cusfloat y, cusfloat ) -> cuscomplex
                                            {
                                                // Calculate Z position over the mesh panel
                                                cusfloat zg = get_z_coord( panel, x, y );

                                                return cuscomplex( y*zg, 0.0 );
                                            };

                auto volume_mom_z_fcn   =   [ panel, get_z_coord ]
                                            ( cusfloat, cusfloat, cusfloat x, cusfloat y, cusfloat ) -> cuscomplex
                                            {
                                                // Calculate Z position over the mesh panel
                                                cusfloat zg = get_z_coord( panel, x, y );

                                                return cuscomplex( zg*zg/2.0, 0.0 );
                                            };
                
                // Calculate submerged volume
                vi                  = adaptive_quadrature_panel(
                                                                    panel_proj,
                                                                    volume_fcn,
                                                                    _vol_abs_eps,
                                                                    _vol_rel_eps
                                                                ).real( );
                _volume             +=  vi * normal_sign;

                // throw std::runtime_error( "" );

                // Calculate x moment submerged volume
                vim_x               = adaptive_quadrature_panel(
                                                                    panel_proj,
                                                                    volume_mom_x_fcn,
                                                                    _vol_abs_eps,
                                                                    _vol_rel_eps
                                                                ).real( );
                _volume_mom_x       +=  vim_x * normal_sign;

                // Calculate x moment submerged volume
                vim_y               = adaptive_quadrature_panel(
                                                                    panel_proj,
                                                                    volume_mom_y_fcn,
                                                                    _vol_abs_eps,
                                                                    _vol_rel_eps
                                                                ).real( );
                _volume_mom_y       +=  vim_y * normal_sign;

                // Calculate z moment submerged volume
                vim_z               = adaptive_quadrature_panel(
                                                                    panel_proj,
                                                                    volume_mom_z_fcn,
                                                                    _vol_abs_eps,
                                                                    _vol_rel_eps
                                                                ).real( );
                _volume_mom_z       +=  vim_z * normal_sign;

                // Integrate panel area
                _ai                  = adaptive_quadrature_panel(
                                                                    panel_proj,
                                                                    wl_area_fcn,
                                                                    _area_abs_eps,
                                                                    _area_rel_eps
                                                                ).real( );
                _wl_area            -= _ai * normal_sign;

                // Integrate panel area moments around x axis
                _amxi                = adaptive_quadrature_panel(
                                                                    panel_proj,
                                                                    wl_area_mom_x_fcn,
                                                                    _area_abs_eps,
                                                                    _area_rel_eps
                                                                ).real( );
                _wl_area_mx         -= _amxi * normal_sign;

                // Integrate panel area moments around y axis
                _amyi                = adaptive_quadrature_panel(
                                                                    panel_proj,
                                                                    wl_area_mom_y_fcn,
                                                                    _area_abs_eps,
                                                                    _area_rel_eps
                                                                ).real( );
                _wl_area_my         -= _amyi * normal_sign;

                // Integrate panel area inertia around x axis
                _aixi                = adaptive_quadrature_panel(
                                                                    panel_proj,
                                                                    wl_area_ixx_fcn,
                                                                    _area_abs_eps,
                                                                    _area_rel_eps
                                                                ).real( );
                _wl_area_ixx        -= _aixi * normal_sign;

                // Integrate panel area inertia around x axis
                _aixyi               = adaptive_quadrature_panel(
                                                                    panel_proj,
                                                                    wl_area_ixy_fcn,
                                                                    _area_abs_eps,
                                                                    _area_rel_eps
                                                                ).real( );
                _wl_area_ixy        -= _aixyi * normal_sign;

                // Integrate panel area inertia around x axis
                _aiyi                = adaptive_quadrature_panel(
                                                                    panel_proj,
                                                                    wl_area_iyy_fcn,
                                                                    _area_abs_eps,
                                                                    _area_rel_eps
                                                                ).real( );
                _wl_area_iyy        -= _aiyi * normal_sign;

            }
            else
            {
                // Calculate submerged volume
                vi                  = panel_proj->area * panel->center[2];
                _volume             +=  vi * normal_sign;

                // throw std::runtime_error( "" );

                // Calculate x moment submerged volume
                vim_x               = vi * panel->center[0];
                _volume_mom_x       +=  vim_x * normal_sign;

                // Calculate x moment submerged volume
                vim_y               = vi * panel->center[1];
                _volume_mom_y       +=  vim_y * normal_sign;

                // Calculate z moment submerged volume
                vim_z               = vi * panel->center[2] / 2.0;
                _volume_mom_z       +=  vim_z * normal_sign;

                // Integrate panel area
                _ai                  = panel_proj->area;
                _wl_area            -= _ai * normal_sign;

                // Integrate panel area moments around x axis
                _amxi                = _ai * ( panel->center[0] - this->cog[0] );
                _wl_area_mx         -= _amxi * normal_sign;

                // Integrate panel area moments around y axis
                _amyi                = _ai * ( panel->center[1] - this->cog[1] );
                _wl_area_my         -= _amyi * normal_sign;

                // Integrate panel area inertia around x axis
                _aixi                = _ai * pow2s( panel->center[1] - this->cog[1] );
                _wl_area_ixx        -= _aixi * normal_sign;

                // Integrate panel area inertia in the XY plane
                _aixyi               = _ai * ( panel->center[0] - this->cog[0] ) * ( panel->center[1] - this->cog[1] );
                _wl_area_ixy        -= _aixyi * normal_sign;

                // Integrate panel area inertia around x axis
                _aiyi                = _ai * pow2s( panel->center[0] - this->cog[0] );
                _wl_area_iyy        -= _aiyi * normal_sign;

            }
        }
    }

    #ifdef MPI_BUILD
    // Sum the volume in all processors
    cusfloat _volume_d = 0.0;
    MPI_Allreduce(
                    &_volume,
                    &_volume_d,
                    1,
                    mpi_cusfloat,
                    MPI_SUM,
                    MPI_COMM_WORLD
                );
    _volume = _volume_d;

    // Sum all the volume X moment
    cusfloat _volume_mom_x_d = 0.0;
    MPI_Allreduce(
                    &_volume_mom_x,
                    &_volume_mom_x_d,
                    1,
                    mpi_cusfloat,
                    MPI_SUM,
                    MPI_COMM_WORLD
                );
    _volume_mom_x   = _volume_mom_x_d;

    // Sum all the volume Y moment
    cusfloat _volume_mom_y_d = 0.0;
    MPI_Allreduce(
                    &_volume_mom_y,
                    &_volume_mom_y_d,
                    1,
                    mpi_cusfloat,
                    MPI_SUM,
                    MPI_COMM_WORLD
                );
    _volume_mom_y   = _volume_mom_y_d;

    // Sum all the volume Z moment
    cusfloat _volume_mom_z_d = 0.0;
    MPI_Allreduce(
                    &_volume_mom_z,
                    &_volume_mom_z_d,
                    1,
                    mpi_cusfloat,
                    MPI_SUM,
                    MPI_COMM_WORLD
                );
    _volume_mom_z   = _volume_mom_z_d;

    // Sum all the area
    cusfloat _wl_area_d = 0.0;
    MPI_Allreduce(
                    &_wl_area,
                    &_wl_area_d,
                    1,
                    mpi_cusfloat,
                    MPI_SUM,
                    MPI_COMM_WORLD
                );
    _wl_area        = _wl_area_d;

    // Sum all the area X moments
    cusfloat _wl_area_mx_d = 0.0;
    MPI_Allreduce(
                    &_wl_area_mx,
                    &_wl_area_mx_d,
                    1,
                    mpi_cusfloat,
                    MPI_SUM,
                    MPI_COMM_WORLD
                );
    _wl_area_mx     = _wl_area_mx_d;

    // Sum all the area Y moments
    cusfloat _wl_area_my_d = 0.0;
    MPI_Allreduce(
                    &_wl_area_my,
                    &_wl_area_my_d,
                    1,
                    mpi_cusfloat,
                    MPI_SUM,
                    MPI_COMM_WORLD
                );
    _wl_area_my     = _wl_area_my_d;

    // Sum all the area interia around X axis
    cusfloat _wl_area_ixx_d = 0.0;
    MPI_Allreduce(
                    &_wl_area_ixx,
                    &_wl_area_ixx_d,
                    1,
                    mpi_cusfloat,
                    MPI_SUM,
                    MPI_COMM_WORLD
                );
    _wl_area_ixx    = _wl_area_ixx_d;

    // Sum all the area interia in XY plane
    cusfloat _wl_area_ixy_d = 0.0;
    MPI_Allreduce(
                    &_wl_area_ixy,
                    &_wl_area_ixy_d,
                    1,
                    mpi_cusfloat,
                    MPI_SUM,
                    MPI_COMM_WORLD
                );
    _wl_area_ixy    = _wl_area_ixy_d;

    // Sum all the area interia around Y axis
    cusfloat _wl_area_iyy_d = 0.0;
    MPI_Allreduce(
                    &_wl_area_iyy,
                    &_wl_area_iyy_d,
                    1,
                    mpi_cusfloat,
                    MPI_SUM,
                    MPI_COMM_WORLD
                );
    _wl_area_iyy    = _wl_area_iyy_d;

    #endif

    // Calculate volume and displacement
    this->volume        = _volume;
    this->displacement  = _volume * rho_water;

    // Calculate KB
    this->cob[0]        = _volume_mom_x / volume;
    this->cob[1]        = _volume_mom_y / volume;
    this->cob[2]        = _volume_mom_z / volume;
    this->kb            = this->cob[2] - mesh->z_min;

    // Storage area properties
    this->wl_area       = _wl_area;
    this->wl_area_mx    = _wl_area_mx;
    this->wl_area_my    = _wl_area_my;
    this->wl_area_ixx   = _wl_area_ixx;
    this->wl_area_ixy   = _wl_area_ixy;
    this->wl_area_iyy   = _wl_area_iyy;

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