
// Include general usage libraries
#include <iostream>

// Include local modules
#include "hydrostatics.hpp"
#include "../math/integration.hpp"
#include "tools.hpp"


void Hydrostatics::_calculate( 
                                Mesh*       mesh
                            )
{
    /*********************************************/
    /****** Calculate volumetric properties ******/
    /*********************************************/
    
    // Define lambda function to calculate the 
    // hydrostatic pressure over each panel
    // auto hydrostat_fcn =        [this]
    //                             ( cusfloat , cusfloat , cusfloat z ) -> cuscomplex
    //                             {
    //                                 // std::cout << "LMB - z: " << z << std::endl;
    //                                 cusfloat hr =  hydrostatic_pressure(
    //                                                                         this->rho_water,
    //                                                                         this->grav_acc,
    //                                                                         z
    //                                                                     );

    //                                 return cuscomplex( hr, 0.0 );
    //                             };

    auto volume_fcn =           [ ]
                                ( cusfloat , cusfloat , cusfloat z ) -> cuscomplex
                                {
                                    // cusfloat hr =  hydrostatic_pressure_mom(
                                    //                                             this->rho_water,
                                    //                                             this->grav_acc,
                                    //                                             z
                                    //                                         );

                                    return cuscomplex( z, 0.0 );
                                };
    
    auto volume_mom_fcn =       [ ]
                                ( cusfloat , cusfloat , cusfloat z ) -> cuscomplex
                                {
                                    // cusfloat hr =  hydrostatic_pressure_mom(
                                    //                                             this->rho_water,
                                    //                                             this->grav_acc,
                                    //                                             z
                                    //                                         );

                                    return cuscomplex( z * z / 2.0, 0.0 );
                                };

    // auto hydrostat_mom_fcn =    [this]
    //                             ( cusfloat , cusfloat , cusfloat z ) -> cuscomplex
    //                             {
    //                                 // cusfloat hr =  hydrostatic_pressure_mom(
    //                                 //                                             this->rho_water,
    //                                 //                                             this->grav_acc,
    //                                 //                                             z
    //                                 //                                         );

    //                                 return cuscomplex( z, 0.0 );
    //                             };
    
    // Loop over panels to calculate total hydrostatic force
    // over the floating object
    // cusfloat fi             = 0.0;
    // cusfloat fim            = 0.0;
    cusfloat force_mag      = 0.0;
    cusfloat force_mag_mom  = 0.0;
    cusfloat force_z        = 0.0;
    cusfloat force_z_mom    = 0.0;
    cusfloat vi             = 0.0;
    cusfloat vim            = 0.0;
    cusfloat volume         = 0.0;
    cusfloat volume_mom     = 0.0;
    for ( int i=0; i<mesh->elems_np; i++ )
    {
        // int node_num = 0;
        // for ( int j=0; j<mesh->elems[i*mesh->enrl]; j++ )
        // {
        //     node_num = mesh->elems[i*mesh->enrl+j+1];
        //     std::cout << "z[" << j << "]: " << mesh->z[node_num] << std::endl;
        // }

        // Calculate submerged volume
        vi              = adaptive_quadrature_panel(
                                                        mesh->panels[i],
                                                        volume_fcn,
                                                        1e-6,
                                                        5
                                                    ).real( );
        volume          +=  vi * mesh->panels[i]->normal_vec[2];

        // Calculate moment submerged volume
        vim             = adaptive_quadrature_panel(
                                                        mesh->panels[i],
                                                        volume_mom_fcn,
                                                        1e-6,
                                                        5
                                                    ).real( );
        volume_mom      +=  vim * mesh->panels[i]->normal_vec[2];

        // // Calculate vertical force over the panel
        // fi              = adaptive_quadrature_panel(
        //                                                 mesh->panels[i],
        //                                                 hydrostat_fcn,
        //                                                 1e-6,
        //                                                 5
        //                                             ).real( );

        // force_mag       +=  fi;
        // force_z         +=  fi * mesh->panels[i]->normal_vec[2];
        // std::cout << "i: " << i << " - fi: " << fi << " - nz:" << mesh->panels[i]->normal_vec[2] << std::endl;

        // // Calculate moments with respect to the free surface
        // fim             = adaptive_quadrature_panel(
        //                                                 mesh->panels[i],
        //                                                 hydrostat_mom_fcn,
        //                                                 1e-6,
        //                                                 5
        //                                             ).real( );

        // force_mag_mom   += fim;
        // force_z_mom     += fim * mesh->panels[i]->normal_vec[2];

        // break;
    }

    // Calculate volume and displacement
    this->volume        = volume;
    this->displacement  = volume * rho_water;

    // Calculate KB
    std::cout << "force_z: " << force_z << std::endl;
    std::cout << "force_z_mom: " << force_z_mom << std::endl;
    std::cout << "force_mag: " << force_mag << std::endl;
    std::cout << "force_mag_mom: " << force_mag_mom << std::endl;
    std::cout << "volume: " << volume << std::endl;
    std::cout << "volume momment: " << volume_mom << std::endl;
    std::cout << "z_min: " << mesh->z_min << std::endl;
    this->kb            = volume_mom / volume - mesh->z_min;

}


Hydrostatics::Hydrostatics( 
                                Mesh*       mesh,
                                cusfloat    rhow,
                                cusfloat    grav_acc
                            )
{
    // Storage required input data
    this->grav_acc  = grav_acc;
    this->rho_water = rhow;

    // Calculate hydrostatics
    this->_calculate( 
                        mesh
                    );
}


void Hydrostatics::print( void )
{
    std::cout << std::fixed << std::setprecision( 6 );
    std::cout << "HYDROSTATICS: " << std::endl;
    std::cout << " - Displaced Volume [m3]: " << this->volume << std::endl;
    std::cout << " - Displacement [kg]:     " << this->displacement << std::endl;
    std::cout << " - KB [m]:                " << this->kb << std::endl;
}