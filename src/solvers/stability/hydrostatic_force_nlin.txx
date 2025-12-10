
// Include local modules
#include "hydrostatic_force_nlin.hpp"


template<typename T>
HydrostaticForcesNLin<T>::HydrostaticForcesNLin( 
                                                    RigidBodyMesh*  mesh_in,
                                                    cusfloat*       pos_in,
                                                    T               force_interf_in,
                                                    cusfloat        weight_in
                                                )
{
    // Storage input arguments
    this->_force_interf = force_interf_in;
    this->_mesh         = mesh_in;
    this->_weight       = weight_in;

    copy_vector( 6, pos_in, this->_pos_init );
}


template<typename T>
cusfloat    HydrostaticForcesNLin<T>::operator()   ( 
                                                        cusfloat heave
                                                    ) const
{
    /*  1.  Move and rotate mesh around centre of gravity to 
            be the state predicted by the stepper */
    this->_mesh->move( 
                        this->_pos_init[0], 
                        this->_pos_init[1], 
                        heave, 
                        this->_pos_init[3], 
                        this->_pos_init[4], 
                        this->_pos_init[5] 
                    );

    /*  2.  Cut mesh along the free surface */
    this->_mesh->check_underwater_panels( );

    /*  3.  Calculate hydrostatic forces and moments around the centre of 
            gravity*/
    
    //  3.1 Clean input RHS vector in order to acount for old spurious data
    cusfloat rhs[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    //  3.2 Calculate hydrostatic forces
    PanelGeom* paneli = nullptr;
    int count_normal = 0;
    cusfloat rprev = 0.0;
    for ( int i=0; i<this->_mesh->get_elems_np( ); i++ )
    {
        paneli = this->_mesh->get_panel( i );
        if ( paneli->location_zone == -1 )
        {
            ( *this->_force_interf )( this->_mesh->get_panel( i ), rhs );
        }
    }

    rhs[2] -= this->_weight;

    return rhs[2];
}