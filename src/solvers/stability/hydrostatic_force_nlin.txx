
// Include local modules
#include "hydrostatic_force_nlin.hpp"
#include "../../math/math_tools.hpp"


template<typename T>
const cusfloat* HydrostaticForcesNLin<T>::get_last_hydrostatic_forces( 
                                                                                void
                                                                        )
{
    return this->_hyd_forces;
}


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
                                                    )
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

    /*  2.  Calculate hydrostatic forces and moments around the centre of 
            gravity*/
    
    //  2.1 Clean input RHS vector in order to acount for old spurious data
    cusfloat rhs[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    //  2.2 Calculate hydrostatic forces
    PanelGeom* paneli = nullptr;
    for ( int i=0; i<this->_mesh->get_elems_np( ); i++ )
    {
        paneli = this->_mesh->get_panel( i );
        if ( paneli->location_zone < 1 )
        {
            ( *this->_force_interf )( this->_mesh->get_panel( i ), rhs );
        }
    }

    // 2.3 Account for weight in the vertical force
    rhs[2] -= this->_weight;
    
    // 3.0 Copy forces intto internal storage to have latera access if needed
    for ( int i=0; i<6; i++ )
    {
        this->_hyd_forces[i] = rhs[i];
    }

    return rhs[2];
}