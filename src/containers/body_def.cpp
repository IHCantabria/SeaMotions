
// Include local modules
#include "body_def.hpp"


BodyDef::~BodyDef( void )
{
    if ( this->is_mesh )
    {
        delete this->mesh;
    }
}


void BodyDef::print( void )
{
    std::cout << "Mesh file name:   " << this->mesh_finame << std::endl;
    std::cout << "Mass:             " << this->mass << std::endl;
    std::cout << "COG_X:            " << this->cog[0] << std::endl;
    std::cout << "COG_Y:            " << this->cog[1] << std::endl;
    std::cout << "COG_Z:            " << this->cog[2] << std::endl;
    std::cout << "RXX:              " << this->rad_inertia[0] << std::endl; 
    std::cout << "RYY:              " << this->rad_inertia[1] << std::endl; 
    std::cout << "RZZ:              " << this->rad_inertia[2] << std::endl; 
}