
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
    std::cout << "Mesh file name: " << this->mesh_finame << std::endl;
    std::cout << "Mass: " << this->mass << std::endl;
    std::cout << "COG_X: " << this->cog_x << std::endl;
    std::cout << "COG_Y: " << this->cog_y << std::endl;
    std::cout << "COG_Z: " << this->cog_z << std::endl;
    std::cout << "RXX: " << this->rxx << std::endl; 
    std::cout << "RYY: " << this->ryy << std::endl; 
    std::cout << "RZZ: " << this->rzz << std::endl; 
}