
#ifndef __body_hpp
#define __body_hpp

// Include general usage libraries
#include <string>

// Include local modules
#include "../mesh/mesh.hpp"


struct BodyDef
{
public:
    // Define class attributes
    cusfloat    cog_x           = 0.0;
    cusfloat    cog_y           = 0.0;
    cusfloat    cog_z           = 0.0;
    bool        is_mesh         = false;
    cusfloat    mass            = 0.0;
    Mesh*       mesh            = nullptr;
    std::string mesh_finame     = "";
    cusfloat    rxx             = 0.0;
    cusfloat    ryy             = 0.0;
    cusfloat    rzz             = 0.0;

    // Define class constructor and destructor
    ~BodyDef( void );

    // Define class methods
    void print( void );

};

#endif