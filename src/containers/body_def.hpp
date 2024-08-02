
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
    cusfloat        cog[3]          = { 0.0, 0.0, 0.0 };
    cusfloat        inertia[6]      = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    bool            interia_by_rad  = false;
    bool            is_fix          = false;
    bool            is_mesh         = false;
    bool            is_mesh_fs_qtf  = false;
    int             lid_type        = 0;
    cusfloat        mass            = 0.0;
    Mesh*           mesh            = nullptr;
    Mesh*           mesh_fs_qtf     = nullptr;
    std::string     mesh_finame     = "";
    std::string     mesh_body_name  = "";
    cusfloat        rad_inertia[3]  = { 0.0, 0.0, 0.0 };

    // Define class constructor and destructor
    ~BodyDef( void );

    // Define class methods
    void print( void );

};

#endif