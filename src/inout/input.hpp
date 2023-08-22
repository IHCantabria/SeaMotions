
#ifndef __input_hpp
#define __input_hpp

// Include general usage libraries
#include <iostream>

// Import local modules
#include "../config.hpp"
#include "../mesh/mesh.hpp"


struct BodyDef
{
public:
    // Define class attributes
    Mesh*       mesh            = nullptr;

    // Define class constructor and destructor


};


struct Input
{
public:
    // Define class attributes
    std::string                 case_fopath     = "";
    std::vector<std::string>    bodies_finame   ;
    std::string                 folder_path     = "";
    int                         bodies_np       = 0;
    cusfloat*                   freqs           = nullptr;
    int                         freqs_np        = 0;
    cusfloat*                   heads           = nullptr;
    int                         heads_np        = 0;
    cusfloat                    grav_acc        = 0.0;
    cusfloat                    water_density   = 0.0;
    cusfloat                    water_depth     = 0.0;

    // Define class constructors and destructors

    // Define class methods
    void print( void );

};

#endif