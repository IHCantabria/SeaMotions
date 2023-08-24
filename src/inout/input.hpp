
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
    std::vector<cusfloat>       angfreqs        ;
    int                         angfreqs_np     = 0;
    std::string                 freqs_unit      = "";
    std::vector<cusfloat>       heads           ;
    int                         heads_np        = 0;
    std::string                 heads_units     = "";
    cusfloat                    grav_acc        = 0.0;
    cusfloat                    water_density   = 0.0;
    cusfloat                    water_depth     = 0.0;

    // Define class constructors and destructors

    // Define class methods
    void    configure( 
                        void
                    );

    void    print(     
                        void 
                );

};

#endif