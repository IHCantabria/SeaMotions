
#ifndef __input_hpp
#define __input_hpp

// Include general usage libraries
#include <iostream>
#include <string>

// Import local modules
#include "../config.hpp"
#include "../containers/body_def.hpp"
#include "../mesh/mesh.hpp"


struct Input
{
public:
    // Define class attributes
    BodyDef**                   bodies          = nullptr;
    std::vector<std::string>    bodies_finame   ;
    int                         bodies_np       = 0;
    std::string                 case_fopath     = "";
    std::string                 folder_path     = "";
    std::vector<cusfloat>       angfreqs        ;
    int                         angfreqs_np     = 0;
    std::string                 freqs_unit      = "";
    std::vector<cusfloat>       heads           ;
    int                         heads_np        = 0;
    std::string                 heads_units     = "";
    bool                        is_bodies       = false;
    cusfloat                    grav_acc        = 0.0;
    cusfloat                    water_density   = 0.0;
    cusfloat                    water_depth     = 0.0;

    // Define class constructors and destructors
    ~Input( void );

    // Define class methods
    void    configure( 
                        void
                    );

    void    print(     
                        void 
                );

};

#endif