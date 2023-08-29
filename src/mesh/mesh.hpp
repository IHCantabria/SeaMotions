
#ifndef __mesh_hpp
#define __mesh_hpp

// Include general usage libraries
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// Include general usage scientific libraries
#include "mkl.h"

// Include local modules
#include "../../src/config.hpp"
#include "../../src/containers/panel_geom.hpp"
#include "../../src/math/math_tools.hpp"
#include "../../src/tools.hpp"


struct Mesh
{
private:
    // Define class attributes
    int valid_elem_type[2]  = { 3, 4 };
    int valid_elem_type_np  = 2;

    // Define class methods
    void    _create_panels(
                                void
                            );
    
    bool    _is_valid_type( 
                                int elem_type 
                           );

    void    _load_poly_mesh( 
                                std::string file_path 
                           );

public:
    // Define class attributes
    int*        elems       = nullptr;
    int         elems_np    = 0;
    int         enrl        = 0;
    int         mnpe        = 0;
    int         nodes_np    = 0;
    PanelGeom** panels      = nullptr;
    cusfloat*   x           = nullptr;
    cusfloat*   y           = nullptr;
    cusfloat*   z           = nullptr;

    // Define class constructor and destructor
    Mesh( std::string file_path );

    ~Mesh( void );

    // Define class methods
    void get_elem_nodes( 
                            int         elem_num, 
                            int&        npe, 
                            cusfloat*   xn, 
                            cusfloat*   yn,
                            cusfloat*   zn
                        );
    
};

#endif