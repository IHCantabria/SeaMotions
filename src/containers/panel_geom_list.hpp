
#ifndef __panel_geom_list_hpp
#define __panel_geom_list_hpp

// Include local modules
#include "panel_geom.hpp"


struct PanelGeomList
{
    // Define class attributes
    int         panel_np;
    PanelGeom** panels;

    // Define constructors and destructor
    PanelGeomList( ) = default;

    PanelGeomList(
                    int         panel_np_in,
                    PanelGeom** panels_in
                )
    {
        this->panel_np  = panel_np_in;
        this->panels    = panels_in;
    }

    ~PanelGeomList( void )
    {
        for ( int i=0; i<panel_np; i++ )
        {
            delete [] this->panels[i];
        }
        delete [] this->panels;
    }
};

#endif