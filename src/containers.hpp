
#ifndef __containers_hpp
#define __containers_hpp

#include "config.hpp"

struct PanelNodes
{
    int num_points = 0;
    cusfloat x[4];
    cusfloat xm = 0;
    cusfloat y[4];
    cusfloat ym = 0;
    cusfloat z[4];
    cusfloat zm = 0;

    // Add method to calculate the geometric propertiess
    void calculateGeomProperties(void)
    {
        // Calculate ceter of the panel
        for (int i=0; i<this->num_points; i++)
        {
            this->xm += this->x[i];
            this->ym += this->y[i];
            this->zm += this->z[i];
        }
        this->xm /= this->num_points;
        this->ym /= this->num_points;
        this->zm /= this->num_points;

    }
};

#endif