
#ifndef __test_panel_geometries_hpp
#define __test_panel_geometries_hpp

#include "../../src/config.hpp"
#include "../../src/containers/containers.hpp"

void define_45_inclined_panel(PanelGeom &panel);
void define_field_points_set_1(int &num_points, cusfloat* &data);
void define_square_panel(PanelGeom &panel, cusfloat scale);
void fill_circle_points(cusfloat* data, cusfloat r, cusfloat z, int &count_point);

#endif