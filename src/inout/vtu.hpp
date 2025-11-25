
#pragma once

// Include local modules
#include "../config.hpp"


bool    write_vtu_binary_appended(
                                    const std::string       &filename,
                                    const std::size_t       nodes_np,
                                    cusfloat*               nodes_x,
                                    cusfloat*               nodes_y,
                                    cusfloat*               nodes_z,
                                    const std::size_t       elems_np,
                                    const int               enrl,
                                    int*                    elements,
                                    int*                    elements_type
                                );


bool    write_vtu_ascii(
                                    const std::string       &filename,
                                    const std::size_t       nodes_np,
                                    cusfloat*               nodes_x,
                                    cusfloat*               nodes_y,
                                    cusfloat*               nodes_z,
                                    const std::size_t       elems_np,
                                    const int               enrl,
                                    int*                    elements,
                                    int*                    elements_type
                        );