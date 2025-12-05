
/*
 * Copyright (c) 2025 Sergio Fernández Ruano / IHCantabria
 *
 * This file is part of SeaMotions Software.
 *
 * SeaMotions is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SeaMotions is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

// Include general usage libraries
#include <fstream>
#include <string>
#include <vector>

// Include local modules
#include "gmsh_element.hpp"
#include "gmsh_node.hpp"


int     gmsh_element_node_count( 
                                    int type 
                                );

bool    read_nodes( 
                                    std::ifstream &in, 
                                    std::vector<GmshNode> &nodes 
                    );

bool    read_nodes_binary( 
                                    std::ifstream &in, 
                                    std::vector<GmshNode> &nodes 
                        );

bool    read_elements( 
                                    std::ifstream &in, 
                                    std::vector<GmshElement> &elements 
                    );

bool    read_elements_binary( 
                                    std::ifstream &in, 
                                    std::vector<GmshElement> &elements 
                            );

bool    read_gmsh_41(
                                    const std::string &filename,
                                    std::vector<GmshNode> &nodes,
                                    std::vector<GmshElement> &elements
                    );

int     gmsh_to_vtk( 
                                    int gmsh 
                    );
