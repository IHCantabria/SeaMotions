
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
