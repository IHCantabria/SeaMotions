
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

// Include general usage libraries
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>

// Include local modules
#include "gmsh_binary_tools.hpp"
#include "read_gmsh.hpp"


/***************************************************************/
/********* Optional: number of nodes for each gmsh type ********/
/***************************************************************/
int gmsh_element_node_count( int type )
{
    switch( type )
    {
        case 1: return 2;    // 2-node line
        case 2: return 3;    // 3-node triangle
        case 3: return 4;    // 4-node quad
        case 4: return 4;    // 4-node tetra
        case 5: return 8;    // 8-node hex
        case 6: return 6;    // 6-node wedge
        case 7: return 5;    // 5-node pyramid
        case 8: return 3;    // 3-node second-order line
        case 9: return 6;    // 6-node second-order triangle
        case 10: return 9;   // 9-node second-order quad
        case 11: return 10;  // 10-node second-order tetra
        case 12: return 27;  // 27-node hex
        case 13: return 18;  // 18-node wedge
        case 14: return 14;  // 14-node pyramid
        case 15: return 1;   // point element
        default:
            return -1;       // means: read from file
    }
}


/***************************************************************/
/************** Read $Nodes block (GMsh 4.1) *******************/
/***************************************************************/
bool read_nodes( std::ifstream &in, std::vector<GmshNode> &nodes )
{
    int num_entity_blocks, num_nodes, min_tag, max_tag;
    in >> num_entity_blocks >> num_nodes >> min_tag >> max_tag;

    nodes.reserve(num_nodes);

    for (int b = 0; b < num_entity_blocks; b++) 
    {
        int entity_dim, entity_tag, parametric, count;
        in >> entity_dim >> entity_tag >> parametric >> count;

        std::vector<int> tags( count );
        for (int i = 0; i < count; i++)
            in >> tags[i];

        for (int i = 0; i < count; i++) 
        {
            GmshNode n;
            n.id = tags[i];
            in >> n.x >> n.y >> n.z;
            nodes.push_back( n );
        }
    }

    std::string end;
    in >> end; // $EndNodes
    return true;
}


bool read_nodes_binary( std::ifstream &in, std::vector<GmshNode> &nodes )
{
    int num_entity_blocks;
    long long num_nodes, min_tag, max_tag;

    read_binary( in, num_entity_blocks );
    read_binary( in, num_nodes );
    read_binary( in, min_tag );
    read_binary( in, max_tag );

    nodes.reserve( num_nodes );

    for (int b = 0; b < num_entity_blocks; b++)
    {
        int entity_dim, entity_tag, parametric;
        long long count;

        read_binary( in, entity_dim );
        read_binary( in, entity_tag );
        read_binary( in, parametric );
        read_binary( in, count );

        std::vector<long long> tags( count );
        read_binary_array( in, tags.data( ), count );

        for (long long i = 0; i < count; i++) {
            GmshNode n;
            n.id = tags[i];
            read_binary( in, n.x );
            read_binary( in, n.y );
            read_binary( in, n.z );
            nodes.push_back( n );

            if ( parametric )
            {
                double u, v, w;
                read_binary( in, u );
                read_binary( in, v );
                read_binary( in, w );
            }
        }
    }

    // Skip $EndNodes
    return true;
}

/***************************************************************/
/*************** Read $Elements block (GMsh 4.1) ***************/
/***************************************************************/
bool read_elements( std::ifstream &in, std::vector<GmshElement> &elements )
{
    int num_entity_blocks, num_elements, min_tag, max_tag;
    in >> num_entity_blocks >> num_elements >> min_tag >> max_tag;

    elements.reserve( num_elements );

    for (int b = 0; b < num_entity_blocks; b++) 
    {
        int entity_dim, entity_tag, element_type, count;
        in >> entity_dim >> entity_tag >> element_type >> count;

        int fixed_count = gmsh_element_node_count( element_type );

        for (int i = 0; i < count; i++) 
        {
            GmshElement e;
            e.type = element_type;

            // read tag
            in >> e.id;

            if (fixed_count > 0) 
            {
                e.nodes.resize( fixed_count );
                for (int k = 0; k < fixed_count; k++)
                    in >> e.nodes[k];
            } 
            else 
            {
                // Unknown type: read until end of line
                std::string rest;
                std::getline(in, rest);
                std::stringstream ss(rest);
                int nid;
                while (ss >> nid)
                    e.nodes.push_back(nid);
            }

            elements.push_back(e);
        }
    }

    std::string end;
    in >> end; // $EndElements
    return true;
}


bool read_elements_binary( std::ifstream &in, std::vector<GmshElement> &elements )
{
    int num_entity_blocks;
    long long numElements, minTag, maxTag;

    read_binary( in, num_entity_blocks );
    read_binary( in, numElements );
    read_binary( in, minTag );
    read_binary( in, maxTag );

    elements.reserve(numElements);

    for ( int b = 0; b < num_entity_blocks; b++ )
    {
        int entity_dim, entity_tag, type;
        long long count;

        read_binary( in, entity_dim );
        read_binary( in, entity_tag );
        read_binary( in, type );
        read_binary( in, count );

        int nPerElem = gmsh_element_node_count(type);

        for (long long i = 0; i < count; i++) 
        {
            GmshElement e;
            e.type = type;

            long long id;
            read_binary(in, id);
            e.id = (int)id;

            e.nodes.resize(nPerElem);
            for (int k = 0; k < nPerElem; k++) 
            {
                long long nid;
                read_binary(in, nid);
                e.nodes[k] = (int)nid;
            }

            elements.push_back(e);
        }
    }

    return true;
}


/***************************************************************/
/***************** Complete GMsh 4.1 Reader ********************/
/***************************************************************/
bool read_gmsh_41(
                    const std::string &filename,
                    std::vector<GmshNode> &nodes,
                    std::vector<GmshElement> &elements
                )
{
    std::ifstream in( filename );
    if ( !in.is_open( ) ) 
    {
        std::cerr << "Cannot open file: " << filename << "\n";
        return false;
    }

    std::string token;
    bool binary = false;
    while ( in >> token )
    {
        if ( token == "$MeshFormat" ) 
        {
            double version; int is_binary, datasize;
            in >> version >> is_binary >> datasize;

            if ( is_binary == 1 ) 
            {
                // skip newline and one int32 check
                int one;
                in.read( (char*)&one, sizeof(int) );
                binary = true;
            }
        }

        if ( token == "$Nodes" ) 
        {
            if ( binary ) read_nodes_binary( in, nodes );
            else read_nodes(in, nodes);
        }
        else if ( token == "$Elements" ) 
        {
            if ( binary ) read_elements_binary( in, elements );
            else read_elements( in, elements );
        }
    }

    return true;
}


/***************************************************************/
/************ Translate GMsh type numbering to VTK *************/
/***************************************************************/
int gmsh_to_vtk( int gmsh )
{
    switch ( gmsh ) 
    {
        case 15: return 1;   // point
        case 1: return 3;    // line
        case 2: return 5;    // triangle
        case 3: return 9;    // quad
        case 4: return 10;   // tetra
        case 5: return 12;   // hex
        case 6: return 13;   // wedge
        case 7: return 14;   // pyramid
        default: return -1;
    }
}
