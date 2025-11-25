
#pragma once

// Include general usage scientific libraries
#include <vector>
#include <iostream>


struct GmshElement
{
public:
    int                 id      = -1;
    int                 type    = -1;
    std::vector<int>    nodes;

    // ------------------------------------------------------------
    // Default constructor
    // ------------------------------------------------------------
    GmshElement()
        : id(-1),
          type(-1),
          nodes()
    {
    }

    // ------------------------------------------------------------
    // Copy constructor
    // ------------------------------------------------------------
    GmshElement(const GmshElement& other)
        : id(other.id),
          type(other.type),
          nodes(other.nodes)   // deep copy
    {
        std::cout << "Copy Constructor..." << std::endl << std::flush;
    }

    // ------------------------------------------------------------
    // Move constructor
    // ------------------------------------------------------------
    GmshElement(GmshElement&& other) noexcept
        : id(other.id),
          type(other.type),
          nodes(std::move(other.nodes))   // steal vector storage
    {
        std::cout << "Move Constructor..." << std::endl << std::flush;

        // optional: reset moved-from object
        other.id = -1;
        other.type = -1;
    }

    // ------------------------------------------------------------
    // Copy assignment operator
    // ------------------------------------------------------------
    GmshElement& operator=(const GmshElement& other)
    {
        std::cout << "Copy assignement operator..." << std::endl << std::flush;
        if (this != &other)
        {
            id    = other.id;
            type  = other.type;
            nodes = other.nodes;  // deep copy
            for ( std::size_t i=0; i<other.nodes.size( ); i++ )
            {
                std::cout << "I: " << this->nodes[i] << " - " << other.nodes[i] << std::endl << std::flush;
            }
        }
        std::cout << "Copy assignement operator... -> Done!" << std::endl << std::flush;
        return *this;
    }

    // ------------------------------------------------------------
    // Move assignment operator
    // ------------------------------------------------------------
    GmshElement& operator=(GmshElement&& other) noexcept
    {
        std::cout << "Move assignement operator..." << std::endl << std::flush;
        if (this != &other)
        {
            id    = other.id;
            type  = other.type;
            nodes = std::move(other.nodes);

            other.id = -1;
            other.type = -1;
        }
        return *this;
    }

    // ------------------------------------------------------------
    // Destructor (optional, but shown for completeness)
    // ------------------------------------------------------------
    ~GmshElement() = default;
};