
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
    }

    // ------------------------------------------------------------
    // Move constructor
    // ------------------------------------------------------------
    GmshElement(GmshElement&& other) noexcept
        : id(other.id),
          type(other.type),
          nodes(std::move(other.nodes))   // steal vector storage
    {
        // optional: reset moved-from object
        other.id = -1;
        other.type = -1;
    }

    // ------------------------------------------------------------
    // Copy assignment operator
    // ------------------------------------------------------------
    GmshElement& operator=(const GmshElement& other)
    {
        if (this != &other)
        {
            id    = other.id;
            type  = other.type;
            nodes = other.nodes;  // deep copy
        }
        return *this;
    }

    // ------------------------------------------------------------
    // Move assignment operator
    // ------------------------------------------------------------
    GmshElement& operator=(GmshElement&& other) noexcept
    {
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