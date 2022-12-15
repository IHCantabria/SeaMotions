
#ifndef __test_containers_hpp
#define __test_containers_hpp

// Include general usage libraries
#include <fstream>
#include <string>

// Include local modules
#include "../src/config.hpp"

struct DataRef
{
    int num_points = 0;
    cusfloat* x = nullptr;
    cusfloat* y = nullptr;

    ~DataRef()
    {
        delete [] this->x;
        delete [] this->y;
    }

    void read_single_channel(std::string file_path)
    {
        // Declare local variables
        std::string s0, s1, s2;

        // Open file unit
        std::ifstream infile(file_path);

        // Read number of points
        infile >> s0 >> s1 >> s2;
        infile >> this->num_points;

        // Allocate heap memory for the data
        this->x = new cusfloat [this->num_points];
        this->y = new cusfloat [this->num_points];

        // Read data from file
        infile >> s0 >> s1;
        for (int i=0; i<this->num_points; i++)
        {
            infile >> this->x[i] >> this->y[i];
        }

        // Close file unit
        infile.close();

    }

};

#endif