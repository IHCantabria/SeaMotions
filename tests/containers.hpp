
#ifndef __test_containers_hpp
#define __test_containers_hpp

// Include general usage libraries
#include <fstream>
#include <string>

// Include local modules
#include "../src/config.hpp"

struct DataRef
{
    int num_channels = 0;
    int num_points = 0;
    cusfloat* x = nullptr;
    cusfloat* y = nullptr;
    cusfloat** data = nullptr;

    ~DataRef()
    {
        // Deallocate memory
        delete [] this->x;
        delete [] this->y;

        if (this->num_channels > 0)
        {
            for (int i=0; i<this->num_channels; i++)
            {
                delete [] this->data[i];
            }
            delete [] this->data;
        }
    }

    void read_multiple_channel(std::string file_path)
    {
        // Declare local variables
        std::string s0, s1, s2;

        // Open file unit
        std::ifstream infile(file_path);

        // Read number of points
        infile >> s0 >> s1 >> s2;
        infile >> this->num_points;

        // Read number of channels
        infile >> s0 >> s1 >> s2;
        infile >> this->num_channels;

        // Allocate heap memory for the data
        this->x = new cusfloat [this->num_points];
        this->data = new cusfloat* [this->num_points];
        for (int i=0; i<this->num_points; i++)
        {
            this->data[i] = new cusfloat [num_channels];
        }

        // Read data header from file
        infile >> s0;
        for (int i=0; i<this->num_channels; i++)
        {
            infile >> s1;
        }

        // Read data from file
        for (int i=0; i<this->num_points; i++)
        {
            infile >> this->x[i];
            for (int j=0; j<this->num_channels; j++)
            {
                infile >> this->data[i][j];
            }
        }

        // Close file unit
        infile.close();

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