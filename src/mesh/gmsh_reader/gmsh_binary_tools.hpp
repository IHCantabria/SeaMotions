
#pragma once

// Include general usage libraries
#include <fstream>


template<typename T>
void read_binary( std::ifstream &in, T &value )
{
    in.read( reinterpret_cast<char*>( &value ), sizeof( T ) );
}

template<typename T>
void read_binary_array( std::ifstream &in, T *data, size_t count ) 
{
    in.read( reinterpret_cast<char*>( data ), sizeof( T ) * count );
}