
// Include general usage libraries
#include <iostream>
#include <sstream>
#include <vector>

// Include local modules
#include "../tools.hpp"


template<typename T>
inline  void    _read_channel_list( 
                                        std::ifstream&  infile,
                                        std::string     target_file,
                                        std::string     target_signal,
                                        int&            line_count,
                                        std::vector<T>& list
                                    )
{
    // Loop until the next header is found
    T           aux_var;
    int         _count = 0;
    int         _first_line = line_count;
    std::string line;
    int         max_iter = 1000;
    std::string read_signal;
    while ( true )
    {
        // Read the current file position
        auto read_pos = infile.tellg( );

        // Check if reading a heder line
        std::getline( infile, line );
        int pos = line.find( "/***" );
        infile.seekg( read_pos );
        if ( !( pos < 0 ) )
        {
            break;
        }

        // Get signal value
        read_signal     = _read_channel_value( infile, aux_var );
        CHECK_SIGNAL_NAME( read_signal, target_signal, target_file, line_count );

        // Storage signal value
        list.push_back( aux_var );

        // Check maximum number of iterations
        _count++;
        if ( _count > max_iter )
        {
            std::cerr << std::endl;
            std::cerr << "Not posible to find the end of the list that begins at line: " << _first_line;
            std::cerr << " of the file: " << target_file << std::endl;
            std::runtime_error( "" );
        }

    }

}


template<typename T>
inline std::string  _read_channel_value( 
                                            std::ifstream& infile, 
                                            T& signal
                                        )
{
    // Generate local auxiliar variables
    std::string aux_str, line, signal_name;

    // Read line from file unit
    std::getline(infile, line);

    // Get signal value and name from the line
    // read
    std::istringstream iss(line);
    if ( typeid( decltype(signal) ) == typeid( bool ) )
    {
        std::string signal_str;
        iss >> signal_str >> signal_name >> aux_str;
        signal = str_to_bool( signal_str );
    }
    else
    {
        iss >> signal >> signal_name >> aux_str;
    }
    
    return signal_name;
}


template<typename T>
inline  void    _read_channel_matrix( 
                                        std::ifstream&  infile, 
                                        int             rows_np, 
                                        int             cols_np, 
                                        T               channel
                                    )
{
    // Declare local auxiliar variables
    std::string aux_str, channel_name, line;

    // Read channel name
    for ( int i=0; i<rows_np; i++ )
    {
        std::getline( infile, line );
        std::istringstream iss(line);
        for ( int j=0; j<cols_np; j++ )
        {
            iss >> channel[i*cols_np+j];
        } 
    }

}