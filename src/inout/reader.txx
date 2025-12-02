
// Include general usage libraries
#include <iostream>
#include <sstream>
#include <type_traits>
#include <vector>

// Include local modules
#include "../tools.hpp"


template<typename T>
inline  void    read_channel_list( 
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
        read_signal     = read_channel_value( infile, aux_var );
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
inline std::string  read_channel_value( 
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
inline  void    read_channel_matrix( 
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


template<typename T>
inline  void    read_list_contraction(
                                        std::ifstream&  infile,
                                        int&            line_count,
                                        std::string     target_file,
                                        std::vector<T>& container
                                    )
{
    // Check input data type
    static_assert( 
                    std::disjunction<
                                        std::is_same<cusfloat, T>,
                                        std::is_same<int, T>
                                    >( ),
                    "Not valid type. Function only accepts: int | cusfloat"
                );
    
    // Loop over lines until the next section header is found
    T                           aux_val         ;
    std::vector<std::string>    aux_vec         ;
    std::vector<cusfloat>       cmp_list        ;
    int                         _count          = 0;
    int                         _first_line     = line_count;
    std::string                 line            ;
    int                         max_count       = 1e4;
    int                         pos_ddot        = 0;
    while ( true )
    {
        // Read the current file position
        auto read_pos = infile.tellg( );

        // Check if there is more lines to read
        if (!( std::getline( infile, line ) ))
        {
            break;
        }

        // Check if reading a header line
        int pos_h = line.find( "/***" );
        int pos_c = line.find( "-->" );
        if ( !( pos_h < 0 ) || !( pos_c < 0 ) )
        {
            infile.seekg( read_pos );
            break;
        }

        // Process line
        aux_vec.clear( );
        squeeze_string( line );
        split_string( 
                        line,
                        aux_vec,
                        ','
                    );

        // Loop over line items
        for ( auto item: aux_vec )
        {
            // Check if the item is a compact list
            pos_ddot = item.find( ":" );
            if ( !( pos_ddot < 0 ) )
            {
                read_compact_list(
                                        item,
                                        container
                                    );
            }
            else
            {
                // Convert from string to the target value
                convert_number( 
                                    item, 
                                    aux_val 
                                );
                
                // Storage in the target container
                container.push_back( aux_val );
            }

        }

        // Check maximum number of iterations
        _count++;
        line_count++;
        if ( _count > max_count )
        {
            std::cerr << std::endl;
            std::cerr << "Not posible to find the end of the list that begins at line: " << _first_line;
            std::cerr << " of the file: " << target_file << std::endl;
            throw std::runtime_error( "" );
        }

    }
}


template<typename T>
inline  void    read_compact_list(
                                    std::string     item,
                                    std::vector<T>& container
                                )
{
    // Check input types
    static_assert(
                    std::disjunction<
                                        std::is_same<T, cusfloat>,
                                        std::is_same<T, int>
                                    >( ),
                    "Input type not allowed!"
                );
    
    // Declare local variables
    int             bound_set_end   = 0;
    int             bound_set_init  = 0;
    std::vector<T>  cmp_list;
    T               cmp_list_val    = 0;
    int             count_cmp_list  = 0;
    int             last_char       = 0;

    // Check boundary type of the set
    if ( item[0] == '(' )
    {
        bound_set_init = 0;
    }
    else if ( item[0] == '[' )
    {
        bound_set_init = 1;
    }
    else
    {
        std::cerr << std::endl;
        std::cerr << "ERROR: " << std::endl;
        std::cerr << "Not possible to find the correct set init boundary symbol";
        std::cerr << std::endl;
        throw std::runtime_error( "" );
    }
    item.erase( 0, 1 );

    last_char = item.length( ) - 1;
    if ( item[last_char] == ')' ) 
    {
        bound_set_end = 0;
    }
    else if ( item[last_char] == ']' )
    {
        bound_set_end = 1;
    }
    else
    {
        std::cerr << std::endl;
        std::cerr << "ERROR: " << std::endl;
        std::cerr << "Not possible to find the correct set end boundary symbol";
        std::cerr << std::endl;
        throw std::runtime_error( "" );
    }
    item.erase( last_char, 1 );

    // Split line using double dot character
    split_string(
                    item,
                    cmp_list,
                    ':'
                );
    
    // Use compact form of the cmp_list to define 
    // a list of floating values attached to the container
    if ( bound_set_init == 0 )
    {
        count_cmp_list = 1;
    }
    else
    {
        count_cmp_list = 0;
    }
    
    while ( true )
    {
        // Generate new list value
        cmp_list_val = cmp_list[0] + count_cmp_list*cmp_list[1];

        // Check end bounds
        if ( !( cmp_list_val < cmp_list[2] ) )
        {
            break;
        }

        // Storage new value
        container.push_back( cmp_list_val );

        // Advance counter
        count_cmp_list++;

    }

    if ( bound_set_end == 1 )
    {
        container.push_back( cmp_list[2] );
    }
}


template<typename T>
inline  void    read_table_row_header(
                                        std::ifstream&                  infile,
                                        int&                            line_count,
                                        std::vector<std::string>&       row_header,
                                        std::vector<std::vector<T>>&    row_values
                                    )
{
    // Check input types
    static_assert(
                    std::disjunction<
                                        std::is_same<T, cusfloat>,
                                        std::is_same<T, int>
                                    >( ),
                    "Input type not allowed!"
                );
    
    // Loop over lines until the next section header is found
    std::string                 astr            ;
    T                           a0              ;
    T                           a1              ;
    T                           a2              ;
    T                           a3              ;
    std::string                 line            ;

    while ( true )
    {
        // Read the current file position
        auto read_pos = infile.tellg( );

        // Check if there is more lines to read
        if (!( std::getline( infile, line ) ))
        {
            break;
        }

        // Check if reading a header line
        int pos = line.find( "/***" );
        if ( !( pos < 0 ) )
        {
            infile.seekg( read_pos );
            break;
        }

        // Process line
        std::istringstream iss( line );
        iss >> astr >> a0 >> a1 >> a2 >> a3;

        // Storage values
        row_header.push_back( astr );

        std::vector<T> rw;
        rw.push_back( a0 );
        rw.push_back( a1 );
        rw.push_back( a2 );
        rw.push_back( a3 );

        row_values.push_back( rw );

        // Increment line count
        line_count++;

    }

}