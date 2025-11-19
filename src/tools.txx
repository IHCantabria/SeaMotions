
// Include general usage libraries
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>


template<class T>
constexpr bool always_false = false;


template<typename T>
inline std::string align_num(T number, int width, int precision, int align, int scientific_flag)
{
    // Check if negative number
    int flag_neg = 0;
    if (number < 0)
    {
        flag_neg = 1;
    }

    // Convert number to string
    std::stringstream ss;
    if (scientific_flag == 1)
    {
        ss << std::fixed << std::setprecision(precision) << std::scientific << number;
    }
    else
    {
        ss << std::fixed << std::setprecision(precision) << number;
    }
    std::string input = ss.str();

    // Check if the length of the column is enought to contain the requiered number
    if ((width-input.length())<2)
    {
        std::string err_message("Column width is shorter than necessary to print the input data.");
        std::cerr << err_message << std::endl;
        throw std::runtime_error(err_message);
    }

    // Compose aligned string
    std::string align_data;
    if (align == 0)
    {
        int len_diff = width-input.length();
        std::string s(len_diff, ' ');
        align_data = s+input;
    }
    else if (align == 1)
    {
        int len_diff = width-input.length();
        int right_space = len_diff/2;
        int left_space = right_space + len_diff%2;

        if (flag_neg == 1)
        {
            right_space++;
            left_space--;
        }

        std::string sl(left_space, ' ');
        std::string sr(right_space, ' ');
        align_data = sl + input + sr;
    }
    else
    {
        std::string err_message("Alignment not valid. Accepted values: 0 - align right | 1 - align center");
        std::cerr << err_message << std::endl;
        throw std::runtime_error(err_message);
    }

    return align_data;
}


template<typename T>
inline  void    convert_number( std::string , T& )
{
    static_assert( always_false<T>, "Not valid type!" );
}


template<>
inline  void    convert_number<std::string>( std::string str, std::string& val )
{
    val = str;
}


template<>
inline  void    convert_number<int>( std::string str, int& val )
{
    val = std::atoi( str.c_str( ) );
}


template<>
inline  void    convert_number<cusfloat>( std::string str, cusfloat& val )
{
    val = std::atof( str.c_str( ) );
}


template<typename T>
inline  bool    is_string( void )
{
    return std::is_same_v<std::remove_cv_t<std::remove_reference_t<T>>, std::string>;
}


template<typename T>
inline  void    split_string( 
                                std::string     str,
                                std::vector<T>& vec,
                                char            sep
                            )
{
    // Declare position indexes
    int pos0 = 0;
    int pos1 = 0;

    // Check if the separators are at the begining and the end of the
    // line are present
    if ( str[0] == sep )
    {
        str.erase( 0, 1 );
    }

    if ( str[str.length( )-1] != sep )
    {
        str = str + sep;
    }

    // Loop over string to be the substring
    std::string substr( "" );
    T           val;
    while ( true )
    {
        // Find position of the next separator character
        pos1 = str.find( sep );

        // Check if there is no more items delimited
        // by the specified separator
        if ( pos1 < 0 )
        {
            break;
        }

        // Read substring
        convert_number<T>( 
                            str.substr( pos0, pos1 - pos0 ), 
                            val 
                        );
        vec.push_back( val );

        // Renew str
        str = str.substr( pos1+1, str.length( ) - (pos1+1) );

    }
}


template<typename T>
std::string vec_to_str( const std::vector<T>& vec )
{
    std::stringstream ss;
    ss << "( " << vec[0];
    for ( size_t i=1; i<vec.size( ); i++ )
    {
        ss << ", " << vec[i];
    }
    ss << " )";

    return ss.str( );
}