
// Include general usage libraries
#include <algorithm>
#include <cctype>
#include <iostream>
#include <sstream>
#include <string>

// Include local modules
#include "tools.hpp"


std::string align_str(std::string input, int width, int align)
{
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
        align_data = s + input;
    }
    else if (align == 1)
    {
        int len_diff = width-input.length();
        int right_space = len_diff/2;
        int left_space = right_space + len_diff%2;

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


bool check_num_cmd_args(int argc, int req_argc)
{
    if (argc < (req_argc+1))
    {
        std::cerr << "Not enough input command line input arguments. ";
        std::cerr << "Received: " << argc-1 << " - Expected :" << req_argc << std::endl;
        return false;
    }
    else if (argc > (req_argc+1))
    {
        std::cerr << "More input arguments than expected." << std::endl;
        std::cerr << "Received: " << argc-1 << " - Expected :" << req_argc << std::endl;
        return false;
    }

    return true;
}


bool is_empty_line( std::string line )
{
    bool is_empty = true;
    for (char& c : line )
    {
        if (
                !(
                    c == ' '
                    ||
                    c == '\t'
                    ||
                    c == '\n'
                )
            )
        {
            is_empty = false;
        }
    }

    return is_empty;
}


void renew_stream( 
                    std::istringstream& iss,
                    std::string         line 
                )
{
    iss.str( std::string() );
    iss.clear( );
    iss.str( line );
}


void squeeze_string( std::string& str )
{
    str.erase(std::remove_if(str.begin( ), str.end( ), static_cast<int(&)(int)>( std::isspace )), str.end( ));
}


bool str_to_bool( std::string v )
{
    bool v_bool = false;
    if ( v.compare( std::string("true") ) == 0 )
    {
        v_bool = true;
    }
    else if ( v.compare( std::string("false") ) == 0 )
    {
        v_bool = false;
    }
    else
    {
        std::cerr << "Not possible to convert: " << v << " to boolean value." << std::endl;
        std::cerr << "Accepted values are: true | false" << std::endl;
        exit(10);
    }

    return v_bool;
}


void str_to_lower( std::string* str )
{
    std::transform(
                        str->begin(), 
                        str->end(), 
                        str->begin(),
                        [](unsigned char c){ return std::tolower(c); }
                    );
}

//  Windows
#ifdef _WIN32
#include <Windows.h>
double get_wall_time()
{
    LARGE_INTEGER time,freq;
    if (!QueryPerformanceFrequency(&freq)){
        //  Handle error
        return 0;
    }
    if (!QueryPerformanceCounter(&time)){
        //  Handle error
        return 0;
    }
    return (double)time.QuadPart / freq.QuadPart;
}
double get_cpu_time()
{
    FILETIME a,b,c,d;
    if (GetProcessTimes(GetCurrentProcess(),&a,&b,&c,&d) != 0){
        //  Returns total user time.
        //  Can be tweaked to include kernel times as well.
        return
            (double)(d.dwLowDateTime |
            ((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
    }else{
        //  Handle error
        return 0;
    }
}

//  Posix/Linux
#else
#include <time.h>
#include <sys/time.h>
double get_wall_time()
{
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time()
{
    return (double)clock() / CLOCKS_PER_SEC;
}
#endif
