
#ifndef __tools_hpp
#define __tools_hpp

// Include general usage libraries
#include <string>

// Include local modules
#include "config.hpp"


///////////////////////////////////////////////
/************** MACRO DEFINITION *************/
///////////////////////////////////////////////
#define GET_PROGRAM_POINT()(                \
    std::stringstream ss;                   \
    ss << "FILE: " << __FILE__;             \
    ss << " - LINE: " << __LINE__ << "\n";  \
    ss.str();                               \
)                                           \


#define CHECK_FILE_UNIT_STATUS( file_id, file_path ){                           \
    if ( !file_id.good( ) )                                                     \
    {                                                                           \
        std::cerr << "\nERROR - INPUT DATA:" << std::endl;                      \
        std::cerr << " -> INPUT FILE PATH: " << file_path;                      \
        std::cerr << " - does not exits. Check input parameters." << std::endl; \
                                                                                \
        if ( _DEBUG_BUILD )                                                     \
        {                                                                       \
            std::cerr << "FILE: " << __FILE__ << " - ";                         \
            std::cerr << "LINE: " << __LINE__;                                  \
        }                                                                       \
        std::cerr << std::endl;                                                 \
        exit(10);                                                               \
    }                                                                           \
}


#define CHECK_INPUT_FILE_VERSION( current_version, file_version, file_path )    \
    if ( current_version.compare( file_version ) != 0 )                         \
    {                                                                           \
        std::cerr << std::endl;                                                 \
        std::cerr << "ERROR - INPUT FILE VERSION" << std::endl;                 \
        std::cerr << " - Input file: " << file_path << std::endl;               \
        std::cerr << "has an unexpected file version: " << file_version;        \
        std::cerr << ". The current program is expecting: " << current_version; \
        std::cerr << std::endl;                                                 \
        exit(12);                                                               \
    }                                                                           \


///////////////////////////////////////////////
/************ FUNCTION DEFINITION ************/
///////////////////////////////////////////////
                                std::string align_str( std::string input, int width, int align );
template<typename T>    inline  std::string align_num( T number, int width, int precision, int align, int scientific_flag );
                                bool        check_num_cmd_args( int argc, int req_argc );
                                double      get_wall_time( );
                                double      get_cpu_time( );
                                bool        is_empty_line( std::string );
                                void        renew_stream( std::istringstream& iss, std::string line );

#include "tools.txx"

#endif