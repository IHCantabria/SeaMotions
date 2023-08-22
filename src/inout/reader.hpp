
#ifndef __reader_hpp
#define __reader_hpp

// Include general usage libraries
#include <fstream>
#include <iostream>
#include <string>

// Include local modules
#include "input.hpp"


#define CHECK_SIGNAL_NAME( signal_name, signal_name_ref, file_name, line_num ) do {         \
    if ( signal_name_ref.compare( signal_name ) != 0 )                                      \
    {                                                                                       \
        std::cerr << std::endl;                                                             \
        std::cerr << "ERROR - Input Reading" << std::endl;                                  \
        std::cerr << "Channel: " << signal_name_ref << " not found on ";                    \
        std::cerr << "line: " << line_num+1 << " of the file: " << file_name << std::endl;  \
        std::cerr << "Instead this is the channel name found: " << signal_name << std::endl;\
        if ( _DEBUG_BUILD )                                                                 \
        {                                                                                   \
            std::cerr << "Program execution error at line: " << __LINE__;                   \
            std::cerr << " of the source code file: " << __FILE__ << std::endl;             \
        }                                                                                   \
        exit(10);                                                                           \
    }                                                                                       \
    line_num++;                                                                             \
} while (0)


                                Input*      read_input_files( 
                                                                    std::string         folder_path
                                                            );

template<typename T>    inline  void        _read_channel_list( 
                                                                    std::ifstream&      infile,
                                                                    std::string         target_file,
                                                                    std::string         target_signal,
                                                                    int&                line_count,
                                                                    std::vector<T>&     list
                                                                );

                                std::string      _read_channel_name( 
                                                                    std::ifstream&      infile 
                                                                );

template<typename T>    inline  std::string      _read_channel_value( 
                                                                    std::ifstream&      infile, 
                                                                    T&                  channel
                                                                );

template<typename T>    inline  void        _read_channel_matrix( 
                                                                    std::ifstream&      infile, 
                                                                    int                 rows_np, 
                                                                    int                 cols_np, 
                                                                    T                   channel
                                                                );

                                void        _skip_header(           
                                                                    std::ifstream&      infile, 
                                                                    int&                line_count, 
                                                                    int                 np 
                                                        );


// Include temaplates definition
#include "reader.txx"


#endif