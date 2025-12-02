
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


                                void            read_body(
                                                                            Input*                          input,
                                                                            std::string                     folder_path,
                                                                            std::string                     target_file,
                                                                            BodyDef*                        body
                                                        );

                                void            read_bodies( 
                                                                            std::string folder_path,
                                                                            Input*      input
                                                            );

template<typename T>    inline  void            read_channel_list( 
                                                                            std::ifstream&                  infile,
                                                                            std::string                     target_file,
                                                                            std::string                     target_signal,
                                                                            int&                            line_count,
                                                                            std::vector<T>&                 list
                                                                    );

                                std::string     read_channel_name( 
                                                                            std::ifstream&                  infile 
                                                                );

template<typename T>    inline  std::string     read_channel_value( 
                                                                            std::ifstream&                  infile, 
                                                                            T&                              channel
                                                                    );

template<typename T>    inline  void            read_channel_matrix( 
                                                                            std::ifstream&                  infile, 
                                                                            int                             rows_np, 
                                                                            int                             cols_np, 
                                                                            T                               channel
                                                                    );

template<typename T>    inline  void            read_list_contraction(
                                                                            std::ifstream&                  infile,
                                                                            int&                            line_count,
                                                                            std::string                     target_file,
                                                                            std::vector<T>&                 container
                                                                    );

template<typename T>    inline  void            read_compact_list(
                                                                            std::string                     item,
                                                                            std::vector<T>&                 container
                                                                    );

                                Input*          read_input_files( 
                                                                            std::string                     folder_path
                                                                );

                                void            read_input_files( 
                                                                            Input*                          input,
                                                                            std::string                     folder_path
                                                                );

template<typename T>    inline  void    read_table_row_header(
                                                                            std::ifstream&                  infile,
                                                                            int&                            line_count,
                                                                            std::vector<std::string>&       row_header,
                                                                            std::vector<std::vector<T>>&    row_values
                                                            );

                                void            skip_header(           
                                                                            std::ifstream&                  infile, 
                                                                            int&                            line_count, 
                                                                            int                             np 
                                                            );


// Include temaplates definition
#include "reader.txx"


#endif