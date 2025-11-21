
#ifndef __tools_hpp
#define __tools_hpp

// Include general usage libraries
#include <string>
#include <vector>

// Include local modules
#include "config.hpp"


///////////////////////////////////////////////
/************** MACRO DEFINITION *************/
///////////////////////////////////////////////

#define CREATE_ATTRIBUTE(                                                       \
                            fid,                                                \
                            attr_name,                                          \
                            data,                                               \
                            data_type                                           \
                        )                                                       \
{                                                                               \
    H5::DataSpace   attr_space(H5S_SCALAR);                                     \
    H5::Attribute   attr_handle    = fid.createAttribute(                       \
                                                        attr_name,              \
                                                        data_type,              \
                                                        attr_space              \
                                                    );                          \
    attr_handle.write( data_type, &data );                                      \
}                                                                               \

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


#define CREATE_DATASET(                                                                             \
                            fid,                                                                    \
                            dset_name,                                                              \
                            shape_np,                                                               \
                            shape,                                                                  \
                            data_type                                                               \
                        )                                                                           \
{                                                                                                   \
    /* Define DataSpace for all the channel */                                                      \
    H5::DataSpace   dspce           = H5::DataSpace( shape_np, shape );                             \
                                                                                                    \
    /* Create Dataset */                                                                            \
    H5::DataSet     dataset         = fid.createDataSet( dset_name, data_type, dspce );             \
                                                                                                    \
    /* Close DataSet */                                                                             \
    dataset.close();                                                                                \
}                                                                                                   \


#define GET_PROGRAM_POINT()(                                                                        \
    std::stringstream ss;                                                                           \
    ss << "FILE: " << __FILE__;                                                                     \
    ss << " - LINE: " << __LINE__ << "\n";                                                          \
    ss.str();                                                                                       \
)                                                                                                   \


#define INFO( message )                                                                             \
if ( mpi_config->is_root( ) ) std::cout << message;                                                 \


#define MPI_TIC( name )                                                                             \
    MPI_Barrier( MPI_COMM_WORLD );                                                                  \
    double name##_tic = MPI_Wtime( );                                                               \


#define MPI_TOC( name )                                                                             \
    MPI_Barrier( MPI_COMM_WORLD );                                                                  \
    double name##_toc       = MPI_Wtime( );                                                         \
    double name##_elapsed   = name##_toc - name##_tic;                                              \


#define SAVE_DATASET_CHUNK(                                                                         \
                                fid,                                                                \
                                dset_name,                                                          \
                                shape_np,                                                           \
                                shape,                                                              \
                                chunk_shape,                                                        \
                                offset_shape,                                                       \
                                data,                                                               \
                                data_type                                                           \
                            )                                                                       \
{                                                                                                   \
    /* Define DataSpace for all the channel */                                                      \
    H5::DataSpace   dspce           = H5::DataSpace( shape_np, shape );                             \
                                                                                                    \
    /* Open dataset */                                                                              \
    H5::DataSet     dataset         = fid.openDataSet( dset_name );                                 \
                                                                                                    \
    /* Get DataSet properties to configure the chunk */                                             \
    H5::DSetCreatPropList prop_list = dataset.getCreatePlist();                                     \
    prop_list.setChunk( shape_np, chunk_shape );                                                    \
                                                                                                    \
    /* Create DataSet for the chunk */                                                              \
    H5::DataSpace   dspce_chunk( shape_np, chunk_shape );                                           \
    dspce.selectHyperslab( H5S_SELECT_SET, chunk_shape, offset_shape );                             \
                                                                                                    \
    /* Save data */                                                                                 \
    dataset.write( data, data_type, dspce_chunk, dspce );                                           \
                                                                                                    \
    /* Close DataSet */                                                                             \
    dataset.close();                                                                                \
}                                                                                                   \


///////////////////////////////////////////////
/************ FUNCTION DEFINITION ************/
///////////////////////////////////////////////
                                std::string align_str( 
                                                                            std::string input, 
                                                                            int width, 
                                                                            int align
                                                        );

template<typename T>    inline  std::string align_num( 
                                                                            T number, 
                                                                            int width, 
                                                                            int precision, 
                                                                            int align, 
                                                                            int scientific_flag 
                                                        );

                                bool        check_num_cmd_args( 
                                                                            int argc, 
                                                                            int req_argc 
                                                                );

template<typename T>    inline  void        convert_number( 
                                                                            std::string str, 
                                                                            T& val 
                                                        );
                                                        
template<>              inline  void        convert_number<std::string>( 
                                                                            std::string str,
                                                                            std::string& val 
                                                                        );
                                                                
template<>              inline  void        convert_number<int>( 
                                                                            std::string str, 
                                                                            int& val 
                                                                );

template<>              inline  void        convert_number<cusfloat>( 
                                                                            std::string str, 
                                                                            cusfloat& val 
                                                                    );

                                std::string get_fipath_extension(           
                                                                            std::string fipath 
                                                                );

                                double      get_wall_time( );

                                double      get_cpu_time( );

                                bool        is_empty_line( 
                                                                            std::string
                                                        );

template<typename T>    inline  bool        is_string( 
                                                                            void 
                                                        );

                                void        renew_stream( 
                                                                            std::istringstream& iss, 
                                                                            std::string line 
                                                        );

template<typename T>    inline  void        split_string( 
                                                                            std::string     str,
                                                                            std::vector<T>& vec,
                                                                            char            sep
                                                        );

                                void        squeeze_string( 
                                                                            std::string& str 
                                                            );

                                bool        str_to_bool( 
                                                                            std::string v 
                                                        );

                                void        str_to_lower( 
                                                                            std::string* str 
                                                        );

template<typename T>            std::string vec_to_str(
                                                                            const std::vector<T> &vec
                                                        );

#include "tools.txx"

#endif