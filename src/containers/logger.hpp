
#pragma once

// Include general usage libraries
#include <iostream>
#include <iomanip>
#include <mutex>
#include "mpi.h"
#include <sstream>
#include <string>

// Include local module
#include "mpi_timer.hpp"



/*******************************************************************************/
/***************************** Module Macros ***********************************/
/*******************************************************************************/
#define LOG_TASK( msg )                         \
    Logger logger( mpi_config->is_root( ) );    \
    logger.task( msg );                         \


#define LOG_TASK_SS( msg )                      \
    std::stringstream ss;                       \
    ss << msg;                                  \
    Logger logger( mpi_config->is_root( ) );    \
    logger.task( ss.str( ) );                   \


#define LOG_TASK_TIME( timer )                  \
logger.task_time( timer );                      \


/*******************************************************************************/
/************************** ANSI Terminal Colors *******************************/
/*******************************************************************************/
namespace TermColor {
                        static constexpr const char* RESET   = "\033[0m";
                        static constexpr const char* RED     = "\033[31m";
                        static constexpr const char* YELLOW  = "\033[33m";
                        static constexpr const char* GREEN   = "\033[32m";
                        static constexpr const char* CYAN    = "\033[36m";
                        static constexpr const char* WHITE   = "\033[37m";
                    }


/*******************************************************************************/
/************************ Define Enum Class Level ******************************/
/*******************************************************************************/
enum class Level { INFO, WARNING, ERROR, TASK };


/*******************************************************************************/
/****************************** Logger Class ***********************************/
/*******************************************************************************/
class Logger
{
private:
    /* Declare private attributes */
    bool                _is_root     = true;
    std::string         _last_msg;
    int                 _time_column = 40;
    std::mutex          _mutex;

    /* Declare private class methods */
    const char*     _color( Level level ) const;

    std::string     _prefix( Level lvl ) const;

public:
    /* Declare class constructors */
    // No arguments contructor. It checks the root process 
    // internally. This allows to use the logger out of 
    // block codes where MpiConfig is present.
    Logger( );

    // Constructor that allows to change root and time column
    Logger( bool is_root );

    // Constructor that allows to change root and time column
    Logger( bool is_root, int time_column );

    /* Declare class methods */
    void log( Level level, const std::string& msg);
    
    void error( const std::string &msg );

    void info( const std::string &msg );

    void set_time_column( int col );
    
    void task( const std::string &msg );
    
    void task_time( const std::string& msg, MpiTimer& timer );

    void task_time( MpiTimer& timer);

    void warn( const std::string &msg );

};
