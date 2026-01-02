
/*
 * Copyright (c) 2025 Sergio Fernández Ruano / IHCantabria
 *
 * This file is part of SeaMotions Software.
 *
 * SeaMotions is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SeaMotions is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */


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


#define LOG_TASK_SS( loglabel, msg )                        \
    Logger loglabel##_logger( mpi_config->is_root( ) );     \
    {                                                       \
        std::stringstream ss;                               \
        ss << msg;                                          \
        loglabel##_logger.task( ss.str( ) );                \
    }                                                       \


#define LOG_TASK_TIME( loglabel, timer )        \
loglabel##_logger.task_time( timer );           \


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
