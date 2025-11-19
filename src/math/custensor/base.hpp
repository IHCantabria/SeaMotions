
#pragma once

// Include local modules
#include "../math_tools.hpp"


/********************************************/
/************** MODULE MACROS ***************/
/********************************************/
#define CUSTENSOR_SHAPE_MISMATCH( sl, sr )                              \
{                                                                       \
    std::stringstream ss;                                               \
    ss << "CusTensor Error:\n";                                         \
    ss << "  Shape mismatch in CusTensor assignment: ";                 \
    ss << vec_to_str( sl ) << " -> " << vec_to_str( sr ) << std::endl;  \
    std::cerr << ss.str( ) << std::endl;                                \
    throw std::runtime_error( ss.str( ) );                              \
}                                                                       