
#ifndef __tools_hpp
#define __tools_hpp


#include <string>

std::string align_str(std::string input, int width, int align);
template<typename T> inline std::string align_num(T number, int width, int precision, int align, int scientific_flag);
bool check_num_cmd_args(int argc, int req_argc);
double get_wall_time();
double get_cpu_time();

#include "tools.txx"

#endif