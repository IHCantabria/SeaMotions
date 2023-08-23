
// Include general usage libraries
#include <iostream>
#include <sstream>
#include <string>


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
inline bool is_string( void )
{
    return std::is_same_v<std::remove_cv_t<std::remove_reference_t<T>>, std::string>;
}