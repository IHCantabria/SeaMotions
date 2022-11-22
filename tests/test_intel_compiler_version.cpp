
// Include general usage libraries
#include <iostream>


int main(void)
{
    std::cout << "This program checks the current version in use of the Intel compiler" << std::endl;
    #if defined(__INTEL_LLVM_COMPILER)
        std::cout << "__INTEL_COMPILER: " << __INTEL_LLVM_COMPILER << std::endl;
        std::cout << "__VERSION: " << __VERSION__ << std::endl;
    #endif 

    #if defined(__INTEL_COMPILER)
        std::cout << "__INTEL_COMPILER: " << __INTEL_COMPILER << std::endl;
    #endif

    return 0;
}