
// Include local modules
#include "../config.hpp"


template<typename Functor>
cusfloat romberg_quadrature(
                                Functor f, 
                                cusfloat a, 
                                cusfloat b, 
                                double precision
                            )
{
    // Define buffers
    const int max_steps = 50;
    cusfloat* rp = new cusfloat [max_steps];
    cusfloat* rc = new cusfloat [max_steps];

    // Calculate maximum step size
    cusfloat h = b - a;

    // Calculate first trapezoidal point
    rp[0] = (f(a)+f(b))*h/2.0;

    // Loop until the maximum step limit
    cusfloat* aux_ptr = nullptr;
    cusfloat local_trap = 0.0;
    int trap_steps = 0;
    for (int i=1; i<max_steps; i++)
    {
        // Calculate i trapezoidal rule value
        local_trap = 0.0;
        trap_steps = 1 << (i-1);
        h /= 2.0;
        for (int j=1; j<=trap_steps; j++)
        {
            // std::cout << " --> Root: " << a+(2*j-1)*h << " - Value: " << f(a+(2*j-1)*h) << std::endl;
            local_trap += f(a+(2*j-1)*h);
        }
        local_trap *= h;

        // Calculate current buffer values
        rc[0] = local_trap + rp[0]/2.0;
        for (int j=1; j<=i; j++)
        {
            rc[j] = rc[j-1] + 1/(pow(4, j)-1)*(rc[j-1]-rp[j-1]);
        }

        // Check for convergence
        if (std::abs(rc[i]-rp[i-1])<precision)
        {
            return rc[i];
        }

        // Change buffers to have the current as the previous
        aux_ptr = rc;
        rc = rp;
        rp = aux_ptr;

    }

    // Get value to return
    cusfloat last_value = rp[max_steps-1];

    // Delete heap memory
    delete [] rp;
    delete [] rc;

    // Return the last value
    std::cout << "WARNING: Romberg quadrature could not find the integral ";
    std::cout << "value with the requested accuracy." << std::endl;
    return last_value;
}