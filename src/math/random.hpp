
#ifndef __random_hpp
#define __random_hpp


void        check_vsl_error( int num );
int64_t     generate_random_seed( void );
void        randn( double mean, double std, int seed, int np, double* data );
void        randu( double a, double b, int seed, int np, double* data );


#endif