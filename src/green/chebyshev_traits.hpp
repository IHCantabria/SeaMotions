
#ifndef __chebyshev_traits_hpp
#define __chebyshev_traits_hpp

// Include local modules
#include "chebyshev_traits_macros.hpp"
#include "./inf_depth_coeffs/R44.hpp"
#include "./inf_depth_coeffs/R44_dX.hpp"


// Generic trait declaration
template<typename NS>
struct ChebyshevTraits;

// Instantiate traits through the corresponding macro
CHEBYSHEV_2D_TRAITS( R44C );
CHEBYSHEV_2D_TRAITS( R44_dXC );

// template<>                                                                                 
// struct ChebyshevTraits<R44C>                                                                 
// {                                                                                          
//     static constexpr            int             max_ref_level       = R44C::max_ref_level;   
//     static constexpr            int             intervals_np        = R44C::intervals_np;    
//     static constexpr            int             max_cheby_order     = R44C::max_cheby_order; 
//     static constexpr const      std::size_t*    blocks_start        = R44C::blocks_start;    
//     static constexpr const      std::size_t*    blocks_coeffs_np    = R44C::blocks_coeffs_np;
//     static constexpr            bool            fcn_log_scale       = R44C::fcn_log_scale;   
//     static constexpr            bool            x_log_scale         = R44C::x_log_scale;     
//     static constexpr            cusfloat        x_max               = R44C::x_max;           
//     static constexpr            cusfloat        x_min               = R44C::x_min;           
//     static constexpr            cusfloat        dx                  = R44C::dx;              
//     static constexpr            bool            y_log_scale         = R44C::y_log_scale;     
//     static constexpr            cusfloat        y_max               = R44C::y_max;           
//     static constexpr            cusfloat        y_min               = R44C::y_min;           
//     static constexpr            cusfloat        dy                  = R44C::dy;              
//     static constexpr            std::size_t     num_c               = R44C::num_c;        
//     static constexpr const      cusfloat*       coeffs              = R44C::c;            
//     static constexpr const      std::size_t*    ncx                 = R44C::ncx; 
//     static constexpr const      std::size_t*    ncy                 = R44C::ncy;          
// };    

#endif