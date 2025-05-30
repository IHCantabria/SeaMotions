
#ifndef __chebyshev_traits_macros_hpp
#define __chebyshev_traits_macros_hpp

#define CHEBYSHEV_2D_TRAITS( NS )                                                               \
template<>                                                                                      \
struct ChebyshevTraits<NS>                                                                      \
{                                                                                               \
    static constexpr            int             max_ref_level       = NS::max_ref_level;        \
    static constexpr            int             intervals_np        = NS::intervals_np;         \
    static constexpr            int             max_cheby_order     = NS::max_cheby_order;      \
    static constexpr const      std::size_t*    blocks_start        = NS::blocks_start;         \
    static constexpr const      std::size_t*    blocks_coeffs_np    = NS::blocks_coeffs_np;     \
    static constexpr            bool            fcn_log_scale       = NS::fcn_log_scale;        \
    static constexpr            bool            x_log_scale         = NS::x_log_scale;          \
    static constexpr            cusfloat        x_max_global        = NS::x_max_global;         \
    static constexpr            cusfloat        x_min_global        = NS::x_min_global;         \
    static constexpr            cusfloat        dx_min_region       = NS::dx_min_region;        \
    static constexpr const      cusfloat*       x_max_region        = NS::x_max_region;         \
    static constexpr const      cusfloat*       x_min_region        = NS::x_min_region;         \
    static constexpr const      cusfloat*       dx_region           = NS::dx_region;            \
    static constexpr            bool            y_log_scale         = NS::y_log_scale;          \
    static constexpr            cusfloat        y_max_global        = NS::y_max_global;         \
    static constexpr            cusfloat        y_min_global        = NS::y_min_global;         \
    static constexpr            cusfloat        dy_min_region       = NS::dy_min_region;        \
    static constexpr const      cusfloat*       y_max_region        = NS::y_max_region;         \
    static constexpr const      cusfloat*       y_min_region        = NS::y_min_region;         \
    static constexpr const      cusfloat*       dy_region           = NS::dy_region;            \
    static constexpr            std::size_t     num_c               = NS::num_cf;               \
    static constexpr const      cusfloat*       coeffs              = NS::cf;                   \
    static constexpr const      std::size_t*    ncx                 = NS::ncxf;                 \
    static constexpr const      std::size_t*    ncy                 = NS::ncyf;                 \
}                                                                                               \

#endif