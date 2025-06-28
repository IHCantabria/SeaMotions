
#pragma once

// Include local modules
#include "../../src/config.hpp"
#include "chebyshev_traits.hpp"
#include "../math/chebyshev.hpp"
#include "../math/math_tools.hpp"


template<typename IDB>
inline void fold_database_1d( cusfloat H )
{
    // Check if z direction is log scaled
    cusfloat Hreg = H;
    if ( IDB::x_log_scale )
    {
        Hreg = std::log10( H );
    }

    // Check Z region bounds
    if ( Hreg > IDB::x_max_global )
    {
        Hreg = IDB::x_max_global;
    }

    if ( Hreg < IDB::x_min_global )
    {
        Hreg = IDB::x_min_global;
    }

    // Get starting region position
    std::size_t start_pos = 0;
    for ( std::size_t i=0; i<IDB::intervals_np; i++ )
    {
        if (
                ( IDB::x_min_region[i] < Hreg )
                &&
                ( IDB::x_max_region[i] > Hreg )
            )
        {
            start_pos = i;
            break;
        }
    }

    // Calculate chebyshev polynomials
    cusfloat    poly_h[ ( IDB::max_cheby_order + 1 ) ];
    cusfloat    h_map = 2.0 * ( Hreg - IDB::x_min_region[start_pos] ) / IDB::dx_region[start_pos] - 1.0;
    
    chebyshev_poly_upto_order( IDB::blocks_max_cheby_order[start_pos], h_map, poly_h );

    // Fold datase into a single scalar
    std::size_t bs = IDB::blocks_start[start_pos];
    std::size_t nt = IDB::blocks_coeffs_np[start_pos];
    ChebyshevTraits<IDB>::coeffs = 0.0;
    for ( std::size_t i=bs; i<bs+nt; i++ )
    {
        ChebyshevTraits<IDB>::coeffs += IDB::c[i] * poly_h[ IDB::ncx[i] ];
    }

}


template<typename IDB>
inline void fold_database_3d( cusfloat H )
{
    // Check if z direction is log scaled
    cusfloat Hreg = H;
    if ( IDB::z_log_scale )
    {
        Hreg = std::log10( H );
    }

    // Check Z region bounds
    if ( Hreg > IDB::z_max_global )
    {
        Hreg = IDB::z_max_global;
    }

    if ( Hreg < IDB::z_min_global )
    {
        Hreg = IDB::z_min_global;
    }

    // Get intersecting regions
    cusfloat    x_min_ints[IDB::blocks_np];
    cusfloat    x_max_ints[IDB::blocks_np];
    cusfloat    y_min_ints[IDB::blocks_np];
    cusfloat    y_max_ints[IDB::blocks_np];
    cusfloat    z_min_ints[IDB::blocks_np];
    cusfloat    z_max_ints[IDB::blocks_np];
    cusfloat    dx_ints[IDB::blocks_np];
    cusfloat    dy_ints[IDB::blocks_np];
    cusfloat    dz_ints[IDB::blocks_np];
    std::size_t blocks_start_ints[IDB::blocks_np];
    std::size_t blocks_coeffs_np_ints[IDB::blocks_np];
    std::size_t blocks_max_cheby_order_ints[IDB::blocks_np];

    std::size_t count_ints = 0;
    for ( std::size_t i=0; i<IDB::blocks_np; i++ )
    {
        if (
                ( IDB::z_min_region[i] < Hreg )
                &&
                ( IDB::z_max_region[i] > Hreg )
            )
        {
            x_min_ints[count_ints]                  = IDB::x_min_region[i];
            x_max_ints[count_ints]                  = IDB::x_max_region[i];
            y_min_ints[count_ints]                  = IDB::y_min_region[i];
            y_max_ints[count_ints]                  = IDB::y_max_region[i];
            z_min_ints[count_ints]                  = IDB::z_min_region[i];
            z_max_ints[count_ints]                  = IDB::z_max_region[i];
            dx_ints[count_ints]                     = IDB::dx_region[i];
            dy_ints[count_ints]                     = IDB::dy_region[i];
            dz_ints[count_ints]                     = IDB::dz_region[i];
            blocks_start_ints[count_ints]           = IDB::blocks_start[i];
            blocks_coeffs_np_ints[count_ints]       = IDB::blocks_coeffs_np[i];
            blocks_max_cheby_order_ints[count_ints] = IDB::blocks_max_cheby_order[i];
            // std::cout << count_ints << " " << i << std::endl;
            // std::cout << x_min_ints[count_ints] << std::endl;
            // std::cout << x_max_ints[count_ints] << std::endl;
            // std::cout << y_min_ints[count_ints] << std::endl;
            // std::cout << y_max_ints[count_ints] << std::endl;
            // std::cout << z_min_ints[count_ints] << std::endl;
            // throw std::exception( );
            count_ints++;
        }
    }

    // std::cout << "block_np: " << IDB::blocks_np;
    // std::cout << "block_np_f: " << IDB::blocks_np_f;
    // std::cout << "count_ints: " << count_ints << std::endl;
    // std::cout << x_min_ints[0] << std::endl;
    // std::cout << x_max_ints[0] << std::endl;
    // std::cout << y_min_ints[0] << std::endl;
    // std::cout << y_max_ints[0] << std::endl;
    // throw std::exception( );

    for ( std::size_t i=0; i<IDB::blocks_np; i++ )
    {
        if ( std::abs( IDB::x_max_region[i] - IDB::x_min_region[i] ) < 1e-6 )
        {
            std::cout << "XMin and XMax are equal!!" << std::endl;
            std::cout << "Index: " << i << " - XMin: " << IDB::x_min_region[i];
            std::cout << " - XMax: " << IDB::x_max_region[i] << std::endl;
            throw std::exception();
        }
    }

    // std::cout << "INTERSECTING REGIONS:" << std::endl;
    // for ( std::size_t i=0; i<count_ints; i++ )
    // {
    //     std::cout << "Index: " << i;
    //     std::cout << " XMin: " << x_min_ints[i];
    //     std::cout << " XMax: " << x_max_ints[i];
    //     std::cout << " YMin: " << y_min_ints[i];
    //     std::cout << " YMax: " << y_max_ints[i];
    //     std::cout << " ZMin: " << z_min_ints[i];
    //     std::cout << " ZMax: " << z_max_ints[i];
    //     std::cout << std::endl;
    // }

    // Copy selected region to an intermediate vector for manipulation
    int block_start_sel[IDB::blocks_np];
    for ( int i=0; i<IDB::blocks_np; i++ )
    {
        block_start_sel[i] = static_cast<int>( blocks_start_ints[i] );
    }

    // Check if there is further repetions of the selected regions
    for ( std::size_t i=0; i<count_ints; i++ )
    {
        if ( block_start_sel[i] > -1 )
        {
            for ( std::size_t j=i+1; j<count_ints; j++ )
            {
                if ( block_start_sel[i] == block_start_sel[j] )
                {
                    block_start_sel[j] = -1;
                }
            }
        }
    }

    // Discard repeated regions
    cusfloat    x_min_rep[IDB::blocks_np_f];
    cusfloat    x_max_rep[IDB::blocks_np_f];
    cusfloat    y_min_rep[IDB::blocks_np_f];
    cusfloat    y_max_rep[IDB::blocks_np_f];
    cusfloat    z_min_rep[IDB::blocks_np_f];
    cusfloat    z_max_rep[IDB::blocks_np_f];
    cusfloat    dx_rep[IDB::blocks_np_f];
    cusfloat    dy_rep[IDB::blocks_np_f];
    cusfloat    dz_rep[IDB::blocks_np_f];
    std::size_t blocks_start_rep[IDB::blocks_np_f];
    std::size_t blocks_coeffs_np_rep[IDB::blocks_np_f];
    std::size_t blocks_max_cheby_order_rep[IDB::blocks_np_f];

    std::size_t count_rep = 0;
    for ( std::size_t i=0; i<count_ints; i++ )
    {
        if ( block_start_sel[i] > -1 )
        {
            blocks_start_rep[count_rep]             = blocks_start_ints[i];
            blocks_coeffs_np_rep[count_rep]         = blocks_coeffs_np_ints[i];
            blocks_max_cheby_order_rep[count_rep]   = blocks_max_cheby_order_ints[i];
            x_min_rep[count_rep]                    = x_min_ints[i];
            x_max_rep[count_rep]                    = x_max_ints[i];
            y_min_rep[count_rep]                    = y_min_ints[i];
            y_max_rep[count_rep]                    = y_max_ints[i];
            z_min_rep[count_rep]                    = z_min_ints[i];
            z_max_rep[count_rep]                    = z_max_ints[i];
            dx_rep[count_rep]                       = dx_ints[i];
            dy_rep[count_rep]                       = dy_ints[i];
            dz_rep[count_rep]                       = dz_ints[i];
            count_rep++;

            if ( count_rep > IDB::blocks_np_f )
            {
                std::cout << "count_rep higher than blocks_np_f" << std::endl;
                throw std::exception( );
            }
        }
    }

    // Generate folded coefficients
    std::size_t blocks_start_fold[IDB::blocks_np_f];
    std::size_t blocks_coeffs_np_fold[IDB::blocks_np_f];
    std::size_t blocks_max_cheby_order_fold[IDB::blocks_np_f];
    cusfloat    poly_h[ ( IDB::max_cheby_order + 1 ) ];

    cusfloat    h_map           = 0.0;
    std::size_t count_fold      = 0;
    std::size_t local_count     = 0;
    std::size_t start_pos       = 0;
    std::size_t max_cheby_xy    = 0;
    for ( std::size_t i=0; i<count_rep; i++ )
    {
        // Map H to local coordinates of the region
        h_map = 2.0 * ( Hreg - z_min_rep[i] ) / dz_rep[i] - 1.0;

        // Calculate chebyshev polynomials up to the order
        chebyshev_poly_upto_order( blocks_max_cheby_order_rep[i], h_map, poly_h );

        // Fold third dimension
        ChebyshevTraits<IDB>::coeffs[count_fold]    = IDB::c[blocks_start_rep[i]+0]*poly_h[IDB::ncz[blocks_start_rep[i]+0]];
        ChebyshevTraits<IDB>::ncx[count_fold]       = IDB::ncx[blocks_start_rep[i]+0];
        ChebyshevTraits<IDB>::ncy[count_fold]       = IDB::ncy[blocks_start_rep[i]+0];
        local_count                                 = 1;
        max_cheby_xy                                = std::max( IDB::ncx[blocks_start_rep[i]+0], IDB::ncy[blocks_start_rep[i]+0] );
        for ( std::size_t j=1; j<blocks_coeffs_np_rep[i]; j++ )
        {
            std::cout << "ncx: " << IDB::ncx[blocks_start_rep[i]+j-1] << " - " << IDB::ncx[blocks_start_rep[i]+j];
            std::cout << " - ncy: " << IDB::ncy[blocks_start_rep[i]+j-1] << " - " << IDB::ncy[blocks_start_rep[i]+j] << std::endl;
            if ( 
                    ( IDB::ncx[blocks_start_rep[i]+j-1] != IDB::ncx[blocks_start_rep[i]+j] )
                    ||
                    ( IDB::ncy[blocks_start_rep[i]+j-1] != IDB::ncy[blocks_start_rep[i]+j] )
                )
            {
                std::cout << "Last Coeffcient Value: " << ChebyshevTraits<IDB>::coeffs[count_fold] << std::endl;
                // Start new coefficient
                count_fold++;
                local_count++;
                ChebyshevTraits<IDB>::coeffs[count_fold]    = IDB::c[blocks_start_rep[i]+j]*poly_h[IDB::ncz[blocks_start_rep[i]+j]];
                ChebyshevTraits<IDB>::ncx[count_fold]       = IDB::ncx[blocks_start_rep[i]+j];
                ChebyshevTraits<IDB>::ncy[count_fold]       = IDB::ncy[blocks_start_rep[i]+j];
            }
            else
            {
                ChebyshevTraits<IDB>::coeffs[count_fold]    += IDB::c[blocks_start_rep[i]+j]*poly_h[IDB::ncz[blocks_start_rep[i]+j]];
            }

            max_cheby_xy = std::max( max_cheby_xy, IDB::ncx[blocks_start_rep[i]+j] );
            max_cheby_xy = std::max( max_cheby_xy, IDB::ncy[blocks_start_rep[i]+j] );
            
        }
        // Storage previous results
        blocks_start_fold[i]            = start_pos;
        blocks_coeffs_np_fold[i]        = local_count;
        blocks_max_cheby_order_fold[i]  = max_cheby_xy;
        count_fold                      += 1;
        start_pos                       = count_fold;
    }

    std::cout << "BLOCK_START_FOLD:" << std::endl;
    for ( std::size_t i=0; i<count_rep; i++ )
    {
        std::cout << "Index: " << i;
        std::cout << " - " << blocks_start_fold[i] << std::endl;
    }

    std::cout << "BLOCK_COEFFS_NP_FOLD:" << std::endl;
    for ( std::size_t i=0; i<count_rep; i++ )
    {
        std::cout << "Index: " << i;
        std::cout << " - " << blocks_coeffs_np_fold[i] << std::endl;
    }

    std::cout << "MAX_CHEBY_ORDER_FOLD:" << std::endl;
    for ( std::size_t i=0; i<count_rep; i++ )
    {
        std::cout << "Index: " << i;
        std::cout << " - " << blocks_max_cheby_order_fold[i] << std::endl;
    }

    // Redistribute 
    cusfloat    dx      = ( IDB::x_max_global - IDB::x_min_global ) / IDB::intervals_np;
    cusfloat    dy      = ( IDB::y_max_global - IDB::y_min_global ) / IDB::intervals_np;
    int         index   = 0;
    cusfloat    xi      = 0.0;
    cusfloat    yi      = 0.0;

    // for ( std::size_t i=0; i<count_rep; i++ )
    // {
    //     std::cout << "Index: " << i << " - x: " << x_min_rep[i] << " - " << x_max_rep[i];
    //     std::cout << " - y: " << y_min_rep[i] << " - " << y_max_rep[i];
    //     std::cout << " - z: " << z_min_rep[i] << " - " << z_max_rep[i] << std::endl;
    // }

    for ( int i=0; i<IDB::intervals_np; i++ )
    {
        for ( int j=0; j<IDB::intervals_np; j++ )
        {
            // Calculate region properties
            index   = i * IDB::intervals_np + j;
            xi      = IDB::x_min_global + dx * ( 2 * i + 1 ) / 2.0;
            yi      = IDB::y_min_global + dy * ( 2 * j + 1 ) / 2.0;
            
            // Find region in between the folded ones
            bool is_found = false;
            for ( std::size_t k=0; k<count_rep; k++ )
            {
                if ( 
                        ( xi > x_min_rep[k] && xi < x_max_rep[k] )
                        &&
                        ( yi > y_min_rep[k] && yi < y_max_rep[k] )
                    )
                {
                    ChebyshevTraits<IDB>::blocks_start[index]           = blocks_start_fold[k];
                    ChebyshevTraits<IDB>::blocks_coeffs_np[index]       = blocks_coeffs_np_fold[k];
                    ChebyshevTraits<IDB>::blocks_max_cheby_order[index] = blocks_max_cheby_order_fold[k];
                    ChebyshevTraits<IDB>::x_min_region[index]           = x_min_rep[k];
                    ChebyshevTraits<IDB>::x_max_region[index]           = x_max_rep[k];
                    ChebyshevTraits<IDB>::dx_region[index]              = dx_rep[k];
                    ChebyshevTraits<IDB>::y_min_region[index]           = y_min_rep[k];
                    ChebyshevTraits<IDB>::y_max_region[index]           = y_max_rep[k];
                    ChebyshevTraits<IDB>::dy_region[index]              = dy_rep[k];
                    is_found                                            = true;
                }
            }
            if ( !is_found )
            {
                std::cout << "Not possible to find region -> x_mean: " << xi << " - y_mean: " << yi << std::endl;
            }
        }
    }
}