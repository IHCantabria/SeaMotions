
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


// Include local modules
#include "farfield_integrator.hpp"
#include "../math/math_tools.hpp"


template<QTFTypeE qtf_type>
void FarFieldIntegrator<qtf_type>::compute_far_field(
                                                        std::size_t freq_index_i,
                                                        std::size_t freq_index_j
                                                    )
{
    // Set frequency indices and derived parameters
    this->set_frequency_indices( freq_index_i, freq_index_j );

    // Clear output vectors
    std::fill( this->qc.begin(), this->qc.end(), cuscomplex( 0.0, 0.0 ) );
    std::fill( this->qs.begin(), this->qs.end(), cuscomplex( 0.0, 0.0 ) );

    // Loop over series size to compute the far field contribution for each term of the expansion
    for ( std::size_t n=0; n<this->_series_size; n++ )
    {
        // Compute the far field contribution for the current term of the expansion
        this->_integrate( n, this->qc.data(), this->qs.data() );

    }

    // Apply scaling factor to the far field contribution
    for ( std::size_t ih0=0; ih0<static_cast<std::size_t>(this->_input->heads_np); ih0++ )
    {
        for ( std::size_t ih1=0; ih1<static_cast<std::size_t>(this->_input->heads_np); ih1++ )
        {
            // Determine output index
            std::size_t ihx = this->_input->heads_np * ih0 + ih1;

            // Apply scaling factor
            this->qc[ihx] *= cuscomplex( 0.0, PI * pow2s( this->_partition_circle ) / 8.0 );
            this->qs[ihx] *= cuscomplex( 0.0, PI * pow2s( this->_partition_circle ) / 8.0 );

        }

    }

}


template<QTFTypeE qtf_type>
FarFieldIntegrator<qtf_type>::FarFieldIntegrator(
                                                    Input*          input,
                                                    std::size_t     freq_i,
                                                    std::size_t     freq_j,
                                                    cusfloat        partition_circle,
                                                    cuscomplex*     Ac,
                                                    cuscomplex*     As,
                                                    cuscomplex*     Bc,
                                                    cuscomplex*     Bs,
                                                    std::size_t     series_size
                                                )
{
    // Check QTF type
    this->_hankel_int_opts.hkind0   = 1;
    if ( qtf_type == QTFTypeE::QTF_DIFF_CODE )
    {
        this->_sign                     = -1.0;
        this->_hankel_int_opts.hkind1   = 1;
    }
    else if ( qtf_type == QTFTypeE::QTF_SUM_CODE )
    {
        this->_sign                     = +1.0;
        this->_hankel_int_opts.hkind1   = 2;
    }
    else
    {
        throw std::invalid_argument("Invalid QTF type. Must be either QTF_DIFF_CODE or QTF_SUM_CODE.");
    }

    // Storage general usage parameters
    this->_ac               = Ac;
    this->_as               = As;
    this->_bc               = Bc;
    this->_bs               = Bs;
    this->_partition_circle = partition_circle;
    this->_series_size      = series_size;

    // Resize output vectors
    this->_hnp2 = this->_input->heads_np * this->_input->heads_np;
    this->qc.resize( this->hnp2, cuscomplex( 0.0, 0.0 ) );
    this->qs.resize( this->hnp2, cuscomplex( 0.0, 0.0 ) );

    // Storage input paramets and recalculation of derived parameters
    this->set_parameters( freq_i, freq_j );

}


template<QTFTypeE qtf_type>
void FarFieldIntegrator<qtf_type>::_integrate(
                                                std::size_t     n,
                                                cuscomplex*     Qc,
                                                cuscomplex*     Qs
                                            )
{
    // Declare local variables for integration
    cuscomplex  acil, acjm, asil, asjm;
    cuscomplex  bcil, bcjm, bsil, bsjm;
    std::size_t ihx, il, jm, m;
    cusfloat    kd;
    cuscomplex  tc1, tc2, tc3;
    cuscomplex  ts1, ts2, ts3;
    cuscomplex  FF, FG, FH;
    auto conj_if = []( const auto& v ) {
        if constexpr ( qtf_type == QTFTypeE::QTF_DIFF_CODE )
        {
            return std::conj( v );
        }

        return v;
    };

    // Sumate expansion series terms

    for ( std::size_t ih0=0; ih0<static_cast<std::size_t>(this->_input->heads_np); ih0++ )
    {
        for ( std::size_t ih1=0; ih1<static_cast<std::size_t>(this->_input->heads_np); ih1++ )
        {
            tc1 = cuscomplex( 0.0, 0.0 ); ts1 = cuscomplex( 0.0, 0.0 );
            tc2 = cuscomplex( 0.0, 0.0 ); ts2 = cuscomplex( 0.0, 0.0 );
            tc3 = cuscomplex( 0.0, 0.0 ); ts3 = cuscomplex( 0.0, 0.0 );

            ihx = this->_input->heads_np * ih0 + ih1;

            // First chunk l=n => M
            for ( std::size_t l=n; l<this->_series_size; l++ )
            {
                m       = l - n;

                il      = this->_freq_i_off + ih0 * QTF_FAR_N + l;
                jm      = this->_freq_j_off + ih1 * QTF_FAR_N + m;

                acil    = this->_ac[ il ];
                asil    = this->_as[ il ];
                acjm    = this->_ac[ jm ];
                asjm    = this->_as[ jm ];
                bcil    = this->_bc[ il ];
                bsil    = this->_bs[ il ];
                bcjm    = this->_bc[ jm ];
                bsjm    = this->_bs[ jm ];

                FF      = this->U( triple_hankel_f, l, m );
                FG      = this->U( triple_hankel_g, l, m );
                FH      = this->U( triple_hankel_h, l, m );
                
                tc1     += FF *( bcil * conj_if( bcjm ) + bsil * conj_if( bsjm ) );
                tc1     += FG *( acil * conj_if( bcjm ) + asil * conj_if( bsjm ) );
                tc1     += FH *( bcil * conj_if( acjm ) + bsil * conj_if( asjm ) );

                ts1     += FF *( bsil * conj_if( bcjm ) - bcil * conj_if( bsjm ) );
                ts1     += FG *( asil * conj_if( bcjm ) - acil * conj_if( bsjm ) );
                ts1     += FH *( bsil * conj_if( acjm ) - bcil * conj_if( asjm ) );

                kd      = static_cast<cusfloat>( l == n );
                Qc[ihx] += ( 1.0 + kd ) * ( tc1 + tc2 + tc3 );
                Qs[ihx] += ( 1.0 + kd ) * ( ts1 + ts2 + ts3 );

            }

            // Second chunk l=0 => M-n
            for ( std::size_t l=0; l<this->_series_size-n; l++ )
            {
                m       = l + n;

                il      = this->_freq_i_off + ih0 * QTF_FAR_N + l;
                jm      = this->_freq_j_off + ih1 * QTF_FAR_N + m;

                acil    = this->_ac[ il ];
                asil    = this->_as[ il ];
                acjm    = this->_ac[ jm ];
                asjm    = this->_as[ jm ];
                bcil    = this->_bc[ il ];
                bsil    = this->_bs[ il ];
                bcjm    = this->_bc[ jm ];
                bsjm    = this->_bs[ jm ];

                FF      = this->V( triple_hankel_f, l, m );
                FG      = this->V( triple_hankel_g, l, m );
                FH      = this->V( triple_hankel_h, l, m );
                
                tc1     += FF *( bcil * conj_if( bcjm ) + bsil * conj_if( bsjm ) );
                tc1     += FG *( acil * conj_if( bcjm ) + asil * conj_if( bsjm ) );
                tc1     += FH *( bcil * conj_if( acjm ) + bsil * conj_if( asjm ) );

                ts1     += FF *( bsil * conj_if( bcjm ) - bcil * conj_if( bsjm ) );
                ts1     += FG *( asil * conj_if( bcjm ) - acil * conj_if( bsjm ) );
                ts1     += FH *( bsil * conj_if( acjm ) - bcil * conj_if( asjm ) );

                kd      = static_cast<cusfloat>( l == 0 );
                Qc[ihx] += ( 1.0 + kd ) * ( tc1 + tc2 + tc3 );
                Qs[ihx] += ( 1.0 + kd ) * ( ts1 + ts2 + ts3 );

            }

            // Second chunk l=1 => n-1
            for ( std::size_t l=1; l<n-1; l++ )
            {
                m       = n - l;

                il      = this->_freq_i_off + ih0 * QTF_FAR_N + l;
                jm      = this->_freq_j_off + ih1 * QTF_FAR_N + m;

                acil    = this->_ac[ il ];
                asil    = this->_as[ il ];
                acjm    = this->_ac[ jm ];
                asjm    = this->_as[ jm ];
                bcil    = this->_bc[ il ];
                bsil    = this->_bs[ il ];
                bcjm    = this->_bc[ jm ];
                bsjm    = this->_bs[ jm ];

                FF      = this->W( triple_hankel_f, l, m );
                FG      = this->W( triple_hankel_g, l, m );
                FH      = this->W( triple_hankel_h, l, m );
                
                // First chunk l=n => M
                tc1     += FF *( bcil * conj_if( bcjm ) - bsil * conj_if( bsjm ) );
                tc1     += FG *( acil * conj_if( bcjm ) - asil * conj_if( bsjm ) );
                tc1     += FH *( bcil * conj_if( acjm ) - bsil * conj_if( asjm ) );

                ts1     += FF *( bsil * conj_if( bcjm ) + bcil * conj_if( bsjm ) );
                ts1     += FG *( asil * conj_if( bcjm ) + acil * conj_if( bsjm ) );
                ts1     += FH *( bsil * conj_if( acjm ) + bcil * conj_if( asjm ) );

                Qc[ihx] += ( tc1 + tc2 + tc3 );
                Qs[ihx] += ( ts1 + ts2 + ts3 );

            }

        }

    }

}


template<QTFTypeE qtf_type>
void FarFieldIntegrator<qtf_type>::set_frequency_indices( 
                                                            std::size_t freq_i,
                                                            std::size_t freq_j
                                                        )
{
    // Storage of input parameters
    this->_freq_i   = freq_i;
    this->_freq_j   = freq_j;
    this->_wi       = this->_input->angfreqs[freq_i];
    this->_wj       = this->_input->angfreqs[freq_j];
    this->_wave_i   = this->_input->wave_numbers[freq_i];
    this->_wave_j   = this->_input->wave_numbers[freq_j];
    this->_wave_k   = w2k( this->_wi + this->_sign * this->_wj, this->_input->grav_acc, this->_input->water_depth );

    // Recalculate derived parameters
    cusfloat nu_i   = pow2s( this->_wi ) / this->_input->grav_acc;
    cusfloat nu_j   = pow2s( this->_wj ) / this->_input->grav_acc;
    this->_omega    = ( 
                            + this->_wi * ( pow2s( this->_wave_j ) - pow2s( nu_j ) )
                            + this->_sign * this->_wj * ( pow2s( this->_wave_i ) - pow2s( nu_i ) )
                            - 2.0 * ( this->_wi + this->_sign * this->_wj ) * nu_i * nu_j
                        );
    this->_alpha    = ( this->_wi + this->_sign * this->_wj ) * this->_wave_i * this->_wave_j;

    // Calculate frequency offsets for expansion series coefficients
    this->_freq_i_off   = this->_freq_i * this->_input->heads_np * QTF_FAR_N;
    this->_freq_j_off   = this->_freq_j * this->_input->heads_np * QTF_FAR_N;

}


template<QTFTypeE qtf_type>
template<typename T>
cuscomplex FarFieldIntegrator<qtf_type>::U( 
                                                T                       F,
                                                int                     l,
                                                int                     m
                                            )
{
    // Set combined order
    int n = l - m;

    // Calculate l left coefficients
    cuscomplex  fl  = F( 
                            l-1,
                            m-1,
                            n,
                            this->_wave_i * this->_partition_circle,
                            this->_wave_j * this->_partition_circle,
                            this->_wave_k * this->_partition_circle,
                            this->_hankel_int_opts
                        );
    
    // Calculate l center coefficients
    cuscomplex  fc  = F( 
                            l,
                            m,
                            n,
                            this->_wave_i * this->_partition_circle,
                            this->_wave_j * this->_partition_circle,
                            this->_wave_k * this->_partition_circle,
                            this->_hankel_int_opts
                        );

    // Calculate l right coefficients
    cuscomplex  fr  = F( 
                            l+1,
                            m+1,
                            n,
                            this->_wave_i * this->_partition_circle,
                            this->_wave_j * this->_partition_circle,
                            this->_wave_k * this->_partition_circle,
                            this->_hankel_int_opts
                        );

    // Calculate composed formula for U coefficient
    return this->_omega * fc - this->_alpha * ( fl + fr );

}


template<QTFTypeE qtf_type>
template<typename T>
cuscomplex FarFieldIntegrator<qtf_type>::V( 
                                                T                       F,
                                                int                     l,
                                                int                     m
                                            )
{
    // Set combined order
    int n = m - l;

    // Calculate l left coefficients
    cuscomplex  fl  = F( 
                            l-1,
                            m-1,
                            n,
                            this->_wave_i * this->_partition_circle,
                            this->_wave_j * this->_partition_circle,
                            this->_wave_k * this->_partition_circle,
                            this->_hankel_int_opts
                        );
    
    // Calculate l center coefficients
    cuscomplex  fc  = F( 
                            l,
                            m,
                            n,
                            this->_wave_i * this->_partition_circle,
                            this->_wave_j * this->_partition_circle,
                            this->_wave_k * this->_partition_circle,
                            this->_hankel_int_opts
                        );

    // Calculate l right coefficients
    cuscomplex  fr  = F( 
                            l+1,
                            m+1,
                            n,
                            this->_wave_i * this->_partition_circle,
                            this->_wave_j * this->_partition_circle,
                            this->_wave_k * this->_partition_circle,
                            this->_hankel_int_opts
                        );

    // Calculate composed formula for U coefficient
    return this->_omega * fc - this->_alpha * ( fl + fr );

}


template<QTFTypeE qtf_type>
template<typename T>
cuscomplex FarFieldIntegrator<qtf_type>::W( 
                                                T                       F,
                                                int                     l,
                                                int                     m
                                            )
{
    // Set combined order
    int n = l + m;

    // Calculate l left coefficients
    cuscomplex  fl  = F( 
                            l-1,
                            m-1,
                            n,
                            this->_wave_i * this->_partition_circle,
                            this->_wave_j * this->_partition_circle,
                            this->_wave_k * this->_partition_circle,
                            this->_hankel_int_opts
                        );
    
    // Calculate l center coefficients
    cuscomplex  fc  = F( 
                            l,
                            m,
                            n,
                            this->_wave_i * this->_partition_circle,
                            this->_wave_j * this->_partition_circle,
                            this->_wave_k * this->_partition_circle,
                            this->_hankel_int_opts
                        );

    // Calculate l right coefficients
    cuscomplex  fr  = F( 
                            l+1,
                            m+1,
                            n,
                            this->_wave_i * this->_partition_circle,
                            this->_wave_j * this->_partition_circle,
                            this->_wave_k * this->_partition_circle,
                            this->_hankel_int_opts
                        );

    // Calculate composed formula for U coefficient
    return this->_omega * fc + this->_alpha * ( fl + fr );

}