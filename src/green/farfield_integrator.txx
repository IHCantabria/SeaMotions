
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
FarFieldIntegrator<qtf_type>::FarFieldIntegrator(
                                                    Input*      input,
                                                    std::size_t freq_i,
                                                    std::size_t freq_j,
                                                    cusfloat    partition_circle
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
    this->_partition_circle = partition_circle;

    // Storage input paramets and recalculation of derived parameters
    set_parameters( freq_i, freq_j );

}


template<QTFTypeE qtf_type>
void FarFieldIntegrator<qtf_type>::integrate(
                                                std::size_t     n,
                                                std::size_t     N,
                                                cuscomplex      Ac,
                                                cuscomplex      As,
                                                cuscomplex      Bc,
                                                cuscomplex      Bs,
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

    // Clear vector 
    int hnp2 = this->_input->heads_np * this->_input->heads_np;
    clear_vector( hnp2, Qc );
    clear_vector( hnp2, Qs );

    // Sumate expansion series terms
    cuscomplex  sum = 0.0;

    for ( std::size_t ih0=0; ih0<this->_input->heads_np; ih0++ )
    {
        for ( std::size_t ih1=0; ih1<this->_input->heads_np; ih1++ )
        {
            tc1 = cuscomplex( 0.0, 0.0 ); ts1 = cuscomplex( 0.0, 0.0 );
            tc2 = cuscomplex( 0.0, 0.0 ); ts2 = cuscomplex( 0.0, 0.0 );
            tc3 = cuscomplex( 0.0, 0.0 ); ts3 = cuscomplex( 0.0, 0.0 );

            ihx = this->_input->heads_np * ih0 + ih1;

            // First chunk l=n => M
            for ( std::size_t l=n; l<N; l++ )
            {
                m       = l - n;

                il      = this->_freq_i_off + ih0 * QTF_FAR_N + l;
                jm      = this->_freq_j_off + ih1 * QTF_FAR_N + m;

                acil    = Ac[ il ];
                asil    = As[ il ];
                acjm    = Ac[ jm ];
                asjm    = As[ jm ];
                bcil    = Bc[ il ];
                bsil    = Bs[ il ];
                bcjm    = Bc[ jm ];
                bsjm    = Bs[ jm ];

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
            for ( std::size_t l=0; l<N-n; l++ )
            {
                m       = l + n;

                il      = this->_freq_i_off + ih0 * QTF_FAR_N + l;
                jm      = this->_freq_j_off + ih1 * QTF_FAR_N + m;

                acil    = Ac[ il ];
                asil    = As[ il ];
                acjm    = Ac[ jm ];
                asjm    = As[ jm ];
                bcil    = Bc[ il ];
                bsil    = Bs[ il ];
                bcjm    = Bc[ jm ];
                bsjm    = Bs[ jm ];

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

                acil    = Ac[ il ];
                asil    = As[ il ];
                acjm    = Ac[ jm ];
                asjm    = As[ jm ];
                bcil    = Bc[ il ];
                bsil    = Bs[ il ];
                bcjm    = Bc[ jm ];
                bsjm    = Bs[ jm ];

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
    cusfloat nu_i   = pow2s( this->_wi ) / this->_input->_grav_acc;
    cusfloat nu_j   = pow2s( this->_wj ) / this->_input->_grav_acc;
    this->_omega    = ( 
                            + this->_wi * ( pow2s( this->_wave_j ) - pow2s( nu_j ) )
                            + this->_sign * this->_wj * ( pow2s( this->_wave_i ) - pow2s( nu_i ) )
                            - 2.0 * ( this->_wi + this->_sign * this->_wj ) * nu_i * nu_j
                        );
    this->_alpha    = ( this->_wi + this->_sign * this->_wj ) * this->_wave_i * this->_wave_j;

    // Calculate frequency offsets for expansion series coefficients
    this->_freq_i_off   = this->_freq_i * this->input->heads_np * QTF_FAR_N;
    this->_freq_j_off   = this->_freq_j * this->input->heads_np * QTF_FAR_N;

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