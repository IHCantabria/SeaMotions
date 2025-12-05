
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

#pragma once

// Include general usage libraries
#include "mpi.h"

// Include general usage scientific libraries
#include <complex>


template<typename T>
inline void cmpi_gather(
                                            T*              sendbuf,
                                            int             sendcount,
                                            T*              recvbuf,
                                            int             recvcount,
                                            int             root,
                                            MPI_Comm        comm
                )
{
    MPI_Gather( 
                    sendbuf,
                    sendcount,
                    MPI_FLOAT,
                    recvbuf,
                    recvcount,
                    MPI_FLOAT,
                    root,
                    comm
                );
}


template<>
inline void cmpi_gather<float>(
                                            float*          sendbuf,
                                            int             sendcount,
                                            float*          recvbuf,
                                            int             recvcount,
                                            int             root,
                                            MPI_Comm        comm
                        )
{
    MPI_Gather( 
                    sendbuf,
                    sendcount,
                    MPI_FLOAT,
                    recvbuf,
                    recvcount,
                    MPI_FLOAT,
                    root,
                    comm
                );
}


template<>
inline void cmpi_gather<double>(
                                            double*         sendbuf,
                                            int             sendcount,
                                            double*         recvbuf,
                                            int             recvcount,
                                            int             root,
                                            MPI_Comm        comm
                        )
{
    MPI_Gather( 
                    sendbuf,
                    sendcount,
                    MPI_DOUBLE,
                    recvbuf,
                    recvcount,
                    MPI_DOUBLE,
                    root,
                    comm
                );
}


template<>
inline void cmpi_gather<std::complex<float>>(
                                            std::complex<float>*    sendbuf,
                                            int                     sendcount,
                                            std::complex<float>*    recvbuf,
                                            int                     recvcount,
                                            int                     root,
                                            MPI_Comm                comm
                        )
{
    MPI_Gather( 
                    sendbuf,
                    sendcount,
                    MPI_COMPLEX8,
                    recvbuf,
                    recvcount,
                    MPI_COMPLEX8,
                    root,
                    comm
                );
}


template<>
inline void cmpi_gather<std::complex<double>>(
                                            std::complex<double>*   sendbuf,
                                            int                     sendcount,
                                            std::complex<double>*   recvbuf,
                                            int                     recvcount,
                                            int                     root,
                                            MPI_Comm                comm
                        )
{
    MPI_Gather( 
                    sendbuf,
                    sendcount,
                    MPI_COMPLEX16,
                    recvbuf,
                    recvcount,
                    MPI_COMPLEX16,
                    root,
                    comm
                );
}
