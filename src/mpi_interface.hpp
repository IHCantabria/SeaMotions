
#ifndef __mpi_interface_hpp
#define __mpi_interface_hpp

// Include general usage libraries
#include "mpi.h"

// Include general usage scientific libraries
#include <complex>


template<typename T>
void cmpi_gather(
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
void cmpi_gather<float>(
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
void cmpi_gather<double>(
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
void cmpi_gather<std::complex<float>>(
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
void cmpi_gather<std::complex<double>>(
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


#endif