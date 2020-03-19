#pragma once
#include "structure/array.hpp"
#include "function/mpi_type.hpp"
#include "function/mpi_err_check.hpp"
#include "function/dtoh.hpp"

#define mpi_isend_world(buf, dest, tag, request) \
        mpi_err_check_(                           \
            mpi_isend_((buf), (dest), (tag), MPI_COMM_WORLD, (request)),  \
                       __FILE__, __LINE__, __func__)

#define mpi_isend(buf, dest, tag, comm, request) \
        mpi_err_check_(mpi_isend_((buf), dest, tag, comm, request),  \
                       __FILE__, __LINE__, __func__)
template <typename Tv>
int mpi_isend_(HostArray<Tv>& buf, int dest, int tag, MPI_Comm comm, 
              MPI_Request& request) {
        return MPI_Isend(buf.u, buf.size, mpi_type(buf.u), dest, tag, comm,
                         &request);
}

template <typename Tv>
int mpi_isend_(DeviceArray<Tv>& buf, int dest, int tag, MPI_Comm comm, 
              MPI_Request& request) {
#ifdef CUDA_AWARE_MPI
        return MPI_Isend(buf.u, buf.size, mpi_type(buf.u), dest, tag, comm,
                         &request);
#else
                cudaMemcpy(buf.v, buf.u, buf.num_bytes, cudaMemcpyDeviceToHost);
        return MPI_Isend(buf.v, buf.size, mpi_type(buf.u), dest, tag, comm,
                         &request);
#endif
}
