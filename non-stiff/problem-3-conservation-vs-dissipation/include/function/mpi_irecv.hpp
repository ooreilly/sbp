#pragma once
#include "structure/array.hpp"
#include "function/mpi_type.hpp"
#include "function/mpi_err_check.hpp"
#include "function/dtoh.hpp"

#define mpi_irecv_world(buf, src, tag, request) \
        mpi_err_check_(                           \
            mpi_irecv_((buf), (src), (tag), MPI_COMM_WORLD, (request)),  \
                       __FILE__, __LINE__, __func__)

#define mpi_irecv(buf, src, tag, comm, request)                            \
        mpi_err_check_(mpi_irecv_((buf), (src), (tag), (comm), (request)), \
                       __FILE__, __LINE__, __func__)
template <typename Tv>
int mpi_irecv_(const HostArray<Tv>& buf, int src, int tag,
              MPI_Comm comm, MPI_Request& request) {
        return MPI_Irecv(buf.u, buf.size, mpi_type(buf.u), src, tag, comm,
                         &request);
}

template <typename Tv>
int mpi_irecv_(const DeviceArray<Tv>& buf, int src, int tag,
              MPI_Comm comm, MPI_Request& request) {
#ifdef CUDA_AWARE_MPI
        return MPI_Irecv(buf.u, buf.size, mpi_type(buf.u), src, tag, comm,
                         &request);
#else
        return MPI_Irecv(buf.v, buf.size, mpi_type(buf.u), src, tag, comm,
                         &request);
#endif
}
