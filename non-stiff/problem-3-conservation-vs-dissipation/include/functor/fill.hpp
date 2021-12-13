#pragma once

#include "structure/array.hpp"
#include "function/malloc.hpp"

template <typename Tv>
class Fill {

        bool hasDeviceData = false;

        Tv *values;
        size_t nx, ny, nz; 

        Tv *u;
        size_t u_line, u_slice;

        public:
         Fill(DeviceArray<Tv>& u_, const Tv *values_, size_t nx_, size_t ny_ = 1,
              size_t nz_ = 1) {
                 size_t num_bytes = sizeof(Tv) * nx_ * ny_ * nz_;
                 cudaErrCheck(cudaMalloc((void **)&values, num_bytes));
                 cudaErrCheck(cudaMemcpy(values, values_, num_bytes,
                                         cudaMemcpyHostToDevice));
                 hasDeviceData = true;

                 nx = nx_;
                 ny = ny_;
                 nz = nz_;

                 u = u_();
                 u_line = u_.line;
                 u_slice = u_.slice;
         }

         Fill(HostArray<Tv> &u_, const Tv *values_, size_t nx_, size_t ny_ = 1,
              size_t nz_ = 1) {
                 size_t num_bytes = sizeof(Tv) * nx_ * ny_ * nz_;
                 hostErrCheck(hostMalloc((void **)&values, num_bytes));
                 for (int i = 0; i < nx_ * ny_ * nz_; i++) values[i] = values_[i];

                 nx = nx_;
                 ny = ny_;
                 nz = nz_;

                 u = u_();
                 u_line = u_.line;
                 u_slice = u_.slice;
         }

         VOID_OPERATOR() const {
                 size_t pos = xi + u_line * yj + u_slice * zk;
                 u[pos] = values[zk + yj * nz + xi * nz * ny]; 
        }

        ~Fill() {
                if (hasDeviceData) cudaErrCheck(cudaFree(values));
        }
};
