#pragma once
#include "structure/array.hpp"
#include "function/cuda_err_check.hpp"

#define dtoh(dest, src) { dtoh_( (dest), (src), __FILE__, __LINE__); }
template <typename Tv>
void dtoh_(HostArray<Tv>& dest, const DeviceArray<Tv>& src, const char* file,
           int line) {
        assert(src.num_bytes >= dest.num_bytes);
        cudaErrCheck_(
            cudaMemcpy(dest.u, src.u, dest.num_bytes, cudaMemcpyDeviceToHost),
            file, line);
}

