#pragma once
#include "function/cuda_err_check.hpp"

#define htod(dest, src) { htod_( (dest), (src), __FILE__, __LINE__); }
template <typename Ta, typename Tb>
void htod_(Ta& dest, Tb& src, const char *file, int line) {
        assert(dest.num_bytes >= src.num_bytes);
        cudaErrCheck_(cudaMemcpy(dest.u, src.u, dest.num_bytes,
                                 cudaMemcpyHostToDevice),
                      file, line);
}


