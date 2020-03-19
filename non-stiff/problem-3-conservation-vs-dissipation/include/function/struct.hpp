#pragma once
#include "function/cuda_err_check.hpp"
#include "function/host_err_check.hpp"

// Copy struct from host to device
#define dstruct(out, in) \
        { dstruct_((out), (in), __FILE__, __LINE__); }
template <typename P>
void dstruct_(P **device, P &host, const char *file, int line) {
        cudaErrCheck_(cudaMalloc((void **)device, sizeof(host)), file, line);
        cudaErrCheck_(
            cudaMemcpy(*device, &host, sizeof(host), cudaMemcpyHostToDevice),
            file, line);
}


