#pragma once
#include "function/cuda_err_check.hpp"
#include "function/host_err_check.hpp"


hostError_t hostMalloc(void **out, size_t size)
{
        *out = malloc(size);
        if (out == NULL) return hostMallocFailure;
        return hostSuccess;
        
}

#define hmalloc(out, size) {hmalloc_((out), (size), __FILE__, __LINE__); }
void hmalloc_(void **out, const size_t size, const char *file, int line) { 
        hostErrCheck_(hostMalloc(out, size), file, line);
}

#define dmalloc(out, size) { dmalloc_((out), (size), __FILE__, __LINE__); }
void dmalloc_(void **out, size_t size, const char *file, int line) { 
        cudaErrCheck_(cudaMalloc(out, size), file, line);
}
