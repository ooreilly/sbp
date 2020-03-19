#pragma once
#include "function/memset.hpp"

#define hzero(out, size) { hzero_((out), (size), __FILE__, __LINE__); }
void hzero_(void *out, size_t size, const char *file, int line) { 
        hostErrCheck_(hostMemset(out, 0, size), file, line);
}

#define dzero(out, size) { dzero_((out), (size), __FILE__, __LINE__); }
void dzero_(void *out, size_t size, const char *file, int line) { 
        cudaErrCheck_(cudaMemset(out, 0, size), file, line);
}

