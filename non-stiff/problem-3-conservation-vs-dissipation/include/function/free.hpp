#pragma once
#include "function/cuda_err_check.hpp"
#include "function/host_err_check.hpp"

hostError_t hostFree(void *out)
{
        if (out == NULL) return hostFreeFailure;
        free(out);
        return hostSuccess;
        
}


#define hfree(out) {hfree_((out), __FILE__, __LINE__); }
void hfree_(void *out, const char *file, int line) { 
        hostErrCheck_(hostFree(out), file, line);
}

#define dfree(stat) { dfree_((stat), __FILE__, __LINE__); }
void dfree_(void *out, const char *file, int line) { 
        if (out != NULL)
        cudaErrCheck_(cudaFree(out), file, line);
}

