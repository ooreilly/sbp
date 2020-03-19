#pragma once
#define cudaErrCheck(stat) { cudaErrCheck_((stat), __FILE__, __LINE__); }
void cudaErrCheck_(cudaError_t stat, const char *file, int line) {
  if( stat == cudaErrorCudartUnloading ) {
          // Ignore this error
  }
  else if( stat != cudaSuccess) {                                                   
        fprintf(stderr, "CUDA error in %s:%i %s.\n",                          
          file, line, cudaGetErrorString(stat) );            
        fflush(stderr);                                                             
  }
}

