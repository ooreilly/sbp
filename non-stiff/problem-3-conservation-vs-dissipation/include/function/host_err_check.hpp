#pragma once

enum hostError_num {
        hostSuccess,
        hostMallocFailure,
        hostFreeFailure,
        hostMemsetFailure
};
typedef enum hostError_num hostError_t;

const char *hostGetErrorString(hostError_t err)
{
        switch (err) {
        case hostSuccess:
                return "";
        case hostMallocFailure:
                return "Memory allocation failed";
        case hostFreeFailure:
                return "Memory deallocation failed";
        case hostMemsetFailure:
                return "Memset failed";
        }
        return "Unrecognized error";
}



#define hostErrCheck(stat) { hostErrCheck_((stat), __FILE__, __LINE__); }
void hostErrCheck_(hostError_t stat, const char *file, int line) {
  if( stat != hostSuccess) {                                                   
        fprintf(stderr, "Host error in %s:%i %s.\n",                          
          file, line, hostGetErrorString(stat) );            
        fflush(stderr);                                                             
  }
}

