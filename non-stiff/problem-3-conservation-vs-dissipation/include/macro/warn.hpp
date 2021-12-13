#pragma once

#ifndef APPROX_DEVICE_NOT_IMPLEMENTED
const bool APPROX_DEVICE_NOT_IMPLEMENTED = true;
#endif
const char *APPROX_DEVICE_NOT_IMPLEMENTED_MSG =
    "Not implemented for device arrays. Copying to host";

#define warn(warn_class, msg) \
        { warn_((warn_class), #warn_class, (msg), __FILE__, __LINE__); }
__inline__ __device__ __host__ void warn_(bool warn_state,
                                          const char *warn_class,
                                          const char *msg, const char *file,
                                          int line) {
        if (warn_state)
        printf("Warning %s %s:%d %s\n", warn_class, file, line, msg);
}

