#pragma once


__inline__ __device__ __host__ void notImplemented(const char *msg, const char *file, int line, int error=1) {
        const char *stat = error == 1 ? "Error" : "Warning";
        if (msg)
                printf("%s %s:%d %s\n", stat, file, line, msg);
        else
                printf("%s %s:%d Not implemented. \n", stat, file, line);
        if (error) assert(0);
}

#define deviceNotImplemented(msg) {deviceNotImplemented_((msg), __FILE__, __LINE__, 1);}
#define warnDeviceNotImplemented(msg) {deviceNotImplemented_((msg), __FILE__, __LINE__, 0);}
__inline__ __device__ __host__ void deviceNotImplemented_(const char *msg, const char *file, int line, int error) {
#ifdef __CUDA_ARCH__
        notImplemented(msg, file, line, error);
#endif
}

#define hostNotImplemented(msg) {hostNotImplemented_((msg), __FILE__, __LINE__, 1);}
#define warnHostNotImplemented(msg) {hostNotImplemented_((msg), __FILE__, __LINE__, 0);}
__inline__ __device__ __host__ void hostNotImplemented_(const char *msg, const char *file, int line, int error) {
#ifndef __CUDA_ARCH__
        notImplemented(msg, file, line, error);
#endif
}

#define notImplemented(msg) {notImplemented_((msg), __FILE__, __LINE__, 1);}
#define warnNotImplemented(msg) {notImplemented_((msg), __FILE__, __LINE__, 0);}
__inline__ __device__ __host__ void notImplemented_(const char *msg, const char *file, int line, int error) {
        deviceNotImplemented_(msg, file, line, error);
        hostNotImplemented_(msg, file, line, error);
}
