#pragma once

#define error(msg) \
        { error_((msg), __FILE__, __LINE__); }
__inline__ __device__ __host__ void error_(const char *msg, const char *file,
                                           int line) {
        printf("Error %s:%d %s\n", file, line, msg);
        assert(0);
}

