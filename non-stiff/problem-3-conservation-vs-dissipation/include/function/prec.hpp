#pragma once

__inline__ __host__ __device__ float defaultprec(const float prec) {
        return 1e-6f;
}
__inline__ __host__ __device__ double defaultprec(const double prec) {
        return 1e-11;
}

