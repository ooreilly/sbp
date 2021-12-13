#pragma once
#include "macro/dim.hpp"

template <int dim, typename T>
__inline__ __host__ __device__ T getdim(T valx, T valy, T valz) {
        if (dim == X)
                return valx;
        if (dim == Y)
                return valy;
        if (dim == Z)
                return valz;
}

