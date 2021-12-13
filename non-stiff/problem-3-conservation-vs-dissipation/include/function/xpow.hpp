#pragma once

// Define pow function such that 0^0 = 1
template <typename Tv>
__inline__ __host__ __device__ Tv xpow(Tv a, Tv p) {
        if (a == (Tv)0 && p == (Tv)0) return (Tv)1.0;
        else return pow(a, p);
}

