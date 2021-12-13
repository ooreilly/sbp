#pragma once

#define CUSTOM_OPERATOR(out_type, in_type, T)                                 \
        __inline__ __host__ __device__ out_type T(                            \
            const size_t xi, const size_t yj, const size_t zk,                \
            const in_type *ref1, const in_type *ref2, const out_type out_val, \
            const in_type in_val1, const in_type in_val2, const size_t nx,    \
            const size_t ny, const size_t nz, const size_t line,              \
            const size_t slice, const bool inplace)
#define OPERATOR() CUSTOM_OPERATOR(Tv, Tv, operator())
#define OPERATOR_TVIO() CUSTOM_OPERATOR(Tvo, Tvi, operator())

#define VOID_OPERATOR()                                 \
        __inline__ __host__ __device__ void operator()( \
            const size_t xi, const size_t yj, const size_t zk)

#define EVAL_OPERATOR(dim)                                                  \
        eval<(dim)>(xi, yj, zk, ref1, ref2, out_val, in_val1, in_val2, nx, ny, nz, \
                    line, slice, inplace)

#define SCALAR_OPERATOR(T)                                                    \
        __inline__ __host__ __device__ Tv T(                                  \
            const size_t xi, const size_t yj, const size_t zk,                \
            const Tv *ref1, const Tv *ref2,                                   \
            const Tv out_val, const Tv in_val1, const Tv in_val2,             \
            const size_t nx, const size_t ny, const size_t nz,                \
            const size_t line, const size_t slice, const bool inplace,        \
            const Tv *left, const Tv *right, const Tv *interior)

#define TENSOR_OPERATOR(T)                                                    \
        __inline__ __host__ __device__ Tv T(                                  \
            const size_t xi, const size_t yj, const size_t zk, const Tv *ref, \
            const Tv out_val, const Tv in_val1, const Tv in_val2,             \
            const size_t nx, const size_t ny, const size_t nz,                \
            const size_t line, const size_t slice, const bool inplace,        \
            const Tv *left1, const Tv *right1, const Tv *interior1,           \
            const Tv *left2, const Tv *right2, const Tv *interior2)

#define REDUCE_OPERATOR(T)                                                    \
        __inline__ __host__ __device__ T operator()(                          \
            const bool state, const size_t xi, const size_t yj,               \
            const size_t zk, const Tv* uref, const Tv uval, const size_t unx, \
            const size_t uny, const size_t unz, const size_t uline,           \
            const size_t uslice, const Tv* vref, const Tv vval,               \
            const size_t vnx, const size_t vny, const size_t vnz,             \
            const size_t vline, const size_t vslice)
