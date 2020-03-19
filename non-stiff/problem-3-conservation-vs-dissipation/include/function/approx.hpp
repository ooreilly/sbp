#pragma once
#include "function/prec.hpp"
#include "function/dtoh.hpp"
#include "transform/reduce.hpp"
#include "functor/approx.hpp"
#include "macro/warn.hpp"


template <typename Tv=double, typename Tc=double>
__inline__ bool approx(const HostArray<Tv>& u, const HostArray<Tv>& v, const Bounds& bounds,
                       Tc prec = 0.0f, const bool verbose = false) {
        Approx<Tv> approx_fcn(prec, verbose);
        return reduce(approx_fcn, bounds, u, v, true);
}

template <typename Tv>
__inline__ bool approx(const DeviceArray<Tv>& u, const DeviceArray<Tv>& v, const Bounds& bounds,
                       Tv prec = 0.0f, const bool verbose = false) {
        warn(APPROX_DEVICE_NOT_IMPLEMENTED, APPROX_DEVICE_NOT_IMPLEMENTED_MSG);
        HostArray<Tv> h_u(u);
        HostArray<Tv> h_v(v);
        dtoh(h_u, u);
        dtoh(h_v, v);
        Approx<Tv> approx_fcn(prec, verbose);
        return reduce(approx_fcn, bounds, h_u, h_v, true);
}

template <typename Tv=double, typename Tc=double>
__inline__ bool approx(const DeviceArray<Tv>& u, const DeviceArray<Tv>& v,
                       Tc prec = 0.0f, const bool verbose = false) {
        return approx(u, v, u.bounds(), prec, verbose);
}


template<typename Tv>
__inline__ __host__ __device__ bool approx(Tv u, Tv v, Tv prec=0.0f, bool verbose=false) {
        return _approx(u, v, prec, verbose);
}
