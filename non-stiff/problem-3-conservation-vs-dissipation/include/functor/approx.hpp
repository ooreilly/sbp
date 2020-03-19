#pragma once

#include "macro/not_implemented.hpp"
#include "macro/operator.hpp"

template <typename Tv>
__inline__ __host__ __device__ bool _approx(Tv u, Tv v, Tv prec = 0.0f,
                                            bool verbose = false) {
        if (prec == 0.0f) prec = defaultprec(prec);

        Tv value = fabs(u - v);
        bool isclose = value <= prec;

        if (!isclose && verbose) printf("approx(u, v) = %g > %g", value, prec);

        return isclose;
}

template <typename Tv>
class Approx {
        bool verbose;
        Tv prec;

       public:
        Approx(Tv prec_ = 0.0f, bool verbose_ = false) {
                if (prec == 0.0f) prec = defaultprec(prec);
                verbose = verbose_;
                prec = prec_;
        }

        REDUCE_OPERATOR(bool) const {
                deviceNotImplemented("");
                bool isclose = _approx(uval, vval, prec, verbose);
                if (!isclose && verbose)
                        printf(" at (%ld, %ld, %ld)\n", xi, yj, zk);
                return state && isclose; 
        }
};

