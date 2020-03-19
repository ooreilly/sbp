#pragma once
#include "macro/not_implemented.hpp"
#include "macro/dim.hpp"
#include "macro/operator.hpp"


// u = a * x_i + b * u + c, x_i = 0, 1, 2, ...
template <typename Tv, typename Tc>
__inline__ __host__ __device__ Tv gridp(Tv val, Tc rp, Tc a, Tc b,
                                         Tc c) {
        return a * rp + b * val + c;
}

template <int dim, template <typename> class Array, typename Tv, typename Tc>
class GridP {
        Tc a, b, c;

        Tv *u;
        size_t u_line, u_slice;

       public:
        GridP(Array<Tv>& u_, Tc a_ = 1.0, Tv b_ = 0.0, Tv c_ = 0.0) {
                a = a_;
                b = b_;
                c = c_;

                u = u_();
                u_line = u_.line;
                u_slice = u_.slice;
        }

        VOID_OPERATOR() const {
                Tc rp = getdim<dim>(xi, yj, zk);
                size_t pos = xi + u_line * yj + u_slice * zk;
                u[pos] = gridp(u[pos], rp, a, b, c);
        }
};

