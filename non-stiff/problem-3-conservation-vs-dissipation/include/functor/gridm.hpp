#pragma once
#include "macro/not_implemented.hpp"
#include "macro/dim.hpp"
#include "macro/operator.hpp"

// u = a * x_i + b * u + c, x_i = 0, (i - 1/2) * h, ..,  n - 1
template <typename Tv, typename Tc>
__inline__ __host__ __device__ Tv gridm(Tv val, Tv rm, size_t nr, Tc a = 1.0,
                                        Tc b = 0.0, Tc c = 1.0) {
        Tv z = 0.0;
        if (rm == 0)
                z = 0.0;
        else if (rm == nr - 1) 
                z = (nr - 2);
        else
                z = (rm - 0.5);
        return a * z + b * val + c;
}

template <int dim, template <typename> class Array, typename Tv, typename Tc>
class GridM {
        Tc a, b, c;

        Tv *u;
        size_t u_line, u_slice;

        size_t nx, ny, nz;

       public:
        GridM(Array<Tv>& u_, Tc a_ = 1.0, Tv b_ = 0.0, Tv c_ = 0.0) {
                a = a_;
                b = b_;
                c = c_;

                u = u_();
                u_line = u_.line;
                u_slice = u_.slice;

                nx = u_.nx;
                ny = u_.ny;
                nz = u_.nz;
        }

        VOID_OPERATOR() const {
                size_t nr = getdim<dim>(nx, ny, nz);
                Tv rm = getdim<dim>(xi, yj, zk);
                size_t pos = xi + u_line * yj + u_slice * zk;
                u[pos] =  gridm(u[pos], rm, nr, a, b, c);
        }
};

