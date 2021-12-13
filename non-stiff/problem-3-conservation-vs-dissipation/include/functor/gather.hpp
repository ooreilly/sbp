#pragma once
#include "macro/operator.hpp"

// v = a * u[bounds] + b * v + c
template <template <typename> class Array, typename Tv, typename Tc>
class Gather {
        Tc a, b, c;

        const Tv *u;
        size_t u_line, u_slice;

        Tv *v;
        size_t v_line, v_slice;

        Bounds bounds;

       public:
        Gather(Array<Tv>& v_, const Array<Tv>& u_, const Bounds& bounds_,
               Tc a_ = 1.0, Tc b_ = 0.0, Tc c_ = 0.0) {
                v = v_();
                v_line = v_.line;
                v_slice = v_.slice;

                u = u_();
                u_line = u_.line;
                u_slice = u_.slice;

                a = a_;
                b = b_;
                c = c_;

                bounds = bounds_;
        }

        VOID_OPERATOR() const {
                size_t u_pos = (xi + bounds.x[0]) +
                               u_line * (yj + bounds.y[0]) +
                               u_slice * (zk + bounds.z[0]);
                size_t v_pos = xi + v_line * yj + v_slice * zk;
                v[v_pos] = a * u[u_pos] + b * v[v_pos] + c;
        }
};

