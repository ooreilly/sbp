#pragma once
#include "macro/operator.hpp"

// v[bounds] = a * u + b * v[bounds] + c
template <template <typename> class Array, typename Tv, typename Tc>
class Scatter {
        Tc a, b, c;

        const Tv *u;
        size_t u_line, u_slice;

        Tv *v;
        size_t v_line, v_slice;

        Bounds bounds;

       public:
        Scatter(Array<Tv>& v_, const Array<Tv>& u_, const Bounds& bounds_,
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
                size_t u_pos = xi + u_line * yj + u_slice * zk;
                size_t v_pos = (xi + bounds.x[0]) +
                               v_line * (yj + bounds.y[0]) +
                               v_slice * (zk + bounds.z[0]);
                v[v_pos] = a * u[u_pos] + b * v[v_pos] + c;
        }
};

