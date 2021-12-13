#pragma once
#include "macro/operator.hpp"

// w = a * u .* v + b * w + c
template <template <typename> class Array, typename To, typename Ti, typename Tc>
class Multiply {
        Tc a, b, c;

        const Ti *u;
        size_t u_line, u_slice;

        const Ti *v;
        size_t v_line, v_slice;

        To *w;
        size_t w_line, w_slice;


       public:
        Multiply(Array<To>& w_, const Array<Ti>& u_, const Array<Ti>& v_,
                 Tc a_ = 1.0, Tc b_ = 0.0, Tc c_ = 0.0) {
                a = a_;
                b = b_;
                c = c_;

                w = w_();
                w_line = w_.line;
                w_slice = w_.slice;

                u = u_();
                u_line = u_.line;
                u_slice = u_.slice;

                v = v_();
                v_line = v_.line;
                v_slice = v_.slice;
        }

        VOID_OPERATOR() const {
                size_t u_pos = xi + u_line * yj + u_slice * zk;
                size_t v_pos = xi + v_line * yj + v_slice * zk;
                size_t w_pos = xi + w_line * yj + w_slice * zk;
                w[w_pos] =  a * u[u_pos] * v[v_pos] + b * w[w_pos] + c;
        }
};
