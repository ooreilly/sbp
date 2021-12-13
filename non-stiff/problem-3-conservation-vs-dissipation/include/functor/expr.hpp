#pragma once
#include "macro/operator.hpp"
#include "function/xpow.hpp"

// u = a * u + b
template <template <typename> class Array, typename Tv, typename Tc>
class ExprInplace {
        Tc a, b;

        Tv *u;
        size_t line, slice;

       public:
        ExprInplace(Array<Tv>& u_, Tc a_ = 1.0, Tc b_ = 1.0) {
                a = a_;
                b = b_;
                u = u_();
                line = u_.line;
                slice = u_.slice;
        }

        VOID_OPERATOR() const {
                size_t pos = xi  + yj * line + zk * slice;
                u[pos] =  a * u[pos] + b;
        }
};

// v = a * u + b * v + c
template <template <typename> class Array, typename Tv, typename Tc>
class Expr2{
        Tc a, b, c;

        const Tv *u; 
        size_t u_line, u_slice;

        Tv *v;
        size_t v_line, v_slice;

       public:
        Expr2(Array<Tv>& v_, const Array<Tv>& u_, Tc a_ = 1.0, Tc b_ = 1.0,
              Tc c_ = 1.0) {
                a = a_;
                b = b_;
                c = c_;
                u = u_();
                v = v_();
                u_line = u_.line;
                u_slice = u_.slice;
                v_line = v_.line;
                v_slice = v_.slice;
        }

        VOID_OPERATOR() const {
                const size_t u_pos = xi  + yj * u_line + zk * u_slice;
                const size_t v_pos = xi  + yj * v_line + zk * v_slice;
               v[v_pos] =  a * u[u_pos] + b * v[v_pos] + c;
        }
};

// w = a * u + b * v + c * w + d
template <template <typename> class Array, typename Tv, typename Tc>
class Expr3{
        Tc a, b, c, d;

        const Tv *u; 
        size_t u_line, u_slice;

        const Tv *v;
        size_t v_line, v_slice;

        Tv *w;
        size_t w_line, w_slice;

       public:
        Expr3(Array<Tv>& w_, const Array<Tv>& u_, const Array<Tv>& v_,
              Tc a_ = 1.0, Tc b_ = 0.0, Tc c_ = 0.0, Tc d_ = 0.0) {
                a = a_;
                b = b_;
                c = c_;
                d = d_;

                u = u_();
                u_line = u_.line;
                u_slice = u_.slice;

                v = v_();
                v_line = v_.line;
                v_slice = v_.slice;

                w = w_();
                w_line = w_.line;
                w_slice = w_.slice;
        }

        VOID_OPERATOR() const {
                const size_t u_pos = xi  + yj * u_line + zk * u_slice;
                const size_t v_pos = xi  + yj * v_line + zk * v_slice;
                const size_t w_pos = xi  + yj * w_line + zk * w_slice;
               w[w_pos] =  a * u[u_pos] + b * v[v_pos] + c * w[w_pos] + d;
        }
};
