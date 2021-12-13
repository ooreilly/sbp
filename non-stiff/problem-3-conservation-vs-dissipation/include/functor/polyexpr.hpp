#pragma once
#include "macro/operator.hpp"
#include "function/xpow.hpp"



// w = a * u^p + c * w + d
template <typename Tv>
class PolyExpr {
        Tv a, p, c, d;

       public:
        PolyExpr(Tv a_ = 1.0, Tv p_ = 0.0, Tv c_ = 1.0, Tv d_ = 0.0) {
                a = a_;
                p = p_;
                c = c_;
                d = d_;
        }

        OPERATOR() const {
                return a * xpow(in_val1, p) + c * out_val + d;
        }
};

// w = a * u^p + b * v^q + c * w
template <typename Tv>
class PolyExprUV {
        Tv a, b, p, q, c;

       public:
        PolyExprUV(Tv a_ = 1.0, Tv b_ = 1.0, Tv p_ = 1.0, Tv q_ = 1.0,
                   Tv c_ = 0.0) {
                a = a_;
                b = b_;
                p = p_;
                q = q_;
                c = c_;
        }

        OPERATOR() const {
                return a * xpow(in_val1, p) + b * xpow(in_val2, q) + c * out_val;
        }
};

// w = a * u^b + c * v^d + e * w + f
template <template <typename> class Array, typename To, typename Ti, typename Tc>
class PolyExpr2 {
        Tc a, b, c, d, e, f;

        const Ti *u;
        size_t u_slice, u_line;

        const Ti *v;
        size_t v_slice, v_line;

        To *w;
        size_t w_slice, w_line;

       public:
        PolyExpr2(Array<To> &w_, const Array<To> &u_, const Array<Ti> &v_,
                  Tc a_ = 1.0, Tc b_ = 1.0, Tc c_ = 0.0, Tc d_ = 0.0,
                  Tc e_ = 0.0, Tc f_ = 0.0) {
                a = a_;
                b = b_;
                c = c_;
                d = d_;
                e = e_;
                f = f_;

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
                size_t u_pos = xi + u_line * yj + u_slice * zk;
                size_t v_pos = xi + v_line * yj + v_slice * zk;
                size_t w_pos = xi + w_line * yj + w_slice * zk;
                w[w_pos] = a * xpow((Tc)u[u_pos], b) +
                           c * xpow((Tc)v[v_pos], d) +
                           e * w[w_pos] + f;
        }
};

// v = a * u^b + c * v + d
template <template <typename> class Array, typename To, typename Ti, typename Tc>
class PolyExpr1 {
        Tc a, b, c, d;

        const Ti *u;
        size_t u_slice, u_line;

        To *v;
        size_t v_slice, v_line;

       public:
        PolyExpr1(Array<To> &v_, const Array<Ti> &u_, Tc a_ = 1.0, Tc b_ = 1.0,
                  Tc c_ = 0.0, Tc d_ = 0.0) {
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
        }

        VOID_OPERATOR() const {
                size_t u_pos = xi + u_line * yj + u_slice * zk;
                size_t v_pos = xi + v_line * yj + v_slice * zk;
                v[v_pos] = a * xpow((Tc)u[u_pos], b) + c * v[v_pos] + d;
        }
};

template <template <typename> class Array, typename Tv, typename Tc>
class PolyExprInplace {
        Tc a, b, c;
        Tv *u;
        size_t u_slice, u_line;

       public:
        PolyExprInplace(Array<Tv> &u_, Tc a_ = 1.0, Tc b_ = 1.0, Tc c_ = 0.0) {
                a = a_;
                b = b_;
                c = c_;

                u = u_();
                u_line = u_.line;
                u_slice = u_.slice;


        }

        VOID_OPERATOR() const {
                size_t pos = xi + u_line * yj + u_slice * zk;
                u[pos] = a * xpow((Tc)u[pos], b) + c;
        }
};
