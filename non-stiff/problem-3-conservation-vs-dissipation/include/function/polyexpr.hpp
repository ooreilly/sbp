#pragma once
#include "transform/apply.hpp"
#include "functor/polyexpr.hpp"

// w = a * u^b + c * v^d + e * w + f
template <template <typename> class Array, typename To = double,
          typename Ti = double, typename Tc = double>
void polyexpr(Array<To>& w, const Array<Ti>& u, const Array<Ti>& v,
              const Bounds& bounds, const Tc a = 1.0, const Tc b = 1.0,
              const Tc c = 0.0, const Tc d = 0.0, const Tc e = 0.0,
              const Tc f = 0.0) {
        PolyExpr2<Array, To, Ti, Tc> expr(w, u, v, a, b, c, d, e, f);
        return apply(expr, bounds, w);
}


// v = a * u^b + c * v + d
template <template <typename> class Array, typename To = double,
          typename Ti = double, typename Tc = double>
void polyexpr(Array<To>& w, const Array<Ti>& u, const Bounds& bounds,
              const Tc a = 1.0, const Tc b = 1.0, const Tc c = 0.0,
              const Tc d = 0.0) {
        PolyExpr1<Array, To, Ti, Tc> expr(w, u, a, b, c);
        apply(expr, bounds, u);
}

// u = a * u^b + c
template <template <typename> class Array, typename Tv = double,
          typename Tc = double>
void polyexpr(Array<Tv>& u, const Bounds& bounds, const Tc a = 1.0,
              const Tc b = 1.0, const Tc c = 0.0) {
        PolyExprInplace<Array, Tv, Tc> expr(u, a, b, c);
        apply(expr, bounds, u);
}


// Without bounds
template <template <typename> class Array, typename To = double,
          typename Ti = double, typename Tc = double>
void polyexpr(Array<To>& w, const Array<Ti>& u, const Array<Ti>& v,
              const Tc a = 1.0, const Tc b = 1.0,
              const Tc c = 0.0, const Tc d = 0.0, const Tc e = 0.0,
              const Tc f = 0.0) {
        polyexpr(w, u, v, w.bounds(), a, b, c, d, e, f);
}

template <template <typename> class Array, typename To = double,
          typename Ti = double, typename Tc = double>
void polyexpr(Array<To>& w, const Array<Ti>& u,
              const Tc a = 1.0, const Tc b = 1.0, const Tc c = 0.0,
              const Tc d = 0.0) {
        polexpr(w, u, w.bounds(), a, b, c, d);
}

template <template <typename> class Array, typename Tv = double,
          typename Tc = double>
void polyexpr(Array<Tv>& u, const Tc a = 1.0, const Tc b = 1.0,
              const Tc c = 0.0) {
        return polyexpr(u, u.bounds(), a, b, c);
}
