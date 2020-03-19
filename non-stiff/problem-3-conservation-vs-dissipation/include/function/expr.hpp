#pragma once
#include "transform/apply.hpp"
#include "structure/array.hpp"
#include "functor/expr.hpp"

// w = a * u + b * v + c * w
template <template <typename> class Array, typename Tv=double, typename Tc=double>
void expr(Array<Tv>& w, const Array<Tv>& u, const Array<Tv>& v, 
          const Bounds& bounds,
          const Tc a = 1.0, const Tc b = 1.0, 
          const Tc c = 0.0,
          const Tc d = 0.0) {
        Expr3<Array, Tv, Tc> expr(w, u, v, a, b, c, d);
        return apply(expr, bounds, w);
}


// v = a * u + b * v + c
template <template <typename> class Array, typename Tv=double, typename Tc=double>
void expr(Array<Tv>& v, const Array<Tv>& u, const Bounds& bounds,
        const Tc a = 1.0, const Tc b = 0.0, const Tc c = 0.0) {
        Expr2<Array, Tv, Tc> expr(v, u, a, b, c);
        apply(expr, bounds, v);
}

// f(u) = a * u + b
template <template <typename> class Array, typename Tv=double, typename Tc=double>
void expr(Array<Tv>& u, const Bounds& bounds,
        const Tc a = 1.0, const Tc b = 0.0) {
        ExprInplace<Array, Tv, Tc> expr(u, a, b);
        apply(expr, bounds, u);
}

template <template <typename> class Array, typename Tv=double, typename Tc=double>
void expr(Array<Tv>& u, const Tc a = 1.0, const Tc b = 0.0) {
        ExprInplace<Array, Tv, Tc> expr(u, a, b);
        apply(expr, u.bounds(), u);
}
