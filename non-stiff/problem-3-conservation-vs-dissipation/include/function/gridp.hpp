#pragma once
#include "transform/apply.hpp"
#include "functor/gridp.hpp"
#include "macro/dim.hpp"

// u = a * x + b * u + c
template <int dim, template <typename> class Array, typename Tv = double,
          typename Tc = double>
void gridp(Array<Tv>& u, const Bounds& bounds, Tc a = 1.0, Tc b = 0.0,
            Tc c = 0.0) {
        GridP<dim, Array, Tv, Tc> gp(u, a, b, c);
        apply(gp, bounds, u);
}

template <template <typename> class Array, typename Tv = double,
          typename Tc = double>
void xp(Array<Tv>& u, const Bounds& bounds, Tc a = 1.0, Tc b = 0.0,
        Tc c = 0.0) {
        gridp<X>(u, bounds, a, b, c);
}

template <template <typename> class Array, typename Tv = double,
          typename Tc = double>
void yp(Array<Tv>& u, const Bounds& bounds, Tc a = 1.0, Tc b = 0.0,
        Tc c = 0.0) {
        gridp<Y>(u, bounds, a, b, c);
}

template <template <typename> class Array, typename Tv = double,
          typename Tc = double>
void zp(Array<Tv>& u, const Bounds& bounds, Tc a = 1.0, Tc b = 0.0,
        Tc c = 0.0) {
        gridp<Z>(u, bounds, a, b, c);
}

// Overloads that do not specify the bounds argument
template <template <typename> class Array, typename Tv = double,
          typename Tc = double>
void xp(Array<Tv>& u, Tc a = 1.0, Tc b = 0.0, Tc c = 0.0) {
        gridp<X>(u, u.bounds(), a, b, c);
}

template <template <typename> class Array, typename Tv = double,
          typename Tc = double>
void yp(Array<Tv>& u, Tc a = 1.0, Tc b = 0.0, Tc c = 0.0) {
        gridp<Y>(u, u.bounds(), a, b, c);
}

template <template <typename> class Array, typename Tv = double,
          typename Tc = double>
void zp(Array<Tv>& u, Tc a = 1.0, Tc b = 0.0, Tc c = 0.0) {
        gridp<Z>(u, u.bounds(), a, b, c);
}
