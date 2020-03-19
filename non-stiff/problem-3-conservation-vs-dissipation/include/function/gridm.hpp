#pragma once
#include "transform/apply.hpp"
#include "functor/gridm.hpp"
#include "macro/dim.hpp"

// u = a * x + b * u + c
template <int dim, template <typename> class Array, typename Tv = double,
          typename Tc = double>
void gridm(Array<Tv>& u, const Bounds& bounds, Tc a = 1.0, Tc b = 0.0,
            Tc c = 0.0) {
        GridM<dim, Array, Tv, Tc> gm(u, a, b, c);
        apply(gm, bounds, u);
}

template <template <typename> class Array, typename Tv = double,
          typename Tc = double>
void xm(Array<Tv>& u, const Bounds& bounds, Tc a = 1.0, Tc b = 0.0,
        Tc c = 0.0) {
        gridm<X>(u, bounds, a, b, c);
}

template <template <typename> class Array, typename Tv = double,
          typename Tc = double>
void ym(Array<Tv>& u, const Bounds& bounds, Tc a = 1.0, Tc b = 0.0,
        Tc c = 0.0) {
        gridm<Y>(u, bounds, a, b, c);
}

template <template <typename> class Array, typename Tv = double,
          typename Tc = double>
void zm(Array<Tv>& u, const Bounds& bounds, Tc a = 1.0, Tc b = 0.0,
        Tc c = 0.0) {
        gridm<Z>(u, bounds, a, b, c);
}

// Overloads that do not specify the bounds argument
template <template <typename> class Array, typename Tv = double,
          typename Tc = double>
void xm(Array<Tv>& u, Tc a = 1.0, Tc b = 0.0, Tc c = 0.0) {
        gridm<X>(u, u.bounds(), a, b, c);
}

template <template <typename> class Array, typename Tv = double,
          typename Tc = double>
void ym(Array<Tv>& u, Tc a = 1.0, Tc b = 0.0, Tc c = 0.0) {
        gridm<Y>(u, u.bounds(), a, b, c);
}

template <template <typename> class Array, typename Tv = double,
          typename Tc = double>
void zm(Array<Tv>& u, Tc a = 1.0, Tc b = 0.0, Tc c = 0.0) {
        gridm<Z>(u, u.bounds(), a, b, c);
}
