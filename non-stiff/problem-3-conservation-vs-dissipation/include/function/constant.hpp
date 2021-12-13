#pragma once
#include "transform/apply.hpp"
#include "structure/array.hpp"
#include "functor/constant.hpp"
#include "macro/dim.hpp"


// v = c
template <template <typename> class Array, typename Tv, typename Tc>
void constant(Array<Tv>& v, const Bounds& bounds, const Tc c = 1.0) {
        Constant<Array, Tv> cnst(v, (Tv)c);
        apply(cnst, bounds, v);
}

template <template <typename> class Array, typename Tv, typename Tc>
void constant(Array<Tv>& v, const Tc c = 1.0) {
        return constant(v, v.bounds(), (Tv) c);
}
