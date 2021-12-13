#pragma once
#include "transform/apply.hpp"
#include "structure/array.hpp"
#include "functor/constant.hpp"
// f(u) = 0

template <template <typename> class Array, typename Tv>
void zero(Array<Tv>& u, const Bounds& bounds)
{
        Constant<Array, Tv> z(u, (Tv)0.0);
        return apply(z, bounds, u);
}

template <template <typename> class Array, typename Tv>
void zero(Array<Tv>& u)
{
        zero(u, u.bounds());
}
