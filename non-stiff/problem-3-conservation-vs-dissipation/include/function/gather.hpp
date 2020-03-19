#pragma once

#include "structure/bounds.hpp"
#include "functor/gather.hpp"
#include "function/cmp_bounds_size.hpp"
#include "transform/apply.hpp"

// v = a * u[bounds] + b
template <template <typename> class Array, typename Tv = double,
          typename Tc = double>
void gather(Array<Tv>& x_gather, const Array<Tv>& x, const Bounds& bounds_x,
             const Tc a = 1.0, const Tc b = 0.0, const Tc c = 0.0) {
        assert(cmp_bounds_size(x_gather.bounds(), bounds_x));
        Gather<Array, Tv, Tc> g(x_gather, x, bounds_x, a, b, c);
        apply(g, x_gather.bounds(), x);
}
