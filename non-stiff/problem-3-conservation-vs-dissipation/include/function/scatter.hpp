#pragma once

#include "structure/bounds.hpp"
#include "functor/scatter.hpp"
#include "function/cmp_bounds_size.hpp"
#include "transform/apply.hpp"


// v[bounds] = a * u + b * v[bounds] + c
template <template <typename> class Array, typename Tv = double,
          typename Tc = double>
void scatter(Array<Tv>& x_scatter, const Array<Tv>& x,
             const Bounds& bounds_scatter, const Tc a = 1.0, const Tc b = 0.0,
             const Tc c = 0.0) {
        assert(cmp_bounds_size(x.bounds(), bounds_scatter));
        Scatter<Array, Tv, Tc> s(x_scatter, x, bounds_scatter, a, b, c);
        apply(s, x.bounds(), x);

}
