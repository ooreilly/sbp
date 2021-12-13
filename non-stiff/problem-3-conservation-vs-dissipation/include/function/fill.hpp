#pragma once
#include "transform/apply.hpp"
#include "structure/bounds.hpp"
#include "functor/fill.hpp"
/* 
 * Take an array of values and put them into a host or device array
 */
template <template <typename> class Array, typename Tv=double>
void fill(Array<Tv>& out, const Tv* values, size_t nx, size_t ny = 1,
           size_t nz = 1) {
        Fill<Tv> f(out, values, nx, ny, nz);
        Bounds bounds({0, nx}, {0, ny}, {0, nz});
        apply(f, bounds, out);
}
