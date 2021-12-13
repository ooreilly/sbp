#pragma once
#include "transform/apply.hpp"
#include "functor/inner.hpp"

template <template <typename> class Array, typename Tv>
void innerf(Array<Tv>& out, 
                          const Array<Tv>& weights_x,
                          const Array<size_t>& indices_x,
                          const Array<Tv>& in) {
        Inner1D<Array, Tv> L(out, weights_x, indices_x, in);
        Bounds bounds({0, weights_x.nx});
        apply(L, bounds, out);
}

template <template <typename> class Array, typename Tv>
void innerf(Array<Tv>& out, 
                          const Array<Tv>& weights_x,
                          const Array<Tv>& weights_y,
                          const Array<size_t>& indices_x,
                          const Array<size_t>& indices_y,
                          const Array<Tv>& in) {
        Inner2D<Array, Tv> L(out, weights_x, weights_y, indices_x,
                                            indices_y, in);
        Bounds bounds({0, weights_x.nx});
        apply(L, bounds, out);
}

template <template <typename> class Array, typename Tv>
void innerf(Array<Tv>& out, 
                          const Array<Tv>& weights_x,
                          const Array<Tv>& weights_y,
                          const Array<Tv>& weights_z,
                          const Array<size_t>& indices_x,
                          const Array<size_t>& indices_y,
                          const Array<size_t>& indices_z, 
                          const Array<Tv>& in) {
        Inner3D<Array, Tv> L(out, weights_x, weights_y, weights_z, 
                                            indices_x,
                                            indices_y, indices_z, in);
        Bounds bounds({0, weights_x.nx});
        apply(L, bounds, out);
}

