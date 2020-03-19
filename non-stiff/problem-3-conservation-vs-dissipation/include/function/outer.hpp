#pragma once
#include "functor/outer.hpp"
#include "transform/apply.hpp"

// v(i,j,k) = w_i * a_i + b * v(i,j,k) for i in Lagrange neighborhood
template <template <typename> class Array, typename Tv>
void outerf(Array<Tv>& out, 
                          const Array<Tv>& weights_x,
                          const Array<size_t>& indices_x,
                          const Array<Tv>& a,
                          Tv b = 1.0
                          ) {
        Outer1D<Array, Tv> L(out, weights_x, indices_x, a, b);
        apply(L, out.bounds(), out);
}

template <template <typename> class Array, typename Tv>
void outerf(Array<Tv>& out, 
                          const Array<Tv>& weights_x,
                          const Array<Tv>& weights_y,
                          const Array<size_t>& indices_x,
                          const Array<size_t>& indices_y,
                          const Array<Tv>& a,
                          Tv b = 1.0
                          ) {
        Outer2D<Array, Tv> L(out, weights_x, weights_y, indices_x, indices_y,
                                     a, b);
        apply(L, out.bounds(), out);
}

template <template <typename> class Array, typename Tv>
void outerf(Array<Tv>& out, 
                          const Array<Tv>& weights_x,
                          const Array<Tv>& weights_y,
                          const Array<Tv>& weights_z,
                          const Array<size_t>& indices_x,
                          const Array<size_t>& indices_y,
                          const Array<size_t>& indices_z,
                          const Array<Tv>& a,
                          Tv b = 1.0
                          ) {
        Outer3D<Array, Tv> L(out, weights_x, weights_y, weights_z, 
                                     indices_x, indices_y, indices_z,
                                     a, b);
        apply(L, out.bounds(), out);
}

