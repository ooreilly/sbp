#pragma once
#include "transform/apply.hpp"
#include "functor/lagrange.hpp"

template <template <typename> class Array, typename Tv>
void lagrange_weights(Array<Tv>& weights, const Array<Tv>& query,
                      const Array<size_t>& indices, const Array<Tv>& grid) {

        LagrangeWeights<Array, Tv> L(weights, query, indices, grid);
        apply(L, weights.bounds(), weights);
}
