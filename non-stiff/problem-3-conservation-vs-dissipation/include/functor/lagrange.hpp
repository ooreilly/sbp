#pragma once

#include "function/argnearest.hpp"

template <template <typename> class Array, typename Tv>
class LagrangeWeights
{
        Tv *weights;
        size_t weights_line;
        const Tv *queries, *grid;
        const size_t *indices;
        size_t num_query, num_basis, indices_line;
        public:
         LagrangeWeights(
                         Array<Tv>& weights_, 
                         const Array<Tv>& queries_, 
                         const Array<size_t>& indices_, 
                         const Array<Tv>& grid_
                         ) {
                 queries = queries_();
                 indices = indices_();
                 grid = grid_();
                 num_query = indices_.nx;
                 num_basis = indices_.ny;
                 indices_line = indices_.line;

                 weights = weights_();
                 weights_line = weights_.line;
        }


        VOID_OPERATOR() const {
                // Compute weight w_j for each query xi, i.e.,
                // xi : query coordinate
                // yj : basis function

                double N = 1.0;
                double denom = 1.0;
                double xa = (double)queries[xi];
                size_t idx_j = indices[xi + indices_line * yj];
                double xj = (double)grid[idx_j];
                // Determine node polynomial
                for (int k = 0; k < num_basis; k++) {
                        size_t idx_k = indices[xi + indices_line * k];
                        double xk = (double)grid[idx_k];
                        N = N * (xa - xk);
                        if (idx_j != idx_k) 
                                denom *= xj - xk;
                }

                Tv weight = xa == xj ? 1 : (Tv) (N / denom / (xa - xj)); 

                size_t pos = xi + weights_line * yj;
                weights[pos] =  weight;
        }

};

