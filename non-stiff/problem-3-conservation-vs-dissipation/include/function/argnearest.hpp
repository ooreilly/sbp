#pragma once
#include "transform/apply.hpp"
#include "structure/array.hpp"
#include "functor/argnearest.hpp"
#include "functor/identity.hpp"

/*
 Argnearest finds the k-nearest indices of a query point within a list of
 points.

        Args:
         out(num_query, num_nearest) : Output array to assign to
         query(num_query, dimension) : Input array of queries
         points(num_points, dimension) : Input array of values to search through
*/
template <template <typename> class Array, typename Tvo, typename Tvi>
void argnearest(Array<Tvo>& out_arg, const Array<Tvi>& query,
                const Array<Tvi>& points) {
        ArgNearest<Array, Tvo, Tvi> arg(out_arg, query, points);
        assert(out_arg.nx == query.nx);
        // Check that query and grid have the same dimensions
        assert(query.ny == points.ny);
        assert(query.ny <= 3);
        assert(points.ny <= 3);
        
        assert(query.nz == 1);
        assert(points.nz == 1);

        apply(arg, out_arg.bounds(), query);
}
