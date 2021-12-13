#pragma once
#include "function/lagrange.hpp"
#include "function/constant.hpp"
#include "function/inner.hpp"
#include "function/outer.hpp"
#include "macro/error.hpp"

template <template <typename> class Array, typename Tv>
class Lagrange {

        public:
                Array<size_t> indices_x, indices_y, indices_z;
                Array<Tv> weights_x, weights_y, weights_z;
                int dim = 1;
                int num_values;

         Lagrange() {
         }

         Lagrange(const Array<Tv>& x, const Array<Tv>& grid_x, 
                  int degree) {
                num_values = x.nx;
                 init(weights_x, indices_x, degree);
                 set_weights(weights_x, indices_x, x, grid_x, degree);

                 dim = 1;
                }

         void init(Array<Tv>& x, Array<size_t>& indices, int degree) {
                        x = Array<Tv>(num_values, degree + 1);
                        indices = Array<size_t>(num_values, degree + 1);
         }

         Lagrange(const Array<Tv>& x, const Array<Tv>& y,
                  const Array<Tv>& grid_x, const Array<Tv>& grid_y,
                  int degree) {
                 num_values = x.nx;
                 assert(x.nx == y.nx);
                 init(weights_x, indices_x, degree);
                 init(weights_y, indices_y, degree);
                 set_weights(weights_x, indices_x, x, grid_x, degree);
                 set_weights(weights_y, indices_y, y, grid_y, degree);
                 dim = 2;
                }

         Lagrange(const Array<Tv>& x, 
                  const Array<Tv>& y, 
                  const Array<Tv>& z,
                  const Array<Tv>& grid_x,        
                  const Array<Tv>& grid_y,        
                  const Array<Tv>& grid_z,        
                  int degree) {
                 num_values = x.nx;
                 assert(x.nx == y.nx);
                 assert(x.nx == z.nx);

                 init(weights_x, indices_x, degree);
                 init(weights_y, indices_y, degree);
                 init(weights_z, indices_z, degree);

                 set_weights(weights_x, indices_x, x, grid_x, degree);
                 set_weights(weights_y, indices_y, y, grid_y, degree);
                 set_weights(weights_z, indices_z, z, grid_z, degree);
                 dim = 3;
                }

         void interpolate(Array<Tv>& out, const Array<Tv>& in)const{
                 assert(out.nx ==  num_values);
                 switch (dim) {
                         case 1:
                                 innerf(out, weights_x, indices_x, in);
                                 break;
                         case 2:
                                 innerf(out, weights_x, weights_y,
                                               indices_x, indices_y, in);
                                 break;
                         case 3:
                                 innerf(out, weights_x, weights_y,
                                               weights_z, indices_x, indices_y,
                                               indices_z, in);
                                 break;
                }
         }

         void outer(Array<Tv>& out, const Array<Tv>& a,
                    Tv b = 1.0) const {
                 assert(a.nx ==  num_values);
                 switch (dim) {
                         case 1:
                         outerf(out, weights_x, indices_x, a, b);
                         break;
                         case 2:
                         outerf(out, weights_x, weights_y, indices_x,
                                        indices_y, a, b);
                         break;
                         case 3:
                         outerf(out, weights_x, weights_y, weights_z,
                                        indices_x, indices_y, indices_z, a, b);
                         break;
                }
         }

        private:
         void set_weights(Array<Tv>& weights, Array<size_t>& indices,
                          const Array<Tv>& query, const Array<Tv>& grid,
                          int degree) const{
                 argnearest(indices, query, grid);
                 lagrange_weights(weights, query, indices, grid);
         }
};
