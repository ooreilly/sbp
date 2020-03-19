#pragma once
#include "function/getdim.hpp"

template <int dim, typename Tvo, typename Tvi, int n_left, int n_right, int n_interior, typename Spec>
__inline__ __host__ __device__ Tvo lineop(const Tvi *u, size_t xi, size_t yj, size_t zk, 
                                          const Spec& spec, 
                                          size_t nx, size_t ny,
                                          size_t nz, size_t line,
                                          size_t slice) {
        Tvo out = 0.0;

        size_t dir = getdim<dim>(xi, yj, zk);
        size_t ndim = getdim<dim>(nx, ny, nz);

        size_t inc = getdim<dim>((size_t)1, line, slice);
        size_t pos = xi + line * yj + slice * zk;

        size_t bnd_pos = getdim<dim>(line * yj + slice * zk,
                                     xi + slice * zk, 
                                     xi + line * yj);



        if (dir < spec.m_left) {
#pragma unroll
                for (int k = 0; k < n_left; ++k) {
                        out += spec.left[dir][k] * u[bnd_pos + inc * k];
                }
        } else if (dir >= ndim - spec.m_right) {
                size_t m_shift = ndim - spec.m_right;
                size_t n_shift = ndim - spec.n_right;
#pragma unroll
                for (int k = 0; k < n_right; ++k) {

                        out += spec.right[(dir - m_shift)][k] *
                               u[bnd_pos + inc * (n_shift + k)];
                                }
        } else {
#pragma unroll
                for (int k = 0; k < n_interior; ++k) {
                        out +=
                            spec.interior[k] * u[pos + inc * (k + spec.offset_interior)];
                }
        }

        return out;
}
