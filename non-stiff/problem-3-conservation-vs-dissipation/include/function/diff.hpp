#pragma once
#include "macro/dim.hpp"

template <typename Tv, int m_left, int n_left, int m_right, int n_right,
          int offset_interior, int n_interior>
SCALAR_OPERATOR(diff_x) {
        Tv Dx_u = 0.0;

        if (xi < m_left) {
#pragma unroll
                for (int k = 0; k < n_left; ++k) {
                        Dx_u += left[k + n_left * xi] *
                                ref1[k + line * yj + slice * zk];
                }
        } else if (xi >= nx - m_right) {
                int m_shift = nx - m_right;
                int n_shift = nx - n_right;
#pragma unroll
                for (int k = 0; k < n_right; ++k)
                        Dx_u += right[k + n_right * (xi - m_shift)] *
                                ref1[n_shift + k + line * yj + slice * zk];
        } else {
#pragma unroll
                for (int k = 0; k < n_interior; ++k) {
                        Dx_u +=
                            interior[k] * ref1[(int)xi + k + offset_interior +
                                              line * yj + slice * zk];
                }
        }

        return Dx_u;
}

template <typename Tv, int m_left, int n_left, int m_right, int n_right,
          int offset_interior, int n_interior>
SCALAR_OPERATOR(diff_y) {
        Tv Dy_u = 0.0;

        if (yj < m_left) {
#pragma unroll
                for (int k = 0; k < n_left; ++k) {
                        Dy_u += left[k + n_left * yj] *
                                ref1[xi + line * k + slice * zk];
                }
        } else if (yj >= ny - m_right) {
                int m_shift = ny - m_right;
                int n_shift = ny - n_right;
#pragma unroll
                for (int k = 0; k < n_right; ++k)
                        Dy_u += right[k + n_right * (yj - m_shift)] *
                                ref1[xi + line * (n_shift + k) + slice * zk];
        } else {
#pragma unroll
                for (int k = 0; k < n_interior; ++k) {
                        Dy_u +=
                            interior[k] *
                            ref1[xi + line * ((int)yj + k + offset_interior) +
                                slice * zk];
                }
        }

        return Dy_u;
}

template <typename Tv, int m_left, int n_left, int m_right, int n_right,
          int offset_interior, int n_interior>
SCALAR_OPERATOR(diff_z) {
        Tv Dz_u = 0.0;

        if (zk < m_left) {
#pragma unroll
                for (int k = 0; k < n_left; ++k) {
                        Dz_u += left[k + n_left * zk] *
                                ref1[xi + line * yj + slice * k];
                }
        } else if (zk >= nz - m_right) {
                int m_shift = nz - m_right;
                int n_shift = nz - n_right;
#pragma unroll
                for (int k = 0; k < n_right; ++k)
                        Dz_u += right[k + n_right * (zk - m_shift)] *
                                ref1[xi + line * yj + slice * (n_shift + k)];
        } else {
#pragma unroll
                for (int k = 0; k < n_interior; ++k) {
                        Dz_u += interior[k] *
                                ref1[xi + line * yj +
                                    slice * ((int)zk + k + offset_interior)];
                }
        }

        return Dz_u;
}

template <int dim, typename Tv, int m_left, int n_left, int m_right,
          int n_right, int offset_interior, int n_interior>
SCALAR_OPERATOR(diff) {
        if (dim == X)
                return diff_x<Tv, m_left, n_left, m_right, n_right,
                              offset_interior, n_interior>(
                    xi, yj, zk, ref1, ref2, out_val, in_val1, in_val2, nx, ny, nz,
                    line, slice, inplace, left, right, interior);
        if (dim == Y)
                return diff_y<Tv, m_left, n_left, m_right, n_right,
                              offset_interior, n_interior>(
                    xi, yj, zk, ref1, ref2, out_val, in_val1, in_val2, nx, ny, nz,
                    line, slice, inplace, left, right, interior);
        if (dim == Z)
                return diff_z<Tv, m_left, n_left, m_right, n_right,
                              offset_interior, n_interior>(
                    xi, yj, zk, ref1, ref2, out_val, in_val1, in_val2, nx, ny, nz,
                    line, slice, inplace, left, right, interior);
}
