#pragma once
#include "structure/array.hpp"
#include "structure/bounds.hpp"

static const size_t reduceBlockDim_x = 2;
static const size_t reduceBlockDim_y = 2;
static const size_t reduceBlockDim_z = 2;

template <size_t blockDim_x, size_t blockDim_y, size_t blockDim_z,
          typename Tout, typename Tv, typename Functor>
Tout _host_reduce(const Tout start, const Functor& fctr, const Bounds& bounds,
                  const size_t blockIdx_x, const size_t blockIdx_y,
                  const size_t blockIdx_z, const Tv* u, const size_t unx,
                  const size_t uny, const size_t unz, const size_t uline,
                  const size_t uslice, const Tv* v, const size_t vnx,
                  const size_t vny, const size_t vnz, const size_t vline,
                  const size_t vslice) {
        size_t idx = reduceBlockDim_x * blockIdx_x + bounds.x[0];
        size_t idy = reduceBlockDim_y * blockIdx_y + bounds.y[0];
        size_t idz = reduceBlockDim_z * blockIdx_z + bounds.z[0];
        Tout out = start;
#pragma unroll
        for (size_t zk = 0; zk < reduceBlockDim_z; ++zk) {
#pragma unroll
                for (size_t yj = 0; yj < reduceBlockDim_y; ++yj) {
#pragma unroll
                        for (size_t xi = 0; xi < reduceBlockDim_x; ++xi) {
                                if (idx + xi < bounds.x[1] && idy + xi < vnx &&
                                    idy + yj < bounds.y[1] && idx + yj < vny &&
                                    idz + zk < bounds.z[1] && idz + zk < vnz) {
                                        size_t uxiyjzk = (idx + xi) +
                                                         uline * (idy + yj) +
                                                         (idz + zk) * uslice;
                                        size_t vxiyjzk = (idx + xi) +
                                                         vline * (idy + yj) +
                                                         (idz + zk) * vslice;
                                        out = fctr(out, idx + xi, idy + yj,
                                                   idz + zk, u, u[uxiyjzk], unx,
                                                   uny, unz, uline, uslice, v,
                                                   v[vxiyjzk], vnx, vny, vnz,
                                                   vline, vslice);
                                }
                        }
                }
        }
        return out;
}

template <typename Tout, typename Tv, typename Functor>
Tout reduce(const Functor& fctr, const Bounds& bounds, const HostArray<Tv>& u,
            const HostArray<Tv>& v, const Tout start = 0) {
        assert(bounds.x[1] <= u.nx);
        assert(bounds.y[1] <= u.ny);
        assert(bounds.z[1] <= u.nz);
        size_t gridDim_x = (bounds.x[1] - 1) / reduceBlockDim_x + 1;
        size_t gridDim_y = (bounds.y[1] - 1) / reduceBlockDim_y + 1;
        size_t gridDim_z = (bounds.z[1] - 1) / reduceBlockDim_z + 1;
        Tout out = start;
        for (size_t bz = 0; bz < gridDim_z; bz++) {
                for (size_t by = 0; by < gridDim_y; by++) {
                        for (size_t bx = 0; bx < gridDim_x; bx++) {
                                out = _host_reduce<
                                    reduceBlockDim_x, reduceBlockDim_y,
                                    reduceBlockDim_z, Tout, Tv, Functor>(
                                    out, fctr, bounds, bx, by, bz, u(), u.nx,
                                    u.ny, u.nz, u.line, u.slice, v(), v.nx,
                                    v.ny, v.nz, v.line, v.slice);
                        }
                }
        }
        return out;
}

