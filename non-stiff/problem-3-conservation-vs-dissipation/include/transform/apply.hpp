#pragma once
#include "function/cmp_bounds_size.hpp"
#include "function/struct.hpp"
#include "structure/array.hpp"
#include "structure/bounds.hpp"

/*----------------------------------------------------------------------------*/
// Device

template <typename GridFcn>
void device_apply(GridFcn& h_fcn, const Bounds& bounds) {
        int numSMs;
        const int numBlocks_y = 2;
        dim3 threads(16, 16, 1);
        GridFcn* d_fcn;
        dstruct(&d_fcn, h_fcn);
        cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, 0);
        dim3 blocks(numSMs, numBlocks_y, 1);

        _device_apply<<<blocks, threads>>>(
            bounds.x[0], bounds.y[0], bounds.z[0],
            bounds.x[1], bounds.y[1], bounds.z[1],
            *d_fcn);
        dfree(d_fcn);
}

template <typename GridFcn>
__global__ void _device_apply(const size_t bounds_x0, const size_t bounds_y0,
                              const size_t bounds_z0, const size_t bounds_x1,
                              const size_t bounds_y1, const size_t bounds_z1,
                              GridFcn& fcn) {
        int idx = threadIdx.x + blockDim.x * blockIdx.x + bounds_x0;
        int idy = threadIdx.y + blockDim.y * blockIdx.y + bounds_y0;
        int idz = threadIdx.z + blockDim.z * blockIdx.z + bounds_z0;

        const int xn = bounds_x1;
        const int yn = bounds_y1;
        const int zn = bounds_z1;

        for (int zk = idz; zk < zn; zk += blockDim.z * gridDim.z) {
                for (int yj = idy; yj < yn; yj += blockDim.y * gridDim.y) {
                        for (int xi = idx; xi < xn;
                             xi += blockDim.x * gridDim.x) {
                                fcn(xi, yj, zk);
                        }
                }
        }
}

template <typename Tv, typename GridFcn>
__inline__ void apply(GridFcn& h_fcn, const Bounds& bounds, const DeviceArray<Tv>& dummy) {
        device_apply(h_fcn, bounds);
}

/*----------------------------------------------------------------------------*/
// Host
const size_t applyBlockDim_x = 2;
const size_t applyBlockDim_y = 2;
const size_t applyBlockDim_z = 2;

template <size_t blockDim_x, size_t blockDim_y, size_t blockDim_z,
          typename Functor>
void _host_apply(const Bounds& bounds,
                 const size_t blockIdx_x, const size_t blockIdx_y,
                 const size_t blockIdx_z, const Functor& fctr
                 ) {
        size_t idx = applyBlockDim_x * blockIdx_x;
        size_t idy = applyBlockDim_y * blockIdx_y;
        size_t idz = applyBlockDim_z * blockIdx_z;

        int xn = bounds.x[1] - bounds.x[0];
        int yn = bounds.y[1] - bounds.y[0];
        int zn = bounds.z[1] - bounds.z[0];

#pragma unroll
        for (size_t zk = 0; zk < applyBlockDim_z; ++zk) {
#pragma unroll
                for (size_t yj = 0; yj < applyBlockDim_y; ++yj) {
#pragma unroll
                        for (size_t xi = 0; xi < applyBlockDim_x; ++xi) {
                                if (idx + xi < xn && idy + yj < yn &&
                                    idz + zk < zn) {
                                            fctr(idx + xi + bounds.x[0],
                                                 idy + yj + bounds.y[0],
                                                 idz + zk + bounds.z[0]);
                                }
                        }
                }
        }
}

template <typename Functor>
void host_apply(const Functor& fctr, 
                const Bounds& bounds) { 
        size_t gridDim_x = (bounds.x[1] - 1) / applyBlockDim_x + 1;
        size_t gridDim_y = (bounds.y[1] - 1) / applyBlockDim_y + 1;
        size_t gridDim_z = (bounds.z[1] - 1) / applyBlockDim_z + 1;
        for (size_t bz = 0; bz < gridDim_z; bz++) {
                for (size_t by = 0; by < gridDim_y; by++) {
                        for (size_t bx = 0; bx < gridDim_x; bx++) {
                                _host_apply<applyBlockDim_x, applyBlockDim_y,
                                            applyBlockDim_z, Functor>(
                                    bounds, bx, by, bz, fctr);
                        }
                }
        }
}

template <typename Tv, typename GridFcn>
__inline__ void apply(GridFcn& h_fcn, const Bounds& bounds, const HostArray<Tv>& dummy) {
        host_apply(h_fcn, bounds);
}
