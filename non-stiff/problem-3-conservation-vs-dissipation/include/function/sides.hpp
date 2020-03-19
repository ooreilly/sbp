#pragma once
#include "macro/sides.hpp"

template <typename Array>
void sides_array(Array *out, int num_sides, size_t nx = 1, size_t ny = 1,
           size_t nz = 1, size_t gx = 1, size_t gy = 1, size_t gz = 1) {
        out[LEFT] = Array(gx, ny, nz);
        out[RIGHT] = Array(gx, ny, nz);
        if (num_sides == 2) return;
        out[BOTTOM] = Array(nx, gy, nz);
        out[TOP] = Array(nx, gy, nz);
        if (num_sides == 4) return;
        out[FRONT] = Array(nx, ny, gz);
        out[BACK] = Array(nx, ny, gz);
}

void sides_bounds(Bounds *out, int num_sides, size_t nx = 1, size_t ny = 1,
                  size_t nz = 1, size_t ax = 1, size_t ay = 1, size_t az = 1) {
        out[LEFT] = Bounds({0, ax}, {0, ny});
        out[RIGHT] = Bounds({nx - ax, nx}, {0, ny});
        if (num_sides == 2) return;
        out[BOTTOM] = Bounds({0, nx}, {0, ay});
        out[TOP] = Bounds({0, nx}, {ny - ay, ny});
        if (num_sides == 4) return;
        out[FRONT] = Bounds({0, nx}, {0, ny}, {nz - az, nz});
        out[BACK] = Bounds({0, nx}, {0, ny}, {nz - ay, nz});
}
