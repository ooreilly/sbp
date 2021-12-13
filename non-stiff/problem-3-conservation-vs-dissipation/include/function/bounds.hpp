#pragma once
#include "structure/bounds.hpp"
#include "macro/dim.hpp"
#include "function/getdim.hpp"
#include "macro/sides.hpp"

Bounds bounds(size_t nx=1, size_t ny=1, size_t nz=1) {
        return Bounds({0, nx}, {0, ny}, {0, nz});
}
Bounds fbounds(size_t nx=1, size_t ny=1, size_t nz=1) {
        return bounds(nx, ny, nz);
}

Bounds* bounds_sides1d(size_t nx, size_t ax=1) {
        static Bounds out[2];
        out[LEFT] = Bounds({0, ax});
        out[RIGHT] = Bounds({nx - ax, nx});
        return out;
}


void bounds_sides2d(Bounds *out, size_t nx, size_t ny, size_t ax=1, size_t ay=1) {
                out[LEFT] = Bounds({0, ax}, {0, ny});
                out[RIGHT] = Bounds({nx-ax, nx}, {0, ny});
                out[BOTTOM] = Bounds({0, nx}, {0, ay});
                out[TOP] = Bounds({0, nx}, {ny - ay, ny});
}

Bounds bounds_point(size_t ix=0, size_t iy=0, size_t iz=0) {
        return Bounds({ix, ix + 1}, {iy, iy + 1}, {iz, iz + 1});
}

// Bounds for operator 
template <int dim>
Bounds leftbounds(size_t left, size_t nx=1, size_t ny=1, size_t nz=1) {
        size_t n = getdim<dim>(nx, ny, nz);
        Bounds bx({0, left}, {0, ny}, {0, nz});
        Bounds by({0, nx}, {0, left}, {0, nz});
        Bounds bz({0, nx}, {0, ny}, {0, left});
        return getdim<dim>(bx, by, bz);
}

template <int dim>
Bounds rightbounds(size_t right, size_t nx=1, size_t ny=1, size_t nz=1) {
        size_t n = getdim<dim>(nx, ny, nz);
        std::initializer_list<size_t> bright = {n - right, n};
        Bounds bx(bright, {0, ny}, {0, nz});
        Bounds by({0, nx}, bright, {0, nz});
        Bounds bz({0, nx}, {0, ny}, bright);
        return getdim<dim>(bx, by, bz);
}

template <int dim>
Bounds intbounds(size_t left, size_t right, size_t nx=1, size_t ny=1, size_t nz=1) {
        size_t n = getdim<dim>(nx, ny, nz);
        std::initializer_list<size_t> bint = {left, n - right};
        Bounds bx(bint, {0, ny}, {0, nz});
        Bounds by({0, nx}, bint, {0, nz});
        Bounds bz({0, nx}, {0, ny}, bint);
        return getdim<dim>(bx, by, bz);
}
