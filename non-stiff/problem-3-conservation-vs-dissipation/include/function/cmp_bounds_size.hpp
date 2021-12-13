#pragma once
#include "structure/bounds.hpp"

// Return true if two bounds objects are of the same size.
bool cmp_bounds_size(const Bounds& b1, const Bounds& b2) {
        bool stat = true;
        stat &= b1.x[1] - b1.x[0] == b2.x[1] - b2.x[0];
        stat &= b1.y[1] - b1.y[0] == b2.y[1] - b2.y[0];
        stat &= b1.z[1] - b1.z[0] == b2.z[1] - b2.z[0];
        return stat;
}
