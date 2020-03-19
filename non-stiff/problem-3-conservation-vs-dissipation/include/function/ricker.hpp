#pragma once

#ifndef M_PI
#define M_PI 3.14159265358979
#endif

// Ricker wavelet
// A * (1 - 2 * s) * exp(-s), s = pi^2 * fp^2 * (t - t0)^2  
template <typename Tv>
Tv ricker(Tv t, Tv t0, Tv fp, Tv A=1.0) {
        Tv s = M_PI * M_PI * fp * fp * pow(t - t0, 2);
        return A * (1 - 2 * s) * exp(-s);
}
