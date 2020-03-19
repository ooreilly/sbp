#pragma once

template <typename T>
__inline__ void copy_initializer_list(T* out,
                                      std::initializer_list<T>& in) {
        int i = 0;
        for (auto xi : in) {
                out[i] = xi;
                i++;
        }
}
