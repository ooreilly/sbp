#pragma once
#include "macro/operator.hpp"

// v = c
template <template <typename> class Array, typename Tv>
class Constant {
        Tv c;
        Tv *v;
        size_t line, slice;

       public:
        Constant(Array<Tv>& v_, Tv c_ = 1.0) {
                v = v_();
                line = v_.line;
                slice = v_.slice;
                c = c_;
        }

        VOID_OPERATOR() const { 
                size_t pos = xi + line * yj + slice * zk;
                v[pos] = c;
        }
};

