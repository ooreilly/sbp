#pragma once
#include "macro/operator.hpp"

// f(u) = 0
template <typename Tv>
class Zero {
       public:
        Zero() {}

        OPERATOR() const {
                return 0.0;
        }
};

