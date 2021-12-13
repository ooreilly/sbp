#pragma once
#include "macro/operator.hpp"

template <typename Tv>
class Identity {

       public:
        Identity() {
        }

        OPERATOR() const {
                return in_val1;
        }
};
