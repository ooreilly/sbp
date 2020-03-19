#pragma once
#include "transform/apply.hpp"
#include "structure/array.hpp"
#include "functor/multiply.hpp"
#include "macro/dim.hpp"


// w = a * u .* v + b * w + c
template <template <typename> class Array, typename To=double, typename Ti=double, typename Tc=double>
void multiply(Array<To>& w, const Array<Ti>& u, const Array<Ti>& v, const Bounds& bounds, Tc a = 1.0,
              Tc b = 0.0, Tc c = 0.0) {
        Multiply<Array, To, Ti, Tc> mul(w, u, v, a, b, c);
        return apply(mul, bounds, w);
}

template <template <typename> class Array, typename To=double, typename Ti=double, typename Tc=double>
void multiply(Array<To>& w, const Array<Ti>& u, const Array<Ti>& v, Tc a = 1.0, Tc b = 0.0,
              Tc c = 0.0) {
        return multiply(w, u, v, w.bounds(), a, b, c);
}
