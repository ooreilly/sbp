#pragma once
#include "functor/lineop.hpp"
#include "structure/bounds.hpp"

/* Differentiate grid function using classical second order SBP finite
 difference operator */
template <typename Tv=double>
class SBP21 {
        public:
        const size_t m_left = 1;
        const size_t n_left = 2;
        const size_t m_right = 1;
        const size_t n_right = 2;
        const size_t n_interior = 3;
        const int offset_interior = -1;
        const int order_left = 1;
        const int order_right = 1;
        const int order_interior = 1;
        const char *label = "D21";

        const Tv left[1][2] = {{-1, 1}};
        const Tv right[1][2] = {{-1, 1}};
        const Tv interior[3] = {-0.5, 0.0, 0.5};
        SBP21() { }
        SBP21 operator=(SBP21& in) {return SBP21(); }
};

template <int dim, template <typename> class Array, typename To = double,
          typename Ti = double, typename Tc = double>
void d21(Array<To>& out, const Array<Ti>& in, const Bounds& bounds, Tc a = 1.0,
         Tc b = 0.0, Tc c = 0.0) {
        const SBP21<Ti> opdata;
        LineOp<dim, Array, 2, 2, 3, SBP21<To>, To, Ti, Tc> D(opdata, out, in, a,
                                                             b, c);
        apply(D, bounds, out);
}

// v = a * Dx * u  + b * v + c
template <template <typename> class Array, typename To = double,
          typename Ti = double, typename Tc = double>
void dx21(Array<To>& out, const Array<Ti>& in, const Bounds& bounds, Tc a = 1.0,
          Tc b = 0.0, Tc c = 0.0) {
        d21<X>(out, in, bounds, a, b, c);
}

template <template <typename> class Array, typename To = double,
          typename Ti = double, typename Tc = double>
void dx21(Array<To>& out, const Array<Ti>& in, Tc a = 1.0, Tc b = 0.0,
          Tc c = 0.0) {
        d21<X>(out, in, out.bounds(), a, b, c);
}

// v = a * Dy * u  + b * v + c
template <template <typename> class Array, typename To = double,
          typename Ti = double, typename Tc = double>
void dy21(Array<To>& out, const Array<Ti>& in, const Bounds& bounds, Tc a = 1.0,
          Tc b = 0.0, Tc c = 0.0) {
        d21<Y>(out, in, bounds, a, b, c);
}

template <template <typename> class Array, typename To = double,
          typename Ti = double, typename Tc = double>
void dy21(Array<To>& out, const Array<Ti>& in, Tc a = 1.0, Tc b = 0.0,
          Tc c = 0.0) {
        d21<Y>(out, in, out.bounds(), a, b, c);
}

// v = a * Dz * u  + b * v + c
template <template <typename> class Array, typename To = double,
          typename Ti = double, typename Tc = double>
void dz21(Array<To>& out, const Array<Ti>& in, const Bounds& bounds, Tc a = 1.0,
          Tc b = 0.0, Tc c = 0.0) {
        d21<Z>(out, in, bounds, a, b, c);
}

template <template <typename> class Array, typename To = double,
          typename Ti = double, typename Tc = double>
void dz21(Array<To>& out, const Array<Ti>& in, Tc a = 1.0, Tc b = 0.0,
          Tc c = 0.0) {
        d21<Z>(out, in, out.bounds(), a, b, c);
}
