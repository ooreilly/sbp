#pragma once
#include "functor/lineop.hpp"
#include "structure/bounds.hpp"

/* Differentiate grid function using fourth order staggered SBP finite
 difference operator (v lives on (+)-grid) and u lives on (-)-grid) */
template <typename Tv=double>
class Dp42 {
        public:
                const int m_left = 4;
                const char *label = "Dp4";
                const int n_left = 6;
                const int m_right = 5;
                const int n_right = 6;
                const int n_interior = 4;
                const int offset_interior = -1;
                const int order_left = 2;
                const int order_interior = 4;
                const int order_right = 2;

        Tv hi;

        const Tv left[4][6] = {
            {-1.7779989465546748, 1.3337480247900155, 0.7775013168066564,
             -0.33325039504199694, 0.0, 0.0},
            {-0.44102173413920587, -0.17308424848898904, 0.4487228323259927,
             0.1653831503022022, 0.0, 0.0},
            {0.17987932138827012, -0.27572572541507884, -0.9597948548284453,
             1.1171892610431817, -0.06154800218792766, 0.0},
            {0.015391138150708755, 0.05688514555035906, -0.19989764645971708,
             -0.8628231468598346, 1.0285385292191949, -0.03809401960071092}};

        const Tv right[5][6] = {
            {0.03809401960071092, -1.0285385292191949, 0.8628231468598346,
             0.19989764645971708, -0.05688514555035906, -0.015391138150708755},
            {-0.0, 0.06154800218792766, -1.1171892610431817, 0.9597948548284453,
             0.27572572541507884, -0.17987932138827012},
            {-0.0, -0.0, -0.1653831503022022, -0.4487228323259927,
             0.17308424848898904, 0.44102173413920587},
            {-0.0, -0.0, 0.33325039504199694, -0.7775013168066564,
             -1.3337480247900155, 1.7779989465546748},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

        const Tv interior[4] = {0.041666666666666664, -1.12500000000000000000,
                                1.1250000000000000, -0.041666666666666664};
        Dp42() { }
        Dp42 operator=(Dp42& in) {return Dp42(); }
};

template <int dim, template <typename> class Array, typename To = double,
          typename Ti = double, typename Tc = double>
void dp42(Array<To>& out, const Array<Ti>& in, const Bounds& bounds, Tc a = 1.0,
         Tc b = 0.0, Tc c = 0.0) {
        const Dp42<Ti> opdata;
        LineOp<dim, Array, 6, 6, 4, Dp42<To>, To, Ti, Tc> D(opdata, out, in, a,
                                                             b, c);
        apply(D, bounds, out);
}

// v = a * Dx * u  + b * v + c
template <template <typename> class Array, typename To = double,
          typename Ti = double, typename Tc = double>
void dpx42(Array<To>& out, const Array<Ti>& in, const Bounds& bounds, Tc a = 1.0,
          Tc b = 0.0, Tc c = 0.0) {
        dp42<X>(out, in, bounds, a, b, c);
}

template <template <typename> class Array, typename To = double,
          typename Ti = double, typename Tc = double>
void dpx42(Array<To>& out, const Array<Ti>& in, Tc a = 1.0, Tc b = 0.0,
          Tc c = 0.0) {
        dp42<X>(out, in, out.bounds(), a, b, c);
}

// v = a * Dy * u  + b * v + c
template <template <typename> class Array, typename To = double,
          typename Ti = double, typename Tc = double>
void dpy42(Array<To>& out, const Array<Ti>& in, const Bounds& bounds, Tc a = 1.0,
          Tc b = 0.0, Tc c = 0.0) {
        dp42<Y>(out, in, bounds, a, b, c);
}

template <template <typename> class Array, typename To = double,
          typename Ti = double, typename Tc = double>
void dpy42(Array<To>& out, const Array<Ti>& in, Tc a = 1.0, Tc b = 0.0,
          Tc c = 0.0) {
        dp42<Y>(out, in, out.bounds(), a, b, c);
}

// v = a * Dz * u  + b * v + c
template <template <typename> class Array, typename To = double,
          typename Ti = double, typename Tc = double>
void dpz42(Array<To>& out, const Array<Ti>& in, const Bounds& bounds, Tc a = 1.0,
          Tc b = 0.0, Tc c = 0.0) {
        dp42<Z>(out, in, bounds, a, b, c);
}

template <template <typename> class Array, typename To = double,
          typename Ti = double, typename Tc = double>
void dpz42(Array<To>& out, const Array<Ti>& in, Tc a = 1.0, Tc b = 0.0,
          Tc c = 0.0) {
        dp42<Z>(out, in, out.bounds(), a, b, c);
}
