#pragma once
#include "functor/lineop.hpp"
#include "structure/bounds.hpp"

/* Differentiate grid function using fourth order staggered SBP finite
 difference operator (v lives on (-)-grid) and u lives on (+)-grid) */
template <typename Tv=double>
class Dm42 {
        public:
        const int m_left = 4;
        const int n_left = 5;
        const int m_right = 4;
        const int n_right = 6;
        const int n_interior = 4;
        const int offset_interior = -2;
        const int order_left = 2;
        const int order_interior = 4;
        const int order_right = 2;
        const char *label = "Dm42";
        
        const Tv left[4][5] = {
            {-1.4511412472637157, 1.853423741791147, -0.3534237417911469,
             -0.04885875273628437, 0.0},
            {-0.8577143189081458, 0.5731429567244373, 0.4268570432755628,
             -0.14228568109185424, 0.0},
            {-0.16745485058828774, -0.4976354482351368, 0.4976354482351368,
             0.16745485058828774, 0.0},
            {0.10270611134051243, -0.26245413264698597, -0.8288742701021167,
             1.0342864927831414, -0.04566420137455131}};

        const Tv right[4][6] = {
            {0.04566420137455131, -1.0342864927831414, 0.8288742701021167,
             0.26245413264698597, -0.10270611134051243, 0.0},
            {-0.0, -0.16745485058828774, -0.4976354482351368,
             0.4976354482351368, 0.16745485058828774, 0.0},
            {-0.0, 0.14228568109185424, -0.4268570432755628,
             -0.5731429567244373, 0.8577143189081458, 0.0},
            {-0.0, 0.04885875273628437, 0.3534237417911469, -1.853423741791147,
             1.4511412472637157, 0.0}};

        const Tv interior[4] = {0.041666666666666664, -1.12500000000000000000,
                                1.1250000000000000, -0.041666666666666664};
        Dm42() { }
        Dm42 operator=(Dm42& in) {return Dm42(); }
};

template <int dim, template <typename> class Array, typename To = double,
          typename Ti = double, typename Tc = double>
void dm42(Array<To>& out, const Array<Ti>& in, const Bounds& bounds, Tc a = 1.0,
         Tc b = 0.0, Tc c = 0.0) {
        const Dm42<Ti> opdata;
        LineOp<dim, Array, 5, 6, 4, Dm42<To>, To, Ti, Tc> D(opdata, out, in, a,
                                                             b, c);
        apply(D, bounds, out);
}

// v = a * Dx * u  + b * v + c
template <template <typename> class Array, typename To = double,
          typename Ti = double, typename Tc = double>
void dmx42(Array<To>& out, const Array<Ti>& in, const Bounds& bounds, Tc a = 1.0,
          Tc b = 0.0, Tc c = 0.0) {
        dm42<X>(out, in, bounds, a, b, c);
}

template <template <typename> class Array, typename To = double,
          typename Ti = double, typename Tc = double>
void dmx42(Array<To>& out, const Array<Ti>& in, Tc a = 1.0, Tc b = 0.0,
          Tc c = 0.0) {
        dm42<X>(out, in, out.bounds(), a, b, c);
}

// v = a * Dy * u  + b * v + c
template <template <typename> class Array, typename To = double,
          typename Ti = double, typename Tc = double>
void dmy42(Array<To>& out, const Array<Ti>& in, const Bounds& bounds, Tc a = 1.0,
          Tc b = 0.0, Tc c = 0.0) {
        dm42<Y>(out, in, bounds, a, b, c);
}

template <template <typename> class Array, typename To = double,
          typename Ti = double, typename Tc = double>
void dmy42(Array<To>& out, const Array<Ti>& in, Tc a = 1.0, Tc b = 0.0,
          Tc c = 0.0) {
        dm42<Y>(out, in, out.bounds(), a, b, c);
}

// v = a * Dz * u  + b * v + c
template <template <typename> class Array, typename To = double,
          typename Ti = double, typename Tc = double>
void dmz42(Array<To>& out, const Array<Ti>& in, const Bounds& bounds, Tc a = 1.0,
          Tc b = 0.0, Tc c = 0.0) {
        dm42<Z>(out, in, bounds, a, b, c);
}

template <template <typename> class Array, typename To = double,
          typename Ti = double, typename Tc = double>
void dmz42(Array<To>& out, const Array<Ti>& in, Tc a = 1.0, Tc b = 0.0,
          Tc c = 0.0) {
        dm42<Z>(out, in, out.bounds(), a, b, c);
}
