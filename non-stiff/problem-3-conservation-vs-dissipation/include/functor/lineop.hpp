#pragma once
#include "macro/operator.hpp"
#include "function/lineop.hpp"

// v = a * D * u + b * v + c
template <int dim, template <typename> class Array, int n_left, int n_right, int n_interior, typename Spec, typename Tvo,
          typename Tvi, typename Tc>
class LineOp {
        Tc a, b, c;
        Tvo *v;
        size_t v_line, v_slice;
        
        const Tvi *u;
        size_t u_line, u_slice;
        size_t u_nx, u_ny, u_nz;

        Spec spec;

       public:
        LineOp(Spec spec_, Array<Tvo>& v_, const Array<Tvi>& u_, 
                Tc a_ = 1.0, Tc b_ = 0.0, Tc c_ = 0.0) {

                spec = spec_;

                v = v_();
                v_line = v_.line;
                v_slice = v_.slice;

                u = u_();
                u_line = u_.line;
                u_slice = u_.slice;
                u_nx = u_.nx;
                u_ny = u_.ny;
                u_nz = u_.nz;

                a = a_;
                b = b_;
                c = c_;
        }

        VOID_OPERATOR() const { 
                size_t pos = xi + v_line * yj + v_slice * zk;
                v[pos] = a * lineop<dim, Tvo, Tvi, n_left, n_right, n_interior,
                                    Spec>(u, xi, yj, zk, spec, u_nx, u_ny, u_nz,
                                          u_line, u_slice) +
                         b * v[pos] + c;
        }
};

