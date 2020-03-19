#pragma once
#include "macro/operator.hpp"
#include "macro/dim.hpp"


// v = alpha * D * u + beta * v + gamma
template <int dim, typename Operator, typename Tv>
class AxpyOp : public Operator {
       public:
        Tv alpha, beta, gamma;
        AxpyOp(Tv hi_, Tv alpha_ = 1.0, Tv beta_ = 0.0, Tv gamma_ = 0.0) {
                this->hi = hi_;
                alpha = alpha_;
                beta = beta_;
                gamma = gamma_;
        }

        OPERATOR() const {
                Tv op_u =0.0;
                if (dim == X)
                        op_u = this->x(xi, yj, zk, ref1, ref2, out_val, in_val1, in_val2, nx, ny, nz, line,
                                         slice, inplace);
                if (dim == Y)
                        op_u = this->y(xi, yj, zk, ref1, ref2, out_val, in_val1, in_val2, nx, ny, nz, line,
                                         slice, inplace);
                if (dim == Z)
                        op_u = this->z(xi, yj, zk, ref1, ref2, out_val, in_val1, in_val2,nx, ny, nz, line,
                                         slice, inplace);
                return alpha * op_u + beta * out_val + gamma;
        }
};
