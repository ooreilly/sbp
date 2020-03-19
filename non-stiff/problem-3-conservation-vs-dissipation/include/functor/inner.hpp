#pragma once

template <template <typename> class Array, typename Tv>
class Inner1D
{
        const Tv *weights_x;
        const size_t *indices_x;
        size_t num_basis, num_query, weights_line, indices_line;
        const Tv *u;

        Tv *v;
        size_t v_line, v_slice;

        public:
         Inner1D(Array<Tv>& out, const Array<Tv>& weights_x_,
                                const Array<size_t>& indices_x_,
                                const Array<Tv>& in) {

                 num_query = weights_x_.nx;
                 num_basis = weights_x_.ny;
                 weights_x = weights_x_();
                 indices_x = indices_x_();
                 weights_line = weights_x_.line;
                 indices_line = indices_x_.line;
                 u = in();

                 v = out();
                 v_line = out.line;
                 v_slice = out.slice;
         
         }

         VOID_OPERATOR() const {
                 // xi: number of queries

                 Tv out = 0.0;
                 for (int r = 0; r < num_basis; ++r) {
                         size_t idx = indices_x[xi + indices_line * r];
                         Tv w = weights_x[xi + weights_line * r];
                         out += u[idx] * w;
                 }
                 size_t pos = xi + v_line * yj + v_slice * zk;
                 v[pos] = out;
         }

};


template <template <typename> class Array, typename Tv>
class Inner2D
{
        const Tv *weights_x, *weights_y;
        const size_t *indices_x, *indices_y;
        size_t num_basis, num_query, weights_line, indices_line, u_line;
        const Tv *u;

        Tv *v;
        size_t v_line, v_slice;

        public:
         Inner2D(Array<Tv>& out, const Array<Tv>& weights_x_,
                 const Array<Tv>& weights_y_, const Array<size_t>& indices_x_,
                 const Array<size_t>& indices_y_, const Array<Tv>& in) {
                 num_query = weights_x_.nx;
                 num_basis = weights_x_.ny;
                 weights_x = weights_x_();
                 weights_y = weights_y_();
                 indices_x = indices_x_();
                 indices_y = indices_y_();
                 weights_line = weights_x_.line;
                 indices_line = indices_x_.line;
                 u_line = in.line; 
                 u = in();
                 u_line = in.line;

                 v = out();
                 v_line = out.line;
                 v_slice = out.slice;
         
         }

         VOID_OPERATOR() const {

                 Tv out = 0.0;
                 for (int p = 0; p < num_basis; ++p) {
                         size_t idy = indices_y[xi + indices_line * p];
                         Tv wy = weights_x[xi + weights_line * p];
                         for (int r = 0; r < num_basis; ++r) {
                                 size_t idx = indices_x[xi + indices_line * r];
                                 Tv wx = weights_x[xi + weights_line * r];
                                 out += u[idx + u_line * idy] * wx * wy;
                         }
                 }

                 size_t pos = xi + v_line * yj + v_slice * zk;
                 v[pos] = out;
                
         }

};

template <template <typename> class Array, typename Tv>
class Inner3D
{
        const Tv *weights_x, *weights_y, *weights_z, *u;
        const size_t *indices_x, *indices_y, *indices_z;
        size_t num_basis, num_query, weights_line, indices_line, u_line, u_slice;

        Tv *v;
        size_t v_line, v_slice;

       public:
        Inner3D(
            Array<Tv>& out, 
            const Array<Tv>& weights_x_, 
            const Array<Tv>& weights_y_,
            const Array<Tv>& weights_z_, 
            const Array<size_t>& indices_x_,
            const Array<size_t>& indices_y_,
            const Array<size_t>& indices_z_,
            const Array<Tv>& in) {
                num_query = weights_x_.nx;
                num_basis = weights_x_.ny;
                weights_x = weights_x_();
                weights_y = weights_y_();
                weights_z = weights_z_();
                indices_x = indices_x_();
                indices_y = indices_y_();
                indices_z = indices_z_();
                weights_line = weights_x_.line;
                indices_line = indices_x_.line;

                u_line = in.line;
                u_slice = in.slice;
                u = in();

                v = out();
                v_line = out.line;
                v_slice = out.slice;
        }

        VOID_OPERATOR() const {
                Tv out = 0.0;
                for (int p = 0; p < num_basis; ++p) {
                        size_t idz = indices_z[xi + indices_line * p];
                        Tv wz = weights_z[xi + weights_line * p];
                        for (int q = 0; q < num_basis; ++q) {
                                size_t idy = indices_y[xi + indices_line * q];
                                Tv wy = weights_x[xi + weights_line * q];
                                for (int r = 0; r < num_basis; ++r) {
                                        size_t idx =
                                            indices_x[xi + indices_line * r];
                                        Tv wx =
                                            weights_x[xi + weights_line * r];
                                        out += u[idx + u_line * idy +
                                                 u_slice * idz] *
                                               wx * wy * wz;
                                }
                        }
                }

                size_t pos = xi + v_line * yj + v_slice * zk;
                v[pos] = out;
         }

};


