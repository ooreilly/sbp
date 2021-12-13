#pragma once

template <template <typename> class Array, typename Tv>
class Outer1D
{
        const Tv *weights_x;
        const size_t *indices_x;
        size_t num_basis, num_query, weights_line, indices_line, scalars_line;

        const Tv *u;
        size_t u_line, u_slice;

        const Tv *a;
        size_t a_line, a_slice;
        Tv b;

        Tv *v;
        size_t v_line, v_slice;

        public:
         Outer1D(Array<Tv>& out, const Array<Tv>& weights_x_,
                         const Array<size_t>& indices_x_,
                         const Array<Tv>& a_, 
                         Tv b_=1.0
                         ) {

                 num_query = weights_x_.nx;
                 num_basis = weights_x_.ny;
                 weights_x = weights_x_();
                 indices_x = indices_x_();
                 weights_line = weights_x_.line;
                 indices_line = indices_x_.line;
                 a = a_();
                 a_line = a_.line;
                 b = b_;

                 v = out();
                 v_line = out.line;
                 v_slice = out.slice;
         }

         VOID_OPERATOR() const {
                 size_t pos = xi + v_line * yj + v_slice * zk;
                 Tv out = b * v[pos];
                 for (int i = 0; i < num_query; ++i) {
                         for (int r = 0; r < num_basis; ++r) {
                                 size_t idx = indices_x[i + indices_line * r];
                                 if (idx != xi) continue;
                                 Tv wx = weights_x[i + weights_line * r];
                                 out += wx * a[i];
                         }
                 }

                 v[pos] = out;
         }

};

template <template <typename> class Array, typename Tv>
class Outer2D
{
        const Tv *weights_x, *weights_y;
        const size_t *indices_x, *indices_y;
        size_t num_basis, num_query, weights_line, indices_line, scalars_line;

        const Tv *u;
        size_t u_line, u_slice;

        const Tv *a;
        size_t a_line, a_slice;
        Tv b;

        Tv * v;
        size_t v_line, v_slice;

        public:
         Outer2D(
                         Array<Tv>& out,
                         const Array<Tv>& weights_x_,
                         const Array<Tv>& weights_y_,
                         const Array<size_t>& indices_x_,
                         const Array<size_t>& indices_y_,
                         const Array<Tv>& a_, 
                         Tv b_=1.0
                         ) {

                 num_query = weights_x_.nx;
                 num_basis = weights_x_.ny;
                 weights_x = weights_x_();
                 weights_y = weights_y_();
                 indices_x = indices_x_();
                 indices_y = indices_y_();
                 weights_line = weights_x_.line;
                 indices_line = indices_x_.line;

                 a = a_();
                 a_line = a_.line;
                 b = b_;

                 v = out();
                 v_line = out.line;
                 v_slice = out.slice;
         }

         VOID_OPERATOR() const {
                 size_t pos = xi + v_line * yj + v_slice * zk;
                 Tv out = b * v[pos];
                 for (int i = 0; i < num_query; ++i) {
                         for (int q = 0; q < num_basis; ++q) {
                                 size_t idy = indices_y[i + indices_line * q];
                                 Tv wy = weights_y[i + weights_line * q];
                                 if (idy != yj) continue;
                         for (int r = 0; r < num_basis; ++r) {
                                 size_t idx = indices_x[i + indices_line * r];
                                 if (idx != xi) continue;
                                 Tv wx = weights_x[i + weights_line * r];
                                 out += wx * wy * a[i];
                         }
                         }
                 }


                 v[pos] = out;
         }

};

template <template <typename> class Array, typename Tv>
class Outer3D
{
        const Tv *weights_x, *weights_y, *weights_z;
        const size_t *indices_x, *indices_y, *indices_z;
        size_t num_basis, num_query, weights_line, indices_line, scalars_line;

        const Tv *u;
        size_t u_line, u_slice;

        const Tv *a;
        size_t a_line, a_slice;
        Tv b;

        Tv * v;
        size_t v_line, v_slice;

        public:
         Outer3D(
                         Array<Tv>& out,
                         const Array<Tv>& weights_x_,
                         const Array<Tv>& weights_y_,
                         const Array<Tv>& weights_z_,
                         const Array<size_t>& indices_x_,
                         const Array<size_t>& indices_y_,
                         const Array<size_t>& indices_z_,
                         const Array<Tv>& a_, 
                         Tv b_=1.0
                         ) {

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

                 a = a_();
                 a_line = a_.line;
                 b = b_;

                 v = out();
                 v_line = out.line;
                 v_slice = out.slice;
         }

         VOID_OPERATOR() const {
                 size_t pos = xi + v_line * yj + v_slice * zk;
                 Tv out = b * v[pos];
                 for (int i = 0; i < num_query; ++i) {
                         for (int p = 0; p < num_basis; ++p) {
                                 size_t idz = indices_z[i + indices_line * p];
                                 Tv wz = weights_z[i + weights_line * p];
                                 if (idz != zk) continue;
                         for (int q = 0; q < num_basis; ++q) {
                                 size_t idy = indices_y[i + indices_line * q];
                                 Tv wy = weights_y[i + weights_line * q];
                                 if (idy != yj) continue;
                         for (int r = 0; r < num_basis; ++r) {
                                 size_t idx = indices_x[i + indices_line * r];
                                 if (idx != xi) continue;
                                 Tv wx = weights_x[i + weights_line * r];
                                 out += wx * wy * wz * a[i];
                         }
                         }
                         }
                 }

                 v[pos] = out;
         }

};

