#pragma once

#include "function/dp42.hpp"
#include "function/dm42.hpp"
#include "function/ricker.hpp"
#include "function/expr.hpp"
#include "function/fill.hpp"
#include "function/hp4.hpp"
#include "function/hm4.hpp"
#include "function/bounds.hpp"
#include "function/sides.hpp"
#include "function/gather.hpp"
#include "function/scatter.hpp"
#include "function/constant.hpp"
#include "function/multiply.hpp"
#include "transform/apply.hpp"
#include "structure/lagrange.hpp"
#include "macro/sides.hpp"

#define NUM_BLOCKS 2
#define NUM_SIDES 4
#define NUM_SOURCES_PER_BLOCK 1


typedef double Tprec;
typedef DeviceArray<Tprec> ComputeArray;
typedef Lagrange<DeviceArray, Tprec> ComputeLagrange;

// Fields
ComputeArray p[NUM_BLOCKS];
ComputeArray vx[NUM_BLOCKS];
ComputeArray vy[NUM_BLOCKS];

// Temporary storage
ComputeArray tmp[NUM_BLOCKS];

// Material properties (wave speed)
Tprec c[NUM_BLOCKS] = {1.0, 4.0};
Tprec rho[NUM_BLOCKS] = {1.0, 800.0};

ComputeArray c_p[NUM_BLOCKS];
ComputeArray c_vx[NUM_BLOCKS];
ComputeArray c_vy[NUM_BLOCKS];

// Boundary values
ComputeArray load_side_p[NUM_BLOCKS][NUM_SIDES], 
             load_side_vx[NUM_BLOCKS][NUM_SIDES], 
             load_side_vy[NUM_BLOCKS][NUM_SIDES];
ComputeArray store_side_p[NUM_BLOCKS][NUM_SIDES], 
             store_side_vx[NUM_BLOCKS][NUM_SIDES], 
             store_side_vy[NUM_BLOCKS][NUM_SIDES];



// Material properties on the boundary
ComputeArray side_Z[NUM_BLOCKS][NUM_SIDES], 
             side_rho[NUM_BLOCKS][NUM_SIDES];

// Bounds
Bounds bounds_p[NUM_BLOCKS];
Bounds bounds_vx[NUM_BLOCKS];
Bounds bounds_vy[NUM_BLOCKS];

Bounds bounds_sides_p[NUM_BLOCKS][NUM_SIDES];
Bounds bounds_sides_vx[NUM_BLOCKS][NUM_SIDES];
Bounds bounds_sides_vy[NUM_BLOCKS][NUM_SIDES];

Bounds bounds_source[NUM_BLOCKS][NUM_SOURCES_PER_BLOCK];

// Grid spacing and reciprocal grid spacing
Tprec hx[NUM_BLOCKS];
Tprec hy[NUM_BLOCKS];
Tprec hxi[NUM_BLOCKS];
Tprec hyi[NUM_BLOCKS];

Tprec normal[NUM_SIDES] = {-1, 1, -1, 1};
Tprec normal_x[NUM_SIDES] = {-1, 1, 0, 0};
Tprec normal_y[NUM_SIDES] = {0, 0, -1, 1};


// Grids
ComputeArray x_p[NUM_BLOCKS];
ComputeArray y_p[NUM_BLOCKS];
ComputeArray z_p[NUM_BLOCKS];

ComputeArray x_vx[NUM_BLOCKS];
ComputeArray y_vx[NUM_BLOCKS];
ComputeArray z_vx[NUM_BLOCKS];

ComputeArray x_vy[NUM_BLOCKS];
ComputeArray y_vy[NUM_BLOCKS];
ComputeArray z_vy[NUM_BLOCKS];

// 1D grids needed for sources and receivers
ComputeArray x1_p[NUM_BLOCKS];
ComputeArray y1_p[NUM_BLOCKS];

int num_blocks;

// Source
const size_t num_source = 1;
int active_source[2][1] = {{1},{0}};
int degree = 3;
Tprec t0[2][1] = {{0.5}, {0.5}};
Tprec fp[2][1] = {{8}, {8}};
Tprec sx[2][num_source] = {{0.75}, {0.5}};
Tprec sy[2][num_source] = {{0.5}, {0.5}};
typedef struct {
        Tprec x[1];
        Tprec y[1];
} source_t;

source_t source_desc[NUM_BLOCKS];

// Hold source time function at each time step for each source
ComputeArray source_f[NUM_BLOCKS];

ComputeArray source_xarr[NUM_BLOCKS];
ComputeArray source_yarr[NUM_BLOCKS];

ComputeLagrange source_L[NUM_BLOCKS];

// Receiver
const size_t num_recv = 1;
Tprec rx[2] = {0.25, 0};
Tprec ry[2] = {0.5, 0};
typedef struct {
        Tprec x[1];
        Tprec y[1];
} recv_t;
recv_t recv_desc[NUM_BLOCKS];
ComputeLagrange recv_L[NUM_BLOCKS];

ComputeArray recv_xarr[NUM_BLOCKS];
ComputeArray recv_yarr[NUM_BLOCKS];

// SAT quadrature weights
Tprec sat_p[NUM_BLOCKS][NUM_SIDES];
Tprec sat_vn[NUM_BLOCKS][NUM_SIDES];

// Boundary conditions
const int BND_PRESSURE = 0;
const int COUPLING = 1;
const int ENERGY_CONSERVATIVE = 2;
const int ENERGY_DISSIPATIVE = 3;
const int NAIVE_COUPLING = 4;
const int ENERGY_CONSERVATIVE_OPT = 5;
const int ENERGY_CONSERVATIVE_SPLIT = 6;
const int NO_COUPLING = -1;

// left, right, bottom, top
int boundary_type[2][NUM_SIDES];

int coupling_block[2][NUM_SIDES];
int coupling_side[2][NUM_SIDES];
int coupling_method[2][NUM_SIDES];


void set_default(int out[NUM_BLOCKS][NUM_SIDES], int value)
{
        for (int block = 0; block < NUM_BLOCKS; ++block)
                for (int side = 0; side < NUM_SIDES; ++side)
                        out[block][side] = value;

}

void couple(int out[NUM_BLOCKS][NUM_SIDES], int side1, int side2, int value1,
            int value2 = NO_COUPLING) {
        out[0][side1] = value1;
        if (value2 == NO_COUPLING) value2 = value1;
        out[1][side2] = value2;
}

void init(size_t nx[], size_t ny[], Tprec hx_[], Tprec hy_[], int coupling, int num_blocks_=2) {
        num_blocks = num_blocks_;

        set_default(boundary_type, BND_PRESSURE);
        set_default(coupling_block, NO_COUPLING);
        set_default(coupling_side, NO_COUPLING);
        set_default(coupling_method, NO_COUPLING);

        couple(boundary_type, RIGHT, LEFT, COUPLING);
        couple(coupling_block, RIGHT, LEFT, 1, 0);
        couple(coupling_side, RIGHT, LEFT, LEFT, RIGHT);
        printf("coupling = %d \n", coupling);
        couple(coupling_method, RIGHT, LEFT, coupling);

        for (int i = 0; i < num_blocks; ++i) { 
                hx[i] = hx_[i];
                hy[i] = hy_[i];
                hxi[i] = 1.0 / hx_[i];
                hyi[i] = 1.0 / hy_[i];

                // Fields
                p[i] = ComputeArray(nx[i], ny[i]);
                vx[i] = ComputeArray(nx[i], ny[i]);
                vy[i] = ComputeArray(nx[i], ny[i]);

                // Temporary storage
                tmp[i] = ComputeArray(nx[i], ny[i]);

                // Material properties
                c_p[i] = ComputeArray(nx[i], ny[i]);
                c_vx[i] = ComputeArray(nx[i], ny[i]);
                c_vy[i] = ComputeArray(nx[i], ny[i]);

                constant(c_p[i], rho[i] * c[i] * c[i]);
                constant(c_vx[i], 1.0 / rho[i]);
                constant(c_vy[i], 1.0 / rho[i]);

                // Bounds
                size_t nxm = nx[i];
                size_t nym = ny[i];
                size_t nxp = nx[i] - 1;
                size_t nyp = ny[i] - 1;

                bounds_p[i] = bounds(nxm, nym);
                bounds_vx[i] = bounds(nxp, nym);
                bounds_vy[i] = bounds(nxm, nyp);

                sides_bounds(bounds_sides_p[i], NUM_SIDES, nxm, nym);
                sides_bounds(bounds_sides_vx[i], NUM_SIDES, nxp, nym);
                sides_bounds(bounds_sides_vy[i], NUM_SIDES, nxm, nyp);

                sides_array(load_side_p[i], NUM_SIDES, nxm, nym);
                sides_array(load_side_vx[i], NUM_SIDES, nxp, nym);
                sides_array(load_side_vy[i], NUM_SIDES, nxm, nyp);
                sides_array(store_side_p[i], NUM_SIDES, nxm, nym);
                sides_array(store_side_vx[i], NUM_SIDES, nxp, nym);
                sides_array(store_side_vy[i], NUM_SIDES, nxm, nyp);
                sides_array(side_Z[i], NUM_SIDES, nxm, nym);
                sides_array(side_rho[i], NUM_SIDES, nxm, nym);

                for (int side = 0; side < NUM_SIDES; ++side) {
                        constant(side_Z[i][side], rho[i] * c[i]);
                        constant(side_rho[i][side], rho[i]);
                }

                // SAT weights
                Tprec hm = hm4<Tprec>(); 
                Tprec hp = hp4<Tprec>(); 
                sat_p[i][LEFT]   = hxi[i] / hm;
                sat_p[i][RIGHT]  = hxi[i] / hm;
                sat_p[i][TOP]    = hyi[i] / hm;
                sat_p[i][BOTTOM] = hyi[i] / hm;

                sat_vn[i][LEFT]   = hxi[i] / hp;
                sat_vn[i][RIGHT]  = hxi[i] / hp;
                sat_vn[i][TOP]    = hyi[i] / hp;
                sat_vn[i][BOTTOM] = hyi[i] / hp;



                // 1D grids
                x1_p[i] = ComputeArray(nx[i]);
                y1_p[i] = ComputeArray(ny[i]);
                xm(x1_p[i], hx[i]);
                xm(y1_p[i], hy[i]);

                // Source
                for (int j = 0; j < num_source; ++j) {
                        source_desc[i].x[j] = sx[i][j];
                        source_desc[i].y[j] = sy[i][j];
                }

                source_xarr[i] = ComputeArray(num_source);
                source_yarr[i] = ComputeArray(num_source);

                fill(source_xarr[i], source_desc[i].x, num_source);
                fill(source_yarr[i], source_desc[i].y, num_source);

                source_f[i] = ComputeArray(num_source);
                source_L[i] = ComputeLagrange(source_xarr[i], source_yarr[i],
                                              x1_p[i], y1_p[i], degree);

                // Receiver
                recv_desc[i].x[0] = rx[i];
                recv_desc[i].y[0] = ry[i];

                recv_xarr[i] = ComputeArray(num_recv);
                recv_yarr[i] = ComputeArray(num_recv);

                fill(recv_xarr[i], recv_desc[i].x, num_recv);
                fill(recv_yarr[i], recv_desc[i].y, num_recv);

                recv_L[i] = ComputeLagrange(recv_xarr[i], recv_yarr[i],
                                            x1_p[i], y1_p[i], degree);

                
        }

        printf("Initialization completed\n");



}

void init_grids(size_t nx[], size_t ny[], Tprec hx[], Tprec hy[], Tprec Lx, Tprec Ly) {
        for (int i = 0; i < num_blocks; ++i) { 
                // Grids
                x_p[i] = ComputeArray(nx[i], ny[i]);
                y_p[i] = ComputeArray(nx[i], ny[i]);
                z_p[i] = ComputeArray(nx[i], ny[i]);

                xm(x_p[i], bounds_p[i], hx[i]);
                ym(y_p[i], bounds_p[i], hy[i]);


                x_vx[i] = ComputeArray(nx[i], ny[i]);
                y_vx[i] = ComputeArray(nx[i], ny[i]);
                z_vx[i] = ComputeArray(nx[i], ny[i]);

                xp(x_vx[i], bounds_vx[i], hx[i]);
                ym(y_vx[i], bounds_vx[i], hy[i]);

                x_vy[i] = ComputeArray(nx[i], ny[i]);
                y_vy[i] = ComputeArray(nx[i], ny[i]);
                z_vy[i] = ComputeArray(nx[i], ny[i]);

                xm(x_vy[i], bounds_vy[i], hx[i]);
                yp(y_vy[i], bounds_vy[i], hy[i]);

        }

        // Shift each grid +1 in the x-direction
        expr(x_p[1], bounds_p[1], 1.0, Lx);
        expr(x_vx[1], bounds_vx[1], 1.0, Lx);
        expr(x_vy[1], bounds_vy[1], 1.0, Lx);

}

std::string vtkfilename(const char *prefix, const char *var, int block,
                        int step) {
        std::string name = std::string(prefix) + std::string(var) +
                           "_block_" + std::to_string(block) + "_step_" +
                           std::to_string(step) + ".vtk";
        return name;
}

void write_fields(const char *prefix, int step=0) {
        for (int i = 0; i < num_blocks; ++i) {
                // p
                std::string filename = vtkfilename(prefix, "p", i, step);
                vtk(filename.c_str(), x_p[i].host(), y_p[i].host(),
                    z_p[i].host(), bounds_p[i]);
                vtk(filename.c_str(), p[i].host(), bounds_p[i], "p");

                // vx
                filename = vtkfilename(prefix, "vx", i, step);
                vtk(filename.c_str(), x_vx[i].host(), y_vx[i].host(),
                    z_vx[i].host(), bounds_vx[i]);
                vtk(filename.c_str(), vx[i].host(), bounds_vx[i], "vx");

                // vy
                filename = vtkfilename(prefix, "vy", i, step);
                vtk(filename.c_str(), x_vy[i].host(), y_vy[i].host(),
                    z_vy[i].host(), bounds_vy[i]);
                vtk(filename.c_str(), vy[i].host(), bounds_vy[i], "vy");
        }
}

// unsigned normal velocity
ComputeArray* normal_velocity(ComputeArray *vx, ComputeArray *vy, int side) {
        ComputeArray *vn = NULL;
        if (side == LEFT || side == RIGHT)
                vn = vx;
        else
                vn = vy;
        return vn;
        return NULL;
}

template <typename Tv>
__inline__ __device__ Tv coupling_param(int coupling, Tv Z1, Tv Z2, Tv rho1, Tv rho2, Tv p1, Tv p2, int eq)
{
        Tv f1 = pow(Z1, 2.0);
        Tv f2 = pow(Z2, 2.0);

        if (coupling == ENERGY_CONSERVATIVE_OPT) {
                if (eq == 1) return f2 / (f1 + f2);
                if (eq == 2) return f1 / (f1 + f2);
                return 0;
        }

        if (coupling == ENERGY_CONSERVATIVE_SPLIT) {
                if (eq == 1) {
                        if (Z2 > Z1) return 1.0;
                        if (Z2 < Z1) return 0.0;
                        if (Z2 == Z1) return 0.5;
                }
                if (eq == 2) {
                        if (Z2 > Z1) return 0.0;
                        if (Z2 < Z1) return 1.0;
                        if (Z2 == Z1) return 0.5;
                }
        }

        return 0;
}

template <typename Tv>
class PressureCoupling {

        Tv *p1_out;
        const Tv *p1, *p2, *vn1, *vn2, *rho1, *rho2, *Z1, *Z2;
        int line;
        int coupling;
        Tv sat, normal;

        public:
         PressureCoupling(ComputeArray &p1_out, 
                         const ComputeArray &p1, 
                         const ComputeArray &p2, const
                         ComputeArray &vn1, const ComputeArray &vn2, const
                         ComputeArray &rho1, 
                         ComputeArray &rho2, 
                         const ComputeArray &Z1, const
                         ComputeArray &Z2, Tv sat, Tv normal, 
                          int coupling = 0)
             : 
               p1_out(p1_out()),
               p1(p1()),
               p2(p2()),
               vn1(vn1()),
               vn2(vn2()),
               rho1(rho1()),
               rho2(rho2()),
               Z1(Z1()),
               Z2(Z2()), 
               line(p1.line), 
               sat(sat), 
               normal(normal), 
               coupling(coupling)
        {}

         // p = sat * (vn1 - vn2) * Z1 * Z2 / (rho1 * (Z1 + Z2))
         VOID_OPERATOR()  {
                size_t pos = xi + line * yj;
                Tv k1 = Z2[pos] / (Z1[pos] + Z2[pos]);

                Tv K = Z1[pos] * Z1[pos] / rho1[pos];

                if (coupling == NAIVE_COUPLING) k1 = 0.5;
                if (coupling == ENERGY_CONSERVATIVE_OPT || coupling == ENERGY_CONSERVATIVE_SPLIT) {
                        k1 = coupling_param(coupling, Z1[pos], Z2[pos], rho1[pos], rho2[pos], p1[pos], p2[pos], 1);

                }
                Tv a = k1;

                p1_out[pos] = K * sat * normal * (vn1[pos] - vn2[pos]) * a;

                if (coupling == ENERGY_DISSIPATIVE) {
                        p1_out[pos] += -K * sat * (p1[pos] - p2[pos]) * Z1[pos] /
                                       (Z1[pos] + Z2[pos]) / rho1[pos];
                }
         }
};

template <typename Tv>
class VelocityCoupling {

        Tv *vn1_out;
        const Tv *p1, *p2, *vn1, *vn2, *rho1, *rho2, *Z1, *Z2;
        int line;
        int coupling;
        Tv sat, normal;

        public:
         VelocityCoupling(ComputeArray &vn1_out, const ComputeArray &p1,
                          const ComputeArray &p2, const ComputeArray &vn1,
                          const ComputeArray &vn2, 
                          const ComputeArray &rho1, const ComputeArray &rho2,
                          const ComputeArray &Z1, const ComputeArray &Z2,
                          Tv sat, Tv normal, int coupling = 0)
             : vn1_out(vn1_out()),
               p1(p1()),
               p2(p2()),
               vn1(vn1()),
               vn2(vn2()),
               rho1(rho1()),
               rho2(rho2()),
               Z1(Z1()),
               Z2(Z2()),
               line(p1.line),
               sat(sat),
               normal(normal),
               coupling(coupling) {}

         // vn = sat * (p1 - p2) * Z1 / (rho1 * (Z1 + Z2))
         VOID_OPERATOR()  {
                size_t pos = xi + line * yj;
                Tv k2 = Z1[pos] / (Z1[pos] + Z2[pos]);
                if (coupling == NAIVE_COUPLING) k2 = 0.5;
                if (coupling == ENERGY_CONSERVATIVE_OPT || coupling == ENERGY_CONSERVATIVE_SPLIT)
                        k2 = coupling_param(coupling, Z1[pos], Z2[pos], rho1[pos], rho2[pos], p1[pos], p2[pos], 2);

                Tv a = k2 / rho1[pos];
                // This might look confusing because of Z1 * p1, but recall that
                // the code treats p as p / Z.
                vn1_out[pos] = sat * normal * a * 
                                (p1[pos] - p2[pos]);
                if (coupling == ENERGY_DISSIPATIVE) {
                        vn1_out[pos] += - sat * a * Z2[pos] * 
                                      (vn1[pos] - vn2[pos]);
                }
         }
};

template <typename Tv>
class PressureBoundaryCondition {

        Tv *v1_out;
        const Tv *p1, *rho1, *Z1;
        int line;
        Tv sat, normal, data;

        public:
         PressureBoundaryCondition(ComputeArray &v1_out, 
                         const ComputeArray &p1, 
                         ComputeArray &rho1, const ComputeArray &Z1, 
                         Tv sat, Tv normal, Tv data)
             : 
               v1_out(v1_out()),
               p1(p1()),
               rho1(rho1()),
               Z1(Z1()),
               line(p1.line), 
               sat(sat), 
               normal(normal)
        {}

         VOID_OPERATOR()  {
                size_t pos = xi + line * yj;
                v1_out[pos] = sat * normal * 1.0 / rho1[pos] * (p1[pos] - data);
         }
};



template <typename Tv>
class Source {
        public:
        void operator()(ComputeArray& dpi, const Tv t, const Tv dt, const int
                        block) {
                // Apply sources
                for (int j = 0; j < num_source; ++j) {
                        if (active_source[block][j]) {
                                Tv g = ricker(t, t0[block][j], fp[block][j]);
                                // dp = dp + dt / (hx * hy) * g(t)
                                Tv data[] = {dt * hxi[block] * hyi[block] * g};
                                fill(source_f[block], data, num_source);
                                source_L[block].outer(dpi, source_f[block]);
                        }
                }
        }

};

template <typename Tv>
class Pressure {

        Source<Tv> source;

        public:
                Pressure(){}

                void operator()(const Tv t, const Tv dt) {
                        for (int i = 0; i < num_blocks; ++i) {
                                gather_data(vx[i], vy[i], i);
                        }
                        
                        for (int i = 0; i < num_blocks; ++i) {
                                update(p[i], vx[i], vy[i], i, t, dt);
                        }
                }

                void update(ComputeArray &dpi, 
                            ComputeArray &dvxi,
                            ComputeArray &dvyi,
                                
                            const int i, const Tv t,
                            const Tv dt, const Tv alpha = 1.0) {
                        // dp = alpha * dp - dt * Dx * vx - dt * Dy * vy + dt *
                        // s(t)
                        // tmp = - dt * Dx * vx
                        dmx42(tmp[i], dvxi, bounds_p[i], -dt * hxi[i]);
                        // tmp = - dt * Dx * vx - dt * Dy * vy
                        dmy42(tmp[i], dvyi, bounds_p[i], -dt * hyi[i],
                              (Tprec)1.0);
                        // p[i] = alpha * dp + c . * tmp
                        multiply(dpi, c_p[i], tmp[i], bounds_p[i], (Tprec)1.0,
                                 alpha);
                        source(dpi, t, dt, i);

                        boundary(i);
                        bscatter(dpi, i, dt);


                }

                void boundary(int block) {
                        // Apply boundary or coupling conditions for block
                        // All types of boundary conditions are called, but the
                        // function exists immediately if the boundary condition
                        // is not supposed to be applied

                        for (int side = 0; side < NUM_SIDES; ++side) {
                                coupling(boundary_type[block][side], block,
                                         side, coupling_block[block][side],
                                         coupling_side[block][side],
                                         sat_p[block][side], normal[side]);
                        }
                }

                void coupling(int type, int block, int side, int coupling_block,
                              int coupling_side, Tv sat, Tv normal) {
                        if (coupling_block == NO_COUPLING) return;

                        ComputeArray *vn1 =
                            normal_velocity(&(load_side_vx[block][side]),
                                            &(load_side_vy[block][side]), side);
                        ComputeArray *vn2 = normal_velocity(
                            &(load_side_vx[coupling_block][coupling_side]),
                            &(load_side_vy[coupling_block][coupling_side]),
                            coupling_side);

                        PressureCoupling<Tv> pc(
                            store_side_p[block][side],
                            load_side_p[block][side], 
                            load_side_p[coupling_block][coupling_side], 
                            *vn1, *vn2,
                            side_rho[block][side], 
                            side_rho[coupling_block][coupling_side], 
                            side_Z[block][side],
                            side_Z[coupling_block][coupling_side], sat, normal,
                            coupling_method[block][side]);

                        apply(pc, store_side_p[block][side].bounds(),
                                        store_side_p[block][side]);
                }

                void gather_data(const ComputeArray& vxi, const ComputeArray& vyi, int block) {
                        for (int side = 0; side < NUM_SIDES; ++side) {
                                gather(load_side_vx[block][side], vxi,
                                       bounds_sides_vx[block][side]);
                                gather(load_side_vy[block][side], vyi,
                                       bounds_sides_vy[block][side]);
                        }
                }

                void bscatter(ComputeArray &dvpi, int block, Tv dt) {
                        for (int side = 0; side < NUM_SIDES; ++side)
                                scatter(dvpi, store_side_p[block][side],
                                        bounds_sides_p[block][side], dt, (Tv)1.0);
                }
};

template <typename Tv>
class Velocity {

        public:
         Velocity() {
                }

        void operator()(const Tv t, const Tv dt, const Tv alpha=1.0) {

                for (int i = 0; i < num_blocks; ++i) {
                        gather_data(p[i], i);
                }

                // vx = alpha * vx - dt * c * Dx * p
                // vy = alpha * vy - dt * c * Dy * p
                for (int i = 0; i < num_blocks; ++i) {
                        update(vx[i], vy[i], p[i], i, t, dt, alpha);
                }

        }

        void update(
                    ComputeArray &dvxi, ComputeArray &dvyi, 
                    ComputeArray &dpi,
                    const int i,
                    const Tv t, const Tv dt, const Tv alpha = 1.0) {

                // tmp = - dt * Dx * p
                dpx42(tmp[i], dpi, bounds_vx[i], -dt * hxi[i]);
                // dvx = alpha * dvx + c .* tmp
                multiply(dvxi, c_vx[i], tmp[i], bounds_vx[i], (Tv)1.0, alpha);

                // tmp = - dt * Dy * p
                dpy42(tmp[i], dpi, bounds_vy[i], -dt * hyi[i]);
                // dvy = alpha * dvy + c .* tmp
                multiply(dvyi, c_vy[i], tmp[i], bounds_vy[i], (Tv)1.0, alpha);

                // Apply boundary conditions
                boundary(i);
                bscatter(dvxi, dvyi, i, dt);
        }

        void gather_data(const ComputeArray& pi, int block) {
                for (int side = 0; side < NUM_SIDES; ++side) {
                        gather(load_side_p[block][side], pi,
                               bounds_sides_p[block][side]);
                }
        }

        void bscatter(ComputeArray &dvxi, ComputeArray &dvyi, int block, Tv dt) {
                scatter(dvxi, store_side_vx[block][LEFT],
                        bounds_sides_vx[block][LEFT], dt, (Tprec)1.0);
                scatter(dvxi, store_side_vx[block][RIGHT],
                        bounds_sides_vx[block][RIGHT], dt, (Tprec)1.0);
                scatter(dvyi, store_side_vy[block][BOTTOM],
                        bounds_sides_vy[block][BOTTOM], dt, (Tprec)1.0);
                scatter(dvyi, store_side_vy[block][TOP],
                        bounds_sides_vy[block][TOP], dt, (Tprec)1.0);
        }



        void boundary(int block) {
                // Apply boundary or coupling conditions for block
                // All types of boundary conditions are called, but the function
                // exists immediately if the boundary condition is not supposed
                // to be applied

                for (int side = 0; side < NUM_SIDES; ++side) {

                        pressure(boundary_type[block][side], block, side,
                                 sat_vn[block][side], normal[side], 0.0);

                        coupling(boundary_type[block][side],
                                 block, side, 
                                 coupling_block[block][side],
                                 coupling_side[block][side],
                                 sat_vn[block][side], normal[side]);
                }
        }

        void pressure(int type, 
                      int block, int side, 
                      Tv sat, Tv normal, Tv data) {
                if (type != BND_PRESSURE) return;

                ComputeArray *side_vn;
                side_vn = normal_velocity(&(store_side_vx[block][side]),
                                          &(store_side_vy[block][side]), side);

                Bounds bounds = (*side_vn).bounds();
                ComputeArray *side_p = &load_side_p[block][side];
                PressureBoundaryCondition<Tv> pbc(
                    *side_vn, *side_p, side_rho[block][side],
                    side_Z[block][side], sat, normal, data);
                apply(pbc, side_vn->bounds(), *side_vn);
                //// vn = sat * side_p
                //expr(*side_vn, *side_p, bounds, normal * sat);
                //// vn = vn - dt * sat * data
                //expr(*side_vn, bounds, (Tv)1.0, -normal * sat * data);
        }

        void coupling(int type, int block, int side, int coupling_block,
                      int coupling_side, Tv sat, Tv normal) {
                if (coupling_block == NO_COUPLING) return;

                ComputeArray *store_vn1 =
                    normal_velocity(&(store_side_vx[block][side]),
                                    &(store_side_vy[block][side]), side);

                ComputeArray *load_vn1 =
                    normal_velocity(&(load_side_vx[block][side]),
                                    &(load_side_vy[block][side]), side);

                ComputeArray *load_vn2 =
                    normal_velocity(&(load_side_vx[coupling_block][coupling_side]),
                                    &(load_side_vy[coupling_block][coupling_side]), side);

                VelocityCoupling<Tv> vc(*store_vn1, 
                                load_side_p[block][side],
                                load_side_p[coupling_block][coupling_side],
                                *load_vn1,
                                *load_vn2,
                                side_rho[block][side], 
                                side_rho[coupling_block][coupling_side], 
                                side_Z[block][side],
                                side_Z[coupling_block][coupling_side], 
                                sat, normal, coupling_method[block][side]);
                apply(vc, store_vn1->bounds(), *store_vn1);
        }
};

template <typename Tv>
class RK4 {
        Pressure<Tv> *pre;
        Velocity<Tv> *vel;
        ComputeArray dp[NUM_BLOCKS];
        ComputeArray dvx[NUM_BLOCKS];
        ComputeArray dvy[NUM_BLOCKS];

        public:
         RK4(Pressure<Tv> &p_, Velocity<Tv> &v_) {
                 pre = &p_;
                 vel = &v_;

                 for (int i = 0; i < num_blocks; ++i) {
                         dp[i] = ComputeArray(p[i]);
                         dvx[i] = ComputeArray(vx[i]);
                         dvy[i] = ComputeArray(vy[i]);
                 }
         }

         void rates(Tv t, Tv alpha) {
                 for (int i = 0; i < num_blocks; ++i) {
                         vel->gather_data(p[i], i);
                         pre->gather_data(vx[i], vy[i], i);
                 }

                 for (int i = 0; i < num_blocks; ++i) {
                        pre->update(dp[i], vx[i], vy[i], i, t, (Tprec)1.0, alpha);
                        vel->update(dvx[i], dvy[i], p[i], i, t, (Tprec)1.0, alpha);
                 }
        }

         // update
         void update(Tv beta){ 
                for (int i = 0; i < num_blocks; ++i) {
                        // p = p + beta * dp
                        expr(p[i], dp[i], bounds_p[i], beta, (Tprec)1.0);
                        // vx = vx + beta * dvx
                        expr(vx[i], dvx[i], bounds_vx[i], beta, (Tprec)1.0);
                        // vy = vy + beta * dvy
                        expr(vy[i], dvy[i], bounds_vy[i], beta, (Tprec)1.0);
                }
         }
};


template <typename Tv>
class RK4StaggeredPressure {
        Pressure<Tv> *pre;
        Velocity<Tv> *vel;
        ComputeArray k1[NUM_BLOCKS];
        ComputeArray k2x[NUM_BLOCKS];
        ComputeArray k2y[NUM_BLOCKS];
        ComputeArray k3[NUM_BLOCKS];
        ComputeArray k4x[NUM_BLOCKS];
        ComputeArray k4y[NUM_BLOCKS];
        ComputeArray k5[NUM_BLOCKS];

       public:
        RK4StaggeredPressure(Pressure<Tv> &pre, Velocity<Tv> &vel)
            : pre(&pre), vel(&vel) {

                for (int i = 0; i < num_blocks; ++i) {
                        k1[i] = ComputeArray(p[i]);
                        k2x[i] = ComputeArray(p[i]);
                        k2y[i] = ComputeArray(p[i]);
                        k3[i] = ComputeArray(p[i]);
                        k4x[i] = ComputeArray(p[i]);
                        k4y[i] = ComputeArray(p[i]);
                        k5[i] = ComputeArray(p[i]);
                }
        }

        void operator()(Tv t, Tv dt) {
                for (int i = 0; i < num_blocks; ++i) {
                        pre->gather_data(vx[i], vy[i], i);
                        vel->gather_data(p[i], i);
                }

                //    k1 = dt * f(t+0.5*dt, v)
                for (int i = 0; i < num_blocks; ++i) {
                        pre->update(k1[i], vx[i], vy[i], i, t + 0.5 * dt, dt,
                                    0.0);
                }

                // k2 = dt * g(t, u)
                for (int i = 0; i < num_blocks; ++i) {
                        vel->update(k2x[i], k2y[i], p[i], i, t, dt, 0.0);
                }
                // k2 = v - k2
                for (int i = 0; i < num_blocks; ++i) {
                        expr(k2x[i], vx[i], k2x[i].bounds(), 1.0, -1.0);
                        expr(k2y[i], vy[i], k2y[i].bounds(), 1.0, -1.0);
                }

                for (int i = 0; i < num_blocks; ++i) {
                        pre->gather_data(k2x[i], k2y[i], i);
                }

                //    k3 = dt * f(t - 0.5 *dt, v - k2)
                for (int i = 0; i < num_blocks; ++i) {
                        pre->update(k3[i], k2x[i], k2y[i], i, t - 0.5 * dt, dt,
                                    0.0);
                }

                // k2 = u + k1
                for (int i = 0; i < num_blocks; ++i) {
                        expr(k2x[i], p[i], k1[i], k2x[i].bounds(), 1.0, 1.0);
                }

                for (int i = 0; i < num_blocks; ++i) {
                        vel->gather_data(k2x[i], i);
                }

                // k4 = dt * g(t + dt, u + k1)
                for (int i = 0; i < num_blocks; ++i) {
                        vel->update(k4x[i], k4y[i], k2x[i], i, t + dt, dt, 0.0);
                }

                // k4 = v + k4
                for (int i = 0; i < num_blocks; ++i) {
                        expr(k4x[i], vx[i], k4x[i], k4x[i].bounds(), 1.0, 1.0,
                             0.0);
                        expr(k4y[i], vy[i], k4y[i], k4y[i].bounds(), 1.0, 1.0,
                             0.0);
                }

                for (int i = 0; i < num_blocks; ++i) {
                        pre->gather_data(k4x[i], k4y[i], i);
                }

                // k5 = dt * f(t + 1.5 * dt, v + k4)
                for (int i = 0; i < num_blocks; ++i) {
                        pre->update(k5[i], k4x[i], k4y[i], i, t + 1.5 * dt, dt,
                                    0.0);
                }

                //    u1 = u + 22/24 * k1 + 1/24 * k3 + 1/24 * k5
                for (int i = 0; i < num_blocks; ++i) {
                        expr(p[i], k1[i], k3[i], p[i].bounds(), 22.0 / 24,
                             1.0 / 24, 1.0);
                        expr(p[i], k5[i], p[i].bounds(), 1.0 / 24, 1.0);
                }
        }
};

template <typename Tv>
class RK4StaggeredVelocity {
        Pressure<Tv> *pre;
        Velocity<Tv> *vel;
        ComputeArray k1x[NUM_BLOCKS];
        ComputeArray k1y[NUM_BLOCKS];
        ComputeArray k2[NUM_BLOCKS];
        ComputeArray k3x[NUM_BLOCKS];
        ComputeArray k3y[NUM_BLOCKS];
        ComputeArray k4[NUM_BLOCKS];
        ComputeArray k5x[NUM_BLOCKS];
        ComputeArray k5y[NUM_BLOCKS];

       public:
        RK4StaggeredVelocity(Pressure<Tv> &pre, Velocity<Tv> &vel)
            : pre(&pre), vel(&vel) {
                for (int i = 0; i < num_blocks; ++i) {
                        k1x[i] = ComputeArray(vx[i]);
                        k1y[i] = ComputeArray(vy[i]);
                        k2[i] = ComputeArray(p[i]);
                        k3x[i] = ComputeArray(vx[i]);
                        k3y[i] = ComputeArray(vy[i]);
                        k4[i] = ComputeArray(p[i]);
                        k5x[i] = ComputeArray(vx[i]);
                        k5y[i] = ComputeArray(vy[i]);
                }
        }

        void operator()(Tv t, Tv dt) {
                for (int i = 0; i < num_blocks; ++i) {
                        vel->gather_data(p[i], i);
                        pre->gather_data(vx[i], vy[i], i);
                }

                //    k1 = dt * g(t + dt, u1)
                for (int i = 0; i < num_blocks; ++i) {
                        vel->update(k1x[i], k1y[i], p[i], i, t + dt, dt, 0.0);
                }
                //    k2 = dt * f(t + 0.5 * dt, v)
                for (int i = 0; i < num_blocks; ++i) {
                        pre->update(k2[i], vx[i], vy[i], i, t + 0.5 * dt, dt,
                                    0.0);
                }

                // k2 = u1 - k2
                for (int i = 0; i < num_blocks; ++i) {
                        expr(k2[i], p[i], k2[i].bounds(), 1.0, -1.0);
                }

                for (int i = 0; i < num_blocks; ++i) {
                        vel->gather_data(k2[i], i);
                }

                //    k3 = dt * g(t, u1 - k2)
                for (int i = 0; i < num_blocks; ++i) {
                        vel->update(k3x[i], k3y[i], k2[i], i, t, dt, 0.0);
                }

                // k5 = v + k1
                for (int i = 0; i < num_blocks; ++i) {
                        expr(k5x[i], vx[i], k1x[i], k5x[i].bounds(), 1.0, 1.0);
                        expr(k5y[i], vy[i], k1y[i], k5y[i].bounds(), 1.0, 1.0);
                }

                for (int i = 0; i < num_blocks; ++i) {
                        pre->gather_data(k5x[i], k5y[i], i);
                }

                //    k4 = dt * f(t + 1.5 * dt, v + k1)
                for (int i = 0; i < num_blocks; ++i) {
                        pre->update(k4[i], k5x[i], k5y[i], i, t + 1.5 * dt, dt,
                                    0.0);
                }

                // k4 = k4 + u1
                for (int i = 0; i < num_blocks; ++i) {
                        expr(k4[i], p[i], k4[i].bounds(), 1.0, 1.0);
                }

                for (int i = 0; i < num_blocks; ++i) {
                        vel->gather_data(k4[i], i);
                }

                //    k5 = dt * g(t + 2 *dt, u1 + k4)
                for (int i = 0; i < num_blocks; ++i) {
                        vel->update(k5x[i], k5y[i], k4[i], i, t + 2.0 * dt, dt,
                                    0.0);
                }

                //    v1 = v + 22/24 * e1 + 1/24 * e3 + 1/24 * e5
                for (int i = 0; i < num_blocks; ++i) {
                        expr(vx[i], k1x[i], k3x[i], vx[i].bounds(), 22.0 / 24,
                             1.0 / 24, 1.0);
                        expr(vx[i], k5x[i], vx[i].bounds(), 1.0 / 24, 1.0);
                        expr(vy[i], k1y[i], k3y[i], vy[i].bounds(), 22.0 / 24,
                             1.0 / 24, 1.0);
                        expr(vy[i], k5y[i], vy[i].bounds(), 1.0 / 24, 1.0);
                }
        }


};

template <typename Tv>
class Recorder {
        int num_steps;
        int stride;
        int frame;
        int len;
        int block;
        Tv *time, *start;
        ComputeArray *in;
        ComputeArray out;

       public:
        Recorder(int num_steps, ComputeArray &in, 
                 int stride = 1, int block=0) :
                num_steps(num_steps), 
                in(&in), 
                stride(stride), 
                block(block)

        {
                len = (num_steps - 1) / stride + 1;
                out = ComputeArray((size_t)len);
                time = (Tv*)malloc(sizeof(Tv) * len);
                frame = 0;
                start = &out.u[0];
        }

        void operator()(const int step, const Tv t, const Tv dt, int final_step=0) {
                if (!final_step) { 
                        if (step % stride != 0) return;
                }
                if (frame >= len) return;

                out.nx = 1;
                out.u = &out.u[frame];
                recv_L[block].interpolate(out, *in);
                out.nx = len;
                time[frame] = t;
                out.u = start;

                frame++;

        }

        void write_solution(const char *name) {
                HostArray<Tprec> host = out.host();
                _write(host(), name);
        }

        void write_time(const char *name) {
                _write(time, name);
        }

        void _write(const Tv *data, const char *name) {
                FILE *fh = fopen(name, "wb");
                if (fh == NULL) return;
                fwrite(data, sizeof(Tv), frame, fh);
                fclose(fh);
                printf("wrote: %s \n", name);

        }


};

template <typename Tv>
class Vtk {

        int stride;
        int frame = 0;
        const char *name;

        public:
        Vtk(const char *name="default", int stride=50) :
                name(name), 
                stride(stride)
        {
        }

        void operator()(const int step, const Tv t, const Tv dt, int final_step = 0) {
                if (stride < 0) return;
                if (final_step || step % stride == 0) {
                     write_fields(name, frame);
                     frame++;
                }
        }

};

template <typename Tv>
class Info {

        int stride;
        int frame = 0;
        const char *name;
        Recorder<Tv> *rec1, *rec2;
        Vtk<Tv> *vtk;
        
        public:

        Info(Recorder<Tv>& rec1, Recorder<Tv>& rec2, Vtk<Tv>& vtk, int stride=50) :
                rec1(&rec1), 
                rec2(&rec2), 
                vtk(&vtk), 
                stride(stride)
        {
        }

        void operator()(const int step, const Tv t, const Tv dt, int final_step = 0) {
                (*rec1)(step, t, dt, final_step);
                (*rec2)(step, t, dt, final_step);
                (*vtk)(step, t, dt, final_step);
                if (final_step || step % stride == 0) {
                     printf("step = %d \t t = %g \n", step, t);
                     fflush(stdout);
                }
        }
};

//
//
//
