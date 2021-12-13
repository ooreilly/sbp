#include <stdio.h>
#include <string>
#include <assert.h>
#include <mpi.h>

#include "structure/array.hpp"
#include "function/vtk.hpp"
#include "function/gridp.hpp"
#include "function/gridm.hpp"
#include "function/zero.hpp"
#include "function/leapfrog.hpp"
#include "function/lsrk4.hpp"
#include "structure/timer.hpp"

#include "acoustic.hpp"

const int LEAP_FROG = 0;
const int STAGGERED_RK4 = 1;
const int TS_RK4 = 2;

int num_steps(Tprec T, Tprec dt) {
        return (int) ( (T - 0.5 * dt) / dt + 1);
}

int main(int argc, char **argv) {
        printf("Acoustic solver\n");
        if (argc != 17) {
                printf(
                    "usage: %s <name> <refinement> <time_stepping_scheme> <nx> "
                    "<ny> <h> <nu> <vtk_stride> <T> <fp> <t0> <sx> <sy> \n",
                    argv[0]);
                printf("Options\n");
                printf(" name:                  Output filename \n");
                printf(" refinement:            Grid refinement number \n");
                printf(" time_stepping_scheme:  0 = leap-frog, 1 = "
                    "staggered RK4, 2 = RK4 \n");
                printf(
                    " coupling method:       2 = energy conserving, 3 = energy "
                    "dissipating, 4 = energy conserving (naive)\n");
                printf(" nx:                    Grid points x \n");
                printf(" ny:                    Grid points y \n");
                printf(" Lx:                    Domain size x \n");
                printf(" Ly:                    Domain size y \n");
                printf(
                    " nu:                    Time step scaling factor, dt = nu "
                    "* h\n");
                printf(
                    " vtk_stride:            VTK output stride factor. \n "
                    "                       Output every `vtk_stride` time "
                    "steps (also influenced by refinement)\n");
                printf(
                    "                        Use negative number to disable\n");
                printf(
                    " stride:            Time series output stride factor. \n ");
                printf(
                    " T:            Simulation end time. \n ");
                printf(
                    " fp:            Source central frequency. \n ");
                printf(
                    " t0:            Source start time. \n ");
                printf(
                    " sx:            Source x-coordinate. \n ");
                printf(
                    " sy:            Source y-coordinate. \n ");
                exit(1);
        }


        std::string name = std::string(argv[1]);
        int ref = atoi(argv[2]);
        int time_stepping_scheme = atoi(argv[3]);
        int coupling = atoi(argv[4]);
        Tprec nx_ = atof(argv[5]);
        Tprec ny_ = atof(argv[6]);
        Tprec Lx = atof(argv[7]);
        Tprec Ly = atof(argv[8]);
        Tprec nu = atof(argv[9]);
        int vtk_stride = atoi(argv[10]);
        int stride = atoi(argv[11]);
        Tprec T = atof(argv[12]);
        Tprec src_fp = atof(argv[13]);
        Tprec src_t0 = atof(argv[14]);
        Tprec src_sx = atof(argv[15]);
        Tprec src_sy = atof(argv[16]);

        size_t nx[] = {(size_t) ((0.5 * nx_ - 2) * pow(2, ref) + 2),(size_t) ((nx_ - 2) * pow(2, ref) + 2)};
        size_t ny[] = {(size_t) ((ny_ - 2) * pow(2, ref) + 2),(size_t) ((ny_ - 2) * pow(2, ref) + 2)};

        // Set source and receiver properties
        fp[0][0] = src_fp;
        t0[0][0] = src_t0;
        sx[0][0] = src_sx;
        sy[0][0] = src_sy;
        // put receivers on the interface
        rx[0] = 0.5;
        ry[0] = 1.0;
        rx[1] = 0.0;
        ry[1] = 1.0;

        Tprec hx[] = {(Tprec)0.5 * Lx / (nx[0] - 2), (Tprec)Lx / (nx[1] - 2)};
        Tprec hy[] = {(Tprec)Ly / (ny[0] - 2), (Tprec)Ly / (ny[1] - 2)};
        Tprec dt = nu * hx[1] / max(c[0], c[1]);
        
        printf("Coupling type: %d \n", coupling);
        for (int i = 0; i < 2; ++i) {
        printf("grid size %d : %ld %ld, spacing : %f %f  \n", i, nx[i], ny[i],
                        hx[i], hy[i]);
        }
        
        printf("Time step: %g \n", dt);
        printf("Number of steps: %d \n", num_steps(T, dt));

        init(nx, ny, hx, hy, coupling);
        if (vtk_stride > 0)
                init_grids(nx, ny, hx, hy, 0.5 * Lx, Ly);



        Pressure<Tprec> pre;
        Velocity<Tprec> vel;
        RK4<Tprec> rk4(pre, vel);
        RK4StaggeredPressure<Tprec> rk4p(pre, vel);
        RK4StaggeredVelocity<Tprec> rk4v(pre, vel);
        Recorder<Tprec> rec0(num_steps(T, dt), p[0], stride * pow(2, ref), 0);
        Recorder<Tprec> rec1(num_steps(T, dt), p[1], stride * pow(2, ref), 1);
        Vtk<Tprec> vtk(name.c_str(), (Tprec) vtk_stride * pow(2, ref));
        Info<Tprec> info(rec0, rec1, vtk);

        size_t  cmemfree, cmemtotal;
        cudaMemGetInfo(&cmemfree, &cmemtotal);
        printf("CUDA memory: free = %f GB \ttotal = %f GB \n",
                   cmemfree/ 1e9, cmemtotal / 1e9);


        HostTimer timer;
        timer.start();
        switch (time_stepping_scheme) {
                case LEAP_FROG:
                        leapfrog(pre, vel, info, T, dt);
                        break;
                case STAGGERED_RK4:
                        leapfrog(rk4p, rk4v, info, T, dt);
                        break;
                case TS_RK4:
                        lsrk4(rk4, info, T, dt);
                        break;
        }
        timer.stop();
        printf("Simulation took: %g s \n", timer.elapsed() / 1000.0);
        fflush(stdout);


        rec0.write_solution((name + std::string("_p0.bin")).c_str());
        rec1.write_solution((name + std::string("_p1.bin")).c_str());
        rec0.write_time((name + std::string("_t.bin")).c_str());
}
