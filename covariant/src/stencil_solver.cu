#include <stdio.h>
#include <cuda.h>
#include <hlmv.hpp>
#include <dlmv.hpp>
#include <hArray.hpp>
#include <dArray.hpp>
#include <dcsr.hpp>
#include <hcsr.hpp>
#include <dict.hpp>
#include <cusparse.h>
#include <string>

#include "helper.hpp"
#include "collocated.cu"
#include "staggered.cu"

int main(int argc, char **argv) {


        if (argc != 2) {
                printf("usage: %s <project_dir> \n", argv[0]);
                return -1;
        }

        cublasHandle_t cublasH;
        cublasErrCheck(cublasCreate(&cublasH));


        std::string project_dir = std::string(argv[1]) + "/";

        printf("Project directory: %s \n", project_dir.c_str());
        Dict cfg = read(project_dir + "config.txt");
        int nt = 2000;
        int nx = cfg["nx"];
        int ny = cfg["ny"];
        int scheme = cfg["scheme"];
        int my = ny;
        int ninfo = 100;
        int nvtk = 100;
        Tv hi1 = cfg["hi1"];
        Tv hi2 = cfg["hi2"];
        Tv dt = cfg["dt"];
        dump(cfg);

        hArray<Tv> p_ = read(project_dir + "p.bin");
        hArray<Tv> v1_ = read(project_dir + "v1.bin");
        hArray<Tv> v2_ = read(project_dir + "v2.bin");

        // Metrics, collocated
        hArray<Tv> J_ = read(project_dir + "J.bin");
        hArray<Tv> g11_ = read(project_dir + "g11.bin");
        hArray<Tv> g12_ = read(project_dir + "g12.bin");
        hArray<Tv> g22_ = read(project_dir + "g22.bin");

        // Metrics, staggered
        hArray<Tv> Jp_ = read(project_dir + "Jp.bin");
        hArray<Tv> J1_ = read(project_dir + "J1.bin");
        hArray<Tv> J2_ = read(project_dir + "J2.bin");
        hArray<Tv> g1_11_ = read(project_dir + "g1_11.bin");
        hArray<Tv> gp_12_ = read(project_dir + "gp_12.bin");
        hArray<Tv> g2_22_ = read(project_dir + "g2_22.bin");

        hArray<Tv> xp_ = read(project_dir + "xp.bin");
        hArray<Tv> yp_ = read(project_dir + "yp.bin");
        hArray<Tv> d_ = read(project_dir + "d.bin");
        
        //Fields
        dArray<Tv> p = htod(p_);
        dArray<Tv> v1 = htod(v1_);
        dArray<Tv> v2 = htod(v2_);

        ////Rates
        dArray<Tv> dp(p.size);
        dArray<Tv> dv1(v1_.size);
        dArray<Tv> dv2(v2_.size);

        // Collocated Metrics
        dArray<Tv> J = htod(J_);
        dArray<Tv> g11 = htod(g11_);
        dArray<Tv> g12 = htod(g12_);
        dArray<Tv> g22 = htod(g22_);

        // Staggered Metrics
        dArray<Tv> Jp = htod(Jp_);
        dArray<Tv> J1 = htod(J1_);
        dArray<Tv> J2 = htod(J2_);
        dArray<Tv> g1_11 = htod(g1_11_);
        dArray<Tv> gp_12 = htod(gp_12_);
        dArray<Tv> g2_22 = htod(g2_22_);

        dArray<Tv> d = htod(d_);

        col_init(a, b);

        Tv t = 0.0;
        cudaEvent_t start, stop;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);

        int vtk_step = 0;
        float elapsed = 0;
        cudaEventRecord(start);
        printf("step \t time \t elapsed (s) \t time per step (ms) \n");
        int run = 1;

        if (run) {
        col_periodic(J, g11, g22, nx, ny, my);
        for (int step = 0; step < nt; ++step) {

                if ( ninfo > 0 && step % ninfo == 0) {
                        cudaEventRecord(stop);
                        cudaEventSynchronize(stop);
                        cudaEventElapsedTime(&elapsed, start, stop);
                        printf("%3d \t %3.4f \t %3.4f \t %3.4f \n", step, t,
                               elapsed * 1e-3, elapsed / step);
                }

                if ( nvtk > 0 && step % nvtk == 0) {
                        dtoh(p_, p);
                        std::string out = project_dir + "p_" +
                                          std::to_string(vtk_step) + ".vtk";
                        write_vtk(out, p_.x, xp_, yp_, nx, ny);
                        vtk_step++;
                }

                for (int k = 0; k < rk4_n; ++k) {
                        t = step * dt + c[k]*dt;

                        if (scheme == 0) {
                                col_periodic(p, v1, v2, nx, ny, my);
                                col_rates(dp, dv1, dv2, p, v1, v2, J, g11, g12,
                                          g22, nx, ny, my, hi1, hi2, k);
                        } else if (scheme == 1) {
                                st_periodic(p, v1, v2, nx, ny, my);
                                st_rates(dp, dv1, dv2, p, v1, v2, Jp, J1, J2,
                                                g1_11, gp_12,
                                         g2_22, nx, ny, my, hi1, hi2, k);

                        }

                        // dp = dp + d * g(t + c[k]*dt)
                        Tv gval = ricker(t);
                        axpy(cublasH, dp, d, gval);

                        if (scheme == 0) {
                                col_update(p, v1, v2, dp, dv1, dv2, nx, ny, my,
                                           dt, k);
                        } else if (scheme == 1) {
                                st_update(p, v1, v2, dp, dv1, dv2, nx, ny, my,
                                          dt, k);
                        }
                }
                t = (1 + step) * dt;

        }
        }  else {
        //dump2d(Jp_, nx, ny);
        //dump2d(J1_, nx, ny);
        //dump2d(J2_, nx, ny);
        //dump2d(g1_11_, nx, ny);
        //dump2d(gp_12_, nx, ny);
        //dump2d(g2_22_, nx, ny);
        dump2d(Jp_, nx, ny);
        dump2d(gp_12_, nx, ny);
        dump2d(p_, nx, ny);
        //Dump2d(v1_, nx, ny);
        //Dump2d(v2_, nx, ny);

        //col_periodic(p, v1, v2, nx, ny, my);
        st_rates(dp, dv1, dv2, p, v1, v2, Jp, J1, J2, g1_11, gp_12, g2_22, nx,
                                  ny, my, hi1, hi2, 0);

        dtoh(p_, dp);
        dtoh(v1_, dv1);
        dtoh(v2_, dv2);
        dump2d(v1_, nx, ny);
        dump2d(v2_, nx, ny);


        }
        //dump2d(v1_, nx, ny);
        //dump2d(v2_, nx, ny);

        //dump2d(g11_, nx, ny);
        //dump2d(g12_, nx, ny);
        //dump2d(g22_, nx, ny);

        cublasErrCheck(cublasDestroy(cublasH));

        return 0;
}
