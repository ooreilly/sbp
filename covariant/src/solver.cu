#include <stdio.h>
#include <hlmv.hpp>
#include <dlmv.hpp>
#include <hArray.hpp>
#include <dArray.hpp>
#include <dcsr.hpp>
#include <hcsr.hpp>
#include <dict.hpp>
#include <cusparse.h>
#include <string>

typedef double Tv;
typedef int Ti;

#include "helper.hpp"

int main(int argc, char **argv) {

        cusparseHandle_t cusparseH;
        cublasHandle_t cublasH;
        cusparseErrCheck(cusparseCreate(&cusparseH));
        cublasErrCheck(cublasCreate(&cublasH));

        if (argc != 2) {
                printf("usage: %s <project_dir> \n", argv[0]);
                return -1;
        }
        std::string project_dir = std::string(argv[1]) + "/";

        printf("Project directory: %s \n", project_dir.c_str());
        Dict cfg = read(project_dir + "config.txt");
        int nt = cfg["nt"];
        Tv dt = cfg["dt"];
        int ninfo = cfg["ninfo"];
        int nxp = cfg["nxp"];
        int nyp = cfg["nyp"];
        int ndump = cfg["ndump"];
        int nvtk = cfg["nvtk"];
        int ncsv = cfg["ncsv"];
        int nrecv = cfg["nrecv"];

        std::string output_dir = cfg["outputdir"];
        output_dir += "/";
        dump(cfg);

        // host variables are denote by `_` at the end. 
        // e.g. `u_` resides on the host, whereas `u` resides on the gpu.
        hCSR<Tv, Ti> A_ = read(project_dir + "A.bin");
        hArray<Tv> u_ = read(project_dir + "q.bin");
        hArray<Tv> d_ = read(project_dir + "d.bin");
        hArray<Tv> xp_ = read(project_dir + "xp.bin");
        hArray<Tv> yp_ = read(project_dir + "yp.bin");
        hArray<Tv> r_ = read(project_dir + "r.bin");

        dCSR<Tv, Ti> A = htod(A_);
        dArray<Tv> u = htod(u_);
        dArray<Tv> d = htod(d_);
        dArray<Tv> du(u.size);
        dArray<Tv> r = htod(r_);

        // Output to csv file
        hArray<Tv> t_csv_(size_t(nt / ncsv)); 
        hArray<Tv> u_csv_(t_csv_.size * nrecv);
        dArray<Tv> u_csv(nrecv * t_csv_.size);

        cudaEvent_t start, stop;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
        Tv t = 0.0;
        int dump_step = 0;
        int csv_step = 0;
        int vtk_step = 0;
        float elapsed = 0;
        cudaEventRecord(start);
        printf("step \t time \t elapsed (s) \t time per step (ms) \n");
        for (int step = 0; step < nt; ++step) {

                if ( ninfo > 0 && step % ninfo == 0) {
                        cudaEventRecord(stop);
                        cudaEventSynchronize(stop);
                        cudaEventElapsedTime(&elapsed, start, stop);
                        printf("%3d \t %3.4f \t %3.4f \t %3.4f \n", step, t,
                               elapsed * 1e-3, elapsed / step);
                }

                if ( ndump > 0 && step % ndump == 0) {
                        dtoh(u_, u);
                        std::string out = output_dir + "u_" +
                                          std::to_string(dump_step) + ".bin";
                        write(out, u_);
                        dump_step++;
                }

                if ( ncsv > 0 && step % ncsv == 0) {
                        // Extract values to save time series for
                        t_csv_[csv_step] = t;
                        u_csv_[csv_step] = dot(cublasH, r, u);
                        csv_step++;
                }

                if ( nvtk > 0 && step % nvtk == 0) {
                        dtoh(u_, u);
                        std::string out = output_dir + "u_" +
                                          std::to_string(vtk_step) + ".vtk";
                        write_vtk(out, u_.x, xp_, yp_, nxp, nyp);
                        vtk_step++;
                }

                for (int k = 0; k < rk4_n; ++k) {
                        t = step * dt + c[k]*dt;

                        // du = A*u + a[k]*du
                        mv(cusparseH, du, A, u, 1.0, a[k]);

                        // du = du + d * g(t + c[k]*dt)
                        Tv gval = ricker(t);
                        axpy(cublasH, du, d, gval);

                        // u = u + dt*b[k]*du;
                        axpy(cublasH, u, du, dt * b[k]);
                }
                t = (1 + step) * dt;
        }

        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&elapsed, start, stop);

        // Write time-series to csv-file
        {
                std::string out = output_dir + "u.csv";
                write_csv(out, u_csv_, t_csv_);
        }

        cudaEventElapsedTime(&elapsed, start, stop);
        printf("Total compute time: %g s \n", elapsed * 1e-3);

        // Write log file
        {
        std::string out = output_dir + "log.txt";
        FILE *fh = fopen(out.c_str(), "w");
        fprintf(fh, "m=%ld\n", A.m);
        fprintf(fh, "n=%ld\n", A.n);
        fprintf(fh, "nnz=%ld\n", A.nnz);
        fprintf(fh, "elapsed=%g\n", elapsed * 1e-3);
        fprintf(fh, "nt=%d\n", nt);
        fclose(fh);
        }

        cusparseErrCheck(cusparseDestroy(cusparseH));
        cublasErrCheck(cublasDestroy(cublasH));

        return 0;
}
