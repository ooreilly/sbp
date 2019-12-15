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

int main(int argc, char **argv) {

        if (argc != 2) {
                printf("usage: %s <project_dir> \n", argv[0]);
                return -1;
        }
        std::string project_dir = std::string(argv[1]) + "/";

        printf("Project directory: %s \n", project_dir.c_str());
        Dict cfg = read(project_dir + "config.txt");
        int nx = cfg["nx"];
        int ny = cfg["ny"];
        int my = ny;
        Tv hi1 = cfg["hi1"];
        Tv hi2 = cfg["hi2"];
        dump(cfg);

        hArray<Tv> p_ = read(project_dir + "p.bin");
        hArray<Tv> v1_ = read(project_dir + "v1.bin");
        hArray<Tv> v2_ = read(project_dir + "v2.bin");
        hArray<Tv> J_ = read(project_dir + "J.bin");


        dump2d(p_, nx, ny);
        dump2d(v1_, nx, ny);
        dump2d(v2_, nx, ny);
        
        //Fields
        dArray<Tv> p = htod(p_);
        dArray<Tv> v1 = htod(v1_);
        dArray<Tv> v2 = htod(v2_);

        ////Rates
        dArray<Tv> dp(p.size);
        dArray<Tv> dv1(v1_.size);
        dArray<Tv> dv2(v2_.size);

        // Metrics
        dArray<Tv> J = htod(J_);

        col_init(a, b);
        col_rates(dp, dv1, dv2, p, v1, v2, J, nx, ny, my, hi1, hi2, 0);

        dtoh(p_, dp);
        dtoh(v1_, dv1);
        dtoh(v2_, dv2);

        dump2d(p_, nx, ny);
        dump2d(v1_, nx, ny);
        dump2d(v2_, nx, ny);
        return 0;
}
