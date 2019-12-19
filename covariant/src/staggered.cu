#define RSTRCT __restrict__

const static int stx = 32;
const static int sty = 1;

const static int snb = 8;

__launch_bounds__(32)
__global__ void st_rates_int(Tv RSTRCT *dp, 
                Tv RSTRCT *dv1, 
                Tv RSTRCT *dv2, 
                const Tv RSTRCT *p, 
                const Tv RSTRCT *v1, 
                const Tv RSTRCT *v2, 
                const Tv RSTRCT *Jp, 
                const Tv RSTRCT *J1, 
                const Tv RSTRCT *J2, 
                const Tv RSTRCT *g11, 
                const Tv RSTRCT *g12, 
                const Tv RSTRCT *g22, 
                const Ti nx, const Ti ny, const
                Ti my, const Tv hi1, const Tv hi2, const int k)
{

        Tv i = threadIdx.x + blockDim.x * blockIdx.x;
        Tv j = threadIdx.y + blockDim.y * blockIdx.y;

        int p0p0 = i + j * my;

        int m3p0 = p0p0 - 3 * my;
        int m2p0 = p0p0 - 2 * my;
        int m1p0 = p0p0 - 1 * my;
        int p1p0 = p0p0 + 1 * my;
        int p2p0 = p0p0 + 2 * my;

        int p0m2 = p0p0 - 2;
        int p0m1 = p0p0 - 1;
        int p0p1 = p0p0 + 1;
        int p0p2 = p0p0 + 2;

        int m1m3 = p0p0 - my - 3;
        int m1m2 = p0p0 - my - 2;
        int m1m1 = p0p0 - my - 1;
        int m1p1 = p0p0 - my + 1;
        int m1p2 = p0p0 - my + 2;
        int m1p3 = p0p0 - my + 3;

        int p0m3 = p0p0 - 3;
        int p0p3 = p0p0 + 3;

        int p1m3 = p0p0 + my - 3;
        int p1m2 = p0p0 + my - 2;
        int p1m1 = p0p0 + my - 1;
        int p1p1 = p0p0 + my + 1;
        int p1p2 = p0p0 + my + 2;
        int p1p3 = p0p0 + my + 3;

        int p2m3 = p0p0 + 2 * my - 3;
        int p2m2 = p0p0 + 2 * my - 2;
        int p2m1 = p0p0 + 2 * my - 1;
        int p2p1 = p0p0 + 2 * my + 1;
        int p2p2 = p0p0 + 2 * my + 2;
        int p2p3 = p0p0 + 2 * my + 3;

        int m3m1 = p0p0 - 3 * my - 1;

        int m2m1 = p0p0 - 2 * my - 1;
        int m2p1 = p0p0 - 2 * my + 1;
        int m2p2 = p0p0 - 2 * my + 2;
        int p3p0 = p0p0 + 3 * my;
        int p3m1 = p0p0 + 3 * my - 1;

        int m3p1 = p0p0 - 3 * my + 1;
        int m3p2 = p0p0 - 3 * my + 2;
        int p3p1 = p0p0 + 3 * my + 1;
        int p3p2 = p0p0 + 3 * my + 2;


        // Load 
        Tv dp_p0p0 = dp[p0p0];
        Tv dv1_p0p0 = dv1[p0p0];
        Tv dv2_p0p0 = dv2[p0p0];

        // p grid

          // Lines needed for treating g12 * dp/dy
          Tv py_m2m3 = p[m1m3];
          Tv py_m2m2 = p[m1m2];
          Tv py_m2m1 = p[m1m1];
          Tv py_m2p1 = p[m1p1];
          Tv py_m2p2 = p[m1p2];
          Tv py_m2p3 = p[m1p3];

          Tv py_m1m3 = p[p0m3];
          Tv py_m1m2 = p[p0m2];
          Tv py_m1m1 = p[p0m1];
          Tv py_m1p1 = p[p0p1];
          Tv py_m1p2 = p[p0p2];
          Tv py_m1p3 = p[p0p3];

          Tv py_p1m3 = p[p1m3];
          Tv py_p1m2 = p[p1m2];
          Tv py_p1m1 = p[p1m1];
          Tv py_p1p1 = p[p1p1];
          Tv py_p1p2 = p[p1p2];
          Tv py_p1p3 = p[p1p3];

          Tv py_p2m3 = p[p2m3];
          Tv py_p2m2 = p[p2m2];
          Tv py_p2m1 = p[p2m1];
          Tv py_p2p1 = p[p2p1];
          Tv py_p2p2 = p[p2p2];
          Tv py_p2p3 = p[p2p3];

          // Lines needed for treating g12 * dp/dx
          Tv px_m3m2 = p[m3m1];
          Tv px_m2m2 = p[m2m1];
          Tv px_m1m2 = p[m1m1];
          Tv px_p1m2 = p[p1m1];
          Tv px_p2m2 = p[p2m1];
          Tv px_p3m2 = p[p3m1];

          Tv px_m3m1 = p[m3p0];
          Tv px_m2m1 = p[m2p0];
          Tv px_m1m1 = p[m1p0];
          Tv px_p1m1 = p[p1p0];
          Tv px_p2m1 = p[p2p0];
          Tv px_p3m1 = p[p3p0];
          
          Tv px_m3p1 = p[m3p1];
          Tv px_m2p1 = p[m2p1];
          Tv px_m1p1 = p[m1p1];
          Tv px_p1p1 = p[p1p1];
          Tv px_p2p1 = p[p2p1];
          Tv px_p3p1 = p[p3p1];
          
          Tv px_m3p2 = p[m3p2];
          Tv px_m2p2 = p[m2p2];
          Tv px_m1p2 = p[m1p2];
          Tv px_p1p2 = p[p1p2];
          Tv px_p2p2 = p[p2p2];
          Tv px_p3p2 = p[p3p2];


        Tv p_m2p0 = p[m1p0];
        Tv p_m1p0 = p[p0p0];
        Tv p_p1p0 = p[p1p0];
        Tv p_p2p0 = p[p2p0];

        Tv p_p0m2 = p[p0m1];
        Tv p_p0m1 = p[p0p0];
        Tv p_p0p1 = p[p0p1];
        Tv p_p0p2 = p[p0p2];

        Tv Jp_m2p0 = Jp[m1p0];
        Tv Jp_m1p0 = Jp[p0p0];
        Tv Jp_p1p0 = Jp[p1p0];
        Tv Jp_p2p0 = Jp[p2p0];

        Tv g12_m2p0 = g12[m1p0];
        Tv g12_m1p0 = g12[p0p0];
        Tv g12_p1p0 = g12[p1p0];
        Tv g12_p2p0 = g12[p2p0];

        Tv Jp_p0m2 = Jp[p0m1];
        Tv Jp_p0m1 = Jp[p0p0];
        Tv Jp_p0p1 = Jp[p0p1];
        Tv Jp_p0p2 = Jp[p0p2];

        Tv g12_p0m2 = g12[p0m1];
        Tv g12_p0m1 = g12[p0p0];
        Tv g12_p0p1 = g12[p0p1];
        Tv g12_p0p2 = g12[p0p2];
        
        Tv Jip = 1.0 / Jp_p0m1;

        // v1 grid
        Tv v1_m2p0 = v1[m2p0];
        Tv v1_m1p0 = v1[m1p0];
        Tv v1_p1p0 = v1[p0p0];
        Tv v1_p2p0 = v1[p1p0];

        Tv J1_m2p0 = J1[m2p0];
        Tv J1_m1p0 = J1[m1p0];
        Tv J1_p1p0 = J1[p0p0];
        Tv J1_p2p0 = J1[p1p0];

        Tv J1i = 1.0 / Jp_p1p0;
        Tv g11_p0p0 = g11[p0p0];

        // v2 grid

        Tv v2_p0m2 = v2[p0m2];
        Tv v2_p0m1 = v2[p0m1];
        Tv v2_p0p1 = v2[p0p0];
        Tv v2_p0p2 = v2[p0p1];

        Tv J2_p0m2 = J2[p0m2];
        Tv J2_p0m1 = J2[p0m1];
        Tv J2_p0p1 = J2[p0p0];
        Tv J2_p0p2 = J2[p0p1];

        Tv J2i = 1.0 / Jp_p0p1;
        Tv g22_p0p0 = g22[p0p0];


        // 4pt Staggered first derivative
        Tv s2 = -0.041666666666666664;
        Tv s1 = 1.125;

        // 7pt Central first derivative
        Tv c3 = 0.0026041666666666665;
        Tv c2 = -0.09375;
        Tv c1 = 0.6796875;

        // 4pt interpolation
        Tv p2 = -1.0/16;
        Tv p1 = 9.0/16;

        // Store
        if (i >= snb && i < ny - snb && j >= snb && j < nx - snb) {

        // Pressure equation
        Tv dp_rhs = ak[k] * dp_p0p0 +
                   - Jip * hi1 *
                        (s2 * (J1_p2p0 * v1_p2p0 - J1_m2p0 * v1_m2p0) +
                         s1 * (J1_p1p0 * v1_p1p0 - J1_m1p0 * v1_m1p0))
                  -  Jip * hi2 *
                        (s2 * (J2_p0p2 * v2_p0p2 - J2_p0m2 * v2_p0m2) 
                       + s1 * (J2_p0p1 * v2_p0p1 - J2_p0m1 * v2_p0m1));


        // Velocity equation 1 
        Tv dpdr1 = hi1 * (s2 * (p_p2p0 - p_m2p0) +
                          s1 * (p_p1p0 - p_m1p0));
        Tv lxm2 = c3 * (py_m2p3 - py_m2m3) +
                  c2 * (py_m2p2 - py_m2m2) +
                  c1 * (py_m2p1 - py_m2m1);
        Tv lxm1 = c3 * (py_m1p3 - py_m1m3) +
                  c2 * (py_m1p2 - py_m1m2) +
                  c1 * (py_m1p1 - py_m1m1);
        Tv lxp1 = c3 * (py_p1p3 - py_p1m3) +
                  c2 * (py_p1p2 - py_p1m2) +
                  c1 * (py_p1p1 - py_p1m1);
        Tv lxp2 = c3 * (py_p2p3 - py_p2m3) +
                  c2 * (py_p2p2 - py_p2m2) +
                  c1 * (py_p2p1 - py_p2m1);

        Tv dpdr12 = hi2 * (p2 * (Jp_p2p0 * g12_p2p0 * lxp2 +
                                 Jp_m2p0 * g12_m2p0 * lxm2) +
                           p1 * (Jp_p1p0 * g12_p1p0 * lxp1 +
                                 Jp_m1p0 * g12_m1p0 * lxm1));


        Tv dv1_rhs = ak[k] * dv1_p0p0 - g11_p0p0 * dpdr1 - J1i * dpdr12;

        //Velocity equation 2
        Tv dpdr2 = hi2 * (s2 * (p_p0p2 - p_p0m2) +
                          s1 * (p_p0p1 - p_p0m1));
        Tv lym2 = c3 * (px_p3m2 - px_m3m2) +
                  c2 * (px_p2m2 - px_m2m2) +
                  c1 * (px_p1m2 - px_m1m2);
        Tv lym1 = c3 * (px_p3m1 - px_m3m1) +
                  c2 * (px_p2m1 - px_m2m1) +
                  c1 * (px_p1m1 - px_m1m1);
        Tv lyp1 = c3 * (px_p3p1 - px_m3p1) +
                  c2 * (px_p2p1 - px_m2p1) +
                  c1 * (px_p1p1 - px_m1p1);
        Tv lyp2 = c3 * (px_p3p2 - px_m3p2) +
                  c2 * (px_p2p2 - px_m2p2) +
                  c1 * (px_p1p2 - px_m1p2);

        Tv dpdr21 = hi1 * (p2 * (Jp_p0p2 * g12_p0p2 * lyp2 +
                                 Jp_p0m2 * g12_p0m2 * lym2) +
                           p1 * (Jp_p0p1 * g12_p0p1 * lyp1 +
                                 Jp_p0m1 * g12_p0m1 * lym1));

        Tv dv2_rhs = ak[k] * dv2_p0p0 - g22_p0p0 * dpdr2 - J2i * dpdr21;

                dp[p0p0] = dp_rhs;
                dv1[p0p0] = dv1_rhs;
                dv2[p0p0] = dv2_rhs;
        }
}

__global__ void st_periodic_x(Tv RSTRCT *p, Tv RSTRCT *v1, Tv RSTRCT *v2,
                             const Ti nx, const Ti ny,
                             const Ti my) {


        // top and bottom boundary
        Tv i = threadIdx.x + blockDim.x * blockIdx.x;
        Tv j = threadIdx.y + blockDim.y * blockIdx.y;
        int p0p0 = j + i * my;
        int p0p1 = p0p0 + snb;
        
        int p0m1 = p0p0 + ny - snb;
        int p0m2 = p0m1 - snb;


        if (j < snb && i < nx) {
                p[p0p0] = p[p0m2];
                p[p0m1] = p[p0p1];
        }

}

__global__ void st_periodic_y(Tv RSTRCT *p, Tv RSTRCT *v1, Tv RSTRCT *v2,
                             const Ti nx, const Ti ny,
                             const Ti my) {


        // left and right boundary
        Tv i = threadIdx.x + blockDim.x * blockIdx.x;
        Tv j = threadIdx.y + blockDim.y * blockIdx.y;
        int p0p0 = i + j * my;
        int p1p0 = p0p0 + snb * my;
        
        int m1p0 = p0p0 + my * nx - snb * my;
        int m2p0 = m1p0 - snb * my;


        if (j < snb && i >= snb && i < ny - snb) {
                p[p0p0] = p[m2p0];
        }
        if (j < snb && i >= snb && i < ny - snb) {
                p[m1p0] = p[p1p0];
        }

}

__global__ void st_update_int(Tv RSTRCT *p, Tv RSTRCT *v1, Tv RSTRCT *v2,
                               const Tv RSTRCT *dp, const Tv RSTRCT *dv1,
                               const Tv RSTRCT *dv2, const Ti nx, const Ti ny,
                               const Ti my, const Tv dt, const int k) {

        Tv i = threadIdx.x + blockDim.x * blockIdx.x;
        Tv j = threadIdx.y + blockDim.y * blockIdx.y;

        int p0p0 = i + j * my;

        if (i >= snb && i < ny - snb && j >= snb && j < nx - snb) {
        p[p0p0]  = p[p0p0]  + bk[k] * dt  * dp[p0p0];   
        v1[p0p0] = v1[p0p0] + bk[k] * dt * dv1[p0p0];   
        v2[p0p0] = v2[p0p0] + bk[k] * dt * dv2[p0p0];   
        }

}

void st_rates(dArray<Tv>& dp, dArray<Tv>& dv1, dArray<Tv>& dv2,  
               dArray<Tv>& p, dArray<Tv>& v1, dArray<Tv>& v2, 
               dArray<Tv>& J, 
               dArray<Tv>& J1, 
               dArray<Tv>& J2, 
               dArray<Tv>& g11, 
               dArray<Tv>& g12, 
               dArray<Tv>& g22, 
               const Ti nx,
               const Ti ny, const Ti my, const Tv hi1, const Tv hi2, const int k)
{
        dim3 threads (stx, sty, 1);
        dim3 blocks ((ny - 1) / stx + 1, (nx - 1) / sty + 1, 1);
        st_rates_int<<<blocks, threads>>>(dp.x, dv1.x, dv2.x, p.x, v1.x, v2.x,
                                           J.x, J1.x, J2.x, g11.x, g12.x, g22.x, nx, ny, my,
                                           hi1, hi2, k);
}

void st_update(
 dArray<Tv>& p, dArray<Tv>& v1, dArray<Tv>& v2,        
 dArray<Tv>& dp, dArray<Tv>& dv1, dArray<Tv>& dv2, 
 const Ti nx, const Ti ny, const Ti my, const Tv dt, const Ti k)
{
        dim3 threads (stx, sty, 1);
        dim3 blocks ((ny - 1) / stx + 1, (nx - 1) / ty + 1, 1);
        st_update_int<<<blocks, threads>>>(p.x, v1.x, v2.x, dp.x, dv1.x, dv2.x,
                                           nx, ny, my, dt, k);

}
 void st_periodic(
 dArray<Tv>& p, dArray<Tv>& v1, dArray<Tv>& v2,        
 const Ti nx, const Ti ny, const Ti my)
{
        {
        dim3 threads (stx, 1, 1);
        dim3 blocks ((nx - 1) / stx + 1, nb, 1);
        st_periodic_x<<<blocks, threads>>>(p.x, v1.x, v2.x, nx, ny, my);
        }
        {
        dim3 threads (stx, 1, 1);
        dim3 blocks ((ny - 1) / stx + 1, nb, 1);
        st_periodic_y<<<blocks, threads>>>(p.x, v1.x, v2.x, nx, ny, my);
        }
}
#undef RSTRCT
