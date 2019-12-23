#define RSTRCT __restrict__

const static int tx = 32;
const static int ty = 1;

const static int nb = 8;

static __device__ __constant__ Tv ak[5];
static __device__ __constant__ Tv bk[5];


void col_init(const Tv *ak_, const Tv *bk_)
{
        cudaErrCheck(cudaMemcpyToSymbol(ak, ak_, sizeof(Tv) * 5));
        cudaErrCheck(cudaMemcpyToSymbol(bk, bk_, sizeof(Tv) * 5));
}

__launch_bounds__(32)
__global__ void col_rates_int(Tv RSTRCT *dp, 
                Tv RSTRCT *dv1, 
                Tv RSTRCT *dv2, 
                const Tv RSTRCT *p, 
                const Tv RSTRCT *v1, 
                const Tv RSTRCT *v2, 
                const Tv RSTRCT *J, 
                const Tv RSTRCT *g11, 
                const Tv RSTRCT *g12, 
                const Tv RSTRCT *g22, 
                const Ti nx, const Ti ny, const
                Ti my, const Tv hi1, const Tv hi2, const int k)
{

        Tv i = threadIdx.x + blockDim.x * blockIdx.x;
        Tv j = threadIdx.y + blockDim.y * blockIdx.y;

        int p0p0 = i + j * my;

        // indices in the grid x-direction
        int m2p0 = p0p0 - 2 * my;
        int m1p0 = p0p0 - 1 * my;
        int p1p0 = p0p0 + 1 * my;
        int p2p0 = p0p0 + 2 * my;

        // indices in the grid y-direction
        int p0m2 = p0p0 - 2;
        int p0m1 = p0p0 - 1;
        int p0p1 = p0p0 + 1;
        int p0p2 = p0p0 + 2;


        // Load 
        Tv dp_p0p0 = dp[p0p0];
        Tv dv1_p0p0 = dv1[p0p0];
        Tv dv2_p0p0 = dv2[p0p0];

        Tv p_m2p0 = p[m2p0];
        Tv p_m1p0 = p[m1p0];
        Tv p_p1p0 = p[p1p0];
        Tv p_p2p0 = p[p2p0];

        Tv p_p0m2 = p[p0m2];
        Tv p_p0m1 = p[p0m1];
        Tv p_p0p1 = p[p0p1];
        Tv p_p0p2 = p[p0p2];

        Tv v1_m2p0 = v1[m2p0];
        Tv v1_m1p0 = v1[m1p0];
        Tv v1_p1p0 = v1[p1p0];
        Tv v1_p2p0 = v1[p2p0];

        Tv v2_p0m2 = v2[p0m2];
        Tv v2_p0m1 = v2[p0m1];
        Tv v2_p0p1 = v2[p0p1];
        Tv v2_p0p2 = v2[p0p2];

        Tv J_p2p0 = J[p2p0];
        Tv J_p1p0 = J[p1p0];
        Tv J_p0p0 = J[p0p0];
        Tv J_m1p0 = J[m1p0];
        Tv J_m2p0 = J[m2p0];

        Tv J_p0p2 = J[p0p2];
        Tv J_p0p1 = J[p0p1];
        Tv J_p0m1 = J[p0m1];
        Tv J_p0m2 = J[p0m2];

        Tv g11_p0p0 = g11[p0p0];
        Tv g12_p0p0 = g12[p0p0];
        Tv g22_p0p0 = g22[p0p0];

        Tv Ji = 1.0 / J_p0p0;

        Tv c2 = -0.083333333333333;
        Tv c1 =  0.666666666666667;

        // Store
        if (i >= nb && i < ny - nb && j >= nb && j < nx - nb) {
                dp[p0p0] =
                    ak[k] * dp_p0p0 +
                   - Ji * hi1 *
                        (c2 * (J_p2p0 * v1_p2p0 - J_m2p0 * v1_m2p0) +
                         c1 * (J_p1p0 * v1_p1p0 - J_m1p0 * v1_m1p0)) 
                   -
                    Ji * hi2 *
                        (c2 * (J_p0p2 * v2_p0p2 - J_p0m2 * v2_p0m2) 
                       + c1 * (J_p0p1 * v2_p0p1 - J_p0m1 * v2_p0m1));

                Tv dpdr1 = hi1 * (c2 * (p_p2p0 - p_m2p0) +
                                  c1 * (p_p1p0 - p_m1p0));
                Tv dpdr2 = hi2 * (c2 * (p_p0p2 - p_p0m2) +
                                  c1 * (p_p0p1 - p_p0m1));
                dv1[p0p0] =
                    ak[k] * dv1_p0p0 - g11_p0p0 * dpdr1 - g12_p0p0 * dpdr2;
                dv2[p0p0] =
                    ak[k] * dv2_p0p0 - g12_p0p0 * dpdr1 - g22_p0p0 * dpdr2;
        }
}


__global__ void col_update_int(Tv RSTRCT *p, Tv RSTRCT *v1, Tv RSTRCT *v2,
                               const Tv RSTRCT *dp, const Tv RSTRCT *dv1,
                               const Tv RSTRCT *dv2, const Ti nx, const Ti ny,
                               const Ti my, const Tv dt, const int k) {

        Tv i = threadIdx.x + blockDim.x * blockIdx.x;
        Tv j = threadIdx.y + blockDim.y * blockIdx.y;

        int p0p0 = i + j * my;

        if (i >= nb && i < ny - nb && j >= nb && j < nx - nb) {
        p[p0p0]  = p[p0p0]  + bk[k] * dt  * dp[p0p0];   
        v1[p0p0] = v1[p0p0] + bk[k] * dt * dv1[p0p0];   
        v2[p0p0] = v2[p0p0] + bk[k] * dt * dv2[p0p0];   
        }

}

void col_rates(dArray<Tv>& dp, dArray<Tv>& dv1, dArray<Tv>& dv2,  
               dArray<Tv>& p, dArray<Tv>& v1, dArray<Tv>& v2, 
               dArray<Tv>& J, 
               dArray<Tv>& g11, 
               dArray<Tv>& g12, 
               dArray<Tv>& g22, 
               const Ti nx,
               const Ti ny, const Ti my, const Tv hi1, const Tv hi2, const int k)
{
        dim3 threads (tx, ty, 1);
        dim3 blocks ((ny - 1) / tx + 1, (nx - 1) / ty + 1, 1);
        col_rates_int<<<blocks, threads>>>(dp.x, dv1.x, dv2.x, p.x, v1.x, v2.x,
                                           J.x, g11.x, g12.x, g22.x, nx, ny, my,
                                           hi1, hi2, k);
}

void col_update(
 dArray<Tv>& p, dArray<Tv>& v1, dArray<Tv>& v2,        
 dArray<Tv>& dp, dArray<Tv>& dv1, dArray<Tv>& dv2, 
 const Ti nx, const Ti ny, const Ti my, const Tv dt, const Ti k)
{
        dim3 threads (tx, ty, 1);
        dim3 blocks ((ny - 1) / tx + 1, (nx - 1) / ty + 1, 1);
        col_update_int<<<blocks, threads>>>(p.x, v1.x, v2.x, dp.x, dv1.x, dv2.x,
                                           nx, ny, my, dt, k);

}
#undef RSTRCT
