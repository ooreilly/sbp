#ifndef HELPER_HPP
#define HELPER_HPP

typedef float Tv;
typedef int Ti;
// Runge-Kutta 4 coefficients
const static int rk4_n = 5;
static Tv a[5] = {0.0,
                             -567301805773.0/1357537059087,
                             -2404267990393.0/2016746695238,
                             -3550918686646.0/2091501179385,
                             -1275806237668.0/842570457699};

static Tv b[5] = {1432997174477.0/9575080441755,
                             5161836677717.0/13612068292357,
                             1720146321549.0/2090206949498,
                             3134564353537.0/4481467310338,
                             2277821191437.0/14882151754819};
static Tv c[5] = {0.0,
                             1432997174477.0/9575080441755,
                             2526269341429.0/6820363962896,
                             2006345519317.0/3224310063776,
                             2802321613138.0/2924317926251};





inline Tv gaussian(Tv t)
{
        return exp( -5e2 * pow(t  - 2.0, 2));
}

inline Tv ricker(Tv t)
{
        Tv pi2 = M_PI*M_PI;
        Tv t0 = 1.7;
        return (1 - 2*pi2*pow(t - t0, 2))* exp( - pi2 * pow(t  - t0, 2));
}

void write_csv(std::string& fname, hArray<Tv>& u, hArray<Tv>& t) {
        FILE *fh;
        fh = fopen(fname.c_str(), "w");

        if (!fh) fprintf(stderr, "Failed to open file\n");

        size_t num_elem = u.size / t.size;

        for (size_t i=0; i < t.size; ++i) {
                        fprintf(fh, "%.17g ", t[i]);
                        for (size_t j=0; j < num_elem; ++j) {
                                fprintf(fh, "%.17g ", u[j + i * num_elem]);
                        }
                        fprintf(fh, "\n");
        }
        fclose(fh);
}

size_t write_vtk(const std::string& fname,
                 const Tv* data,
                 hArray<Tv>& x,
                 hArray<Tv>& y,
                         const int nx, const int ny
                         )
{
        size_t count = 0;
        int numpts = nx * ny;
        
        FILE *fh = fopen(fname.c_str(), "w");
        
        if (!fh) {
                return count;
        }

        // Header
        fprintf(fh, "# vtk DataFile Version 4.2\n");
        fprintf(fh, "vtk output\n");
        fprintf(fh, "ASCII\n");
        fprintf(fh, "DATASET STRUCTURED_GRID\n");
        fprintf(fh, "DIMENSIONS %d %d %d\n", nx, ny, 1);
        fprintf(fh, "POINTS %d float\n", numpts);
        
        // Coordinates
        for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
                int idx = j + i*ny;
                fprintf(fh, "%f %f %f\n", (float)x[idx], (float)y[idx], 0.0);
                count++;
        }
        }

        fprintf(fh, "POINT_DATA %d \n", numpts);
        fprintf(fh, "FIELD scalar 1\n");
        fprintf(fh, "data 1 %d float\n", numpts);
        for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
                int idx = j + i*ny;
                fprintf(fh, "%f \n", (float) data[idx]);
                count++;
        }
        }

        fclose(fh);

        return count;
}

#define dump2d(x, m, n) dump2d_(x, #x, m, n)

template<typename T>
void dump2d_(hArray<T>& x, const char *name, int m, int n, FILE *stream=NULL)
{
        if (!stream) stream = stdout;
        fprintf(stream, "%s = [\n", name);
                for (int j = n - 1; j >= 0; --j) {
                fprintf(stream, "       ");
        for (int i = 0; i < m; ++i) {
                        fprintf(stream, " %3.3f", x[j + i * n]);
        }
                fprintf(stream, "\n");
        }
        fprintf(stream, "     ]\n");
}

#endif
