#pragma once
#include "structure/bounds.hpp"

template <typename Array>
size_t vtk(const char *fname, 
           const Array& x,
           const Array& y,
           const Array& z,
           const Bounds& bounds)
{
        FILE *fh;
        size_t count = 0;
        fh = fopen(fname, "w");
        
        if (!fh) {
                return count;
        }

        size_t nx = bounds.x[1] - bounds.x[0];
        size_t ny = bounds.y[1] - bounds.y[0];
        size_t nz = bounds.z[1] - bounds.z[0];

        size_t numpts = nx * ny * nz;

        // Header
        fprintf(fh, "# vtk DataFile Version 4.2\n");
        fprintf(fh, "vtk output\n");
        fprintf(fh, "ASCII\n");
        fprintf(fh, "DATASET STRUCTURED_GRID\n");
        fprintf(fh, "DIMENSIONS %ld %ld %ld\n", 
                        nx, ny, nz);
        fprintf(fh, "POINTS %ld float\n", numpts);
        
        // Coordinates
        for (size_t k = bounds.z[0]; k < bounds.z[1]; ++k) {
        for (size_t j = bounds.y[0]; j < bounds.y[1]; ++j) {
        for (size_t i = bounds.x[0]; i < bounds.x[1]; ++i) {
                fprintf(fh, "%f %f %f\n", 
                        x(i, j, k), y(i, j, k), z(i, j, k));
                count++;
        }
        }
        }
        assert(count == numpts);

        fclose(fh);

        return count;
}

template <typename Array>
size_t vtk(const char *fname,
           const Array& data,
           const Bounds bounds,
           const char *label="u"
           )
{
        size_t count = 0;
        size_t nx = bounds.x[1] - bounds.x[0];
        size_t ny = bounds.y[1] - bounds.y[0];
        size_t nz = bounds.z[1] - bounds.z[0];

        size_t numpts = nx * ny * nz;
        
        
        FILE *fh = fopen(fname, "a");
        
        if (!fh) {
                return count;
        }

        fprintf(fh, "POINT_DATA %ld \n", numpts);
        fprintf(fh, "FIELD scalar 1\n");
        fprintf(fh, "%s 1 %ld float\n", label, numpts);

        for (size_t k = bounds.z[0]; k < bounds.z[1]; ++k) {
        for (size_t j = bounds.y[0]; j < bounds.y[1]; ++j) {
        for (size_t i = bounds.x[0]; i < bounds.x[1]; ++i) {
                fprintf(fh, "%f\n", 
                        data(i, j, k));
                count++;
        }
        }
        }
        assert(count == numpts);

        fclose(fh);

        return count;
}
