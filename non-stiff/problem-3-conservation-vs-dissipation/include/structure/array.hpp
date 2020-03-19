#pragma once
#include "function/malloc.hpp"
#include "function/free.hpp"
#include "function/pitch.hpp"
#include "function/memzero.hpp"
#include "structure/bounds.hpp"
#include "function/bounds.hpp"

template <typename Tv>
class DeviceArray;

void _dump_shape(size_t nx, size_t ny, size_t nz, const char *msg = "",
                FILE *stream = NULL) {
        printf("%s nx = %ld ny = %ld nz = %ld \n", msg, nx, ny, nz);
}


template <typename Tv>
class HostArray {

        public:
                size_t nx=1, ny=1, nz=1, line, slice, size, num_bytes = 0;
                Tv *u;

                HostArray() {
                        u = NULL;
                }

                HostArray(size_t nx_, size_t ny_=1, size_t nz_=1) {
                        init(nx_, ny_, nz_);
                }
                
                HostArray<Tv>& operator=(const HostArray<Tv>& rhs) {
                        if (rhs.num_bytes > 0) {
                                init(rhs.nx, rhs.ny, rhs.nz);
                                memcpy(u, rhs.u, num_bytes);
                        }
                        return *this;
                }

                template <typename Array>
                HostArray(const Array& in) {
                        init(in.nx, in.ny, in.nz);
                }

                HostArray(const Bounds& in) {
                        init(in.x[1] - in.x[0], in.y[1] - in.y[0],
                             in.z[1] - in.z[0]);
                }

                void init(size_t nx_=1, size_t ny_=1, size_t nz_=1) {
                        assert(nx_ > 0);
                        assert(ny_ > 0);
                        assert(nz_ > 0);
                        nx = nx_;
                        ny = ny_;
                        nz = nz_;
                        size_t mx = pitch<Tv>(nx);
                        size = mx * ny * nz;
                        line = mx;
                        slice = mx * ny;
                        num_bytes = size * sizeof(Tv);
                        hmalloc((void**)&u, num_bytes);
                        hzero(u, num_bytes);
                }


                HostArray<Tv> host() {
                        HostArray<Tv> out(nx, ny, nz);
                        memcpy(out.u, u, num_bytes);
                        return out;
                }

                DeviceArray<Tv> device() const {
                        DeviceArray<Tv> dev(nx, ny, nz);
                        cudaErrCheck(cudaMemcpy(dev(), u, num_bytes,
                                                cudaMemcpyHostToDevice));
                        return dev;
                }


                const Tv * operator()() const{
                        return u;
                }
                Tv * operator()(){
                        return u;
                }

                __inline__ Tv operator()(size_t xi, size_t yj=0,
                                         size_t zk=0) {
                        assert(xi < nx);
                        assert(yj < ny);
                        assert(zk < nz);
                        return u[xi + line * yj + slice * zk];
                }

                __inline__ Tv operator() (size_t xi, size_t yj=0,
                                         size_t zk=0) const {
                        assert(xi < nx);
                        assert(yj < ny);
                        assert(zk < nz);
                        return u[xi + line * yj + slice * zk];
                }

                void dump(const char *msg="", FILE *stream=NULL) const{
                        if (stream == NULL) stream = stdout;
                        fprintf(stream, "%s", msg);
                        for (int k = 0; k < nz; ++k) {
                        fprintf(stream, "(:,:, %d) = [\n", k);
                        for (int j = ny - 1; j >= 0; --j) {
                                fprintf(stream, "       ");
                                for (int i = 0; i < nx; ++i) {
                                        Tv val = operator()(i, j, k);
                                        if (typeid(Tv) == typeid(float) ||
                                            typeid(Tv) == typeid(double))
                                                fprintf(stream, " %3.3f", (float)val);
                                        if (typeid(Tv) == typeid(int))
                                                fprintf(stream, " %d", (int)val);
                                        if (typeid(Tv) == typeid(size_t))
                                                fprintf(stream, " %ld", (size_t)val);
                                }
                                fprintf(stream, "\n");
                        }
                        fprintf(stream, "     ]\n");
                        }
                }

                void dump_shape(const char *msg="", FILE *stream=NULL) const {
                        _dump_shape(nx, ny, nz, msg, stream);
                }

                Bounds bounds() const {
                        return fbounds(nx, ny, nz);
                }

                ~HostArray(void) {
                        if (u != NULL) hfree((void*)u);
                }
};

template <typename Tv>
class DeviceArray {

        public:
                size_t nx=1, ny=1, nz=1, line, slice, size, num_bytes=0;
                Tv *u, *v;

                DeviceArray() {
                        u = NULL;
                }

                DeviceArray(size_t nx_, size_t ny_=1, size_t nz_=1) {
                        init(nx_, ny_, nz_);
                }
                
                DeviceArray<Tv>& operator=(const DeviceArray<Tv>& rhs) {
                        if (rhs.num_bytes > 0) {
                                init(rhs.nx, rhs.ny, rhs.nz);
                                cudaErrCheck(
                                    cudaMemcpy(u, rhs(), num_bytes,
                                               cudaMemcpyDeviceToDevice));
                        }
                        return *this;
                }

                template <typename Array>
                DeviceArray(Array& in) {
                        init(in.nx, in.ny, in.nz);
                }

                DeviceArray(const Bounds& in) {
                        init(in.x[1] - in.x[0], in.y[1] - in.y[0],
                             in.z[1] - in.z[0]);
                }

                HostArray<Tv> host() const {
                        HostArray<Tv> host(nx, ny, nz);
                        cudaErrCheck(cudaMemcpy(host(), u, num_bytes,
                                                cudaMemcpyDeviceToHost));
                        return host;
                }

                void init(size_t nx_=1, size_t ny_=1, size_t nz_=1) {
                        assert(nx_ > 0);
                        assert(ny_ > 0);
                        assert(nz_ > 0);
                        nx = nx_;
                        ny = ny_;
                        nz = nz_;
                        size_t mx = pitch<Tv>(nx);
                        size = mx * ny * nz;
                        line = mx;
                        slice = mx * ny;
                        num_bytes = size * sizeof(Tv);
                        dmalloc((void**)&u, num_bytes);
                        hmalloc((void**)&v, num_bytes);
                        dzero(u, num_bytes);
                        hzero(v, num_bytes);
                }


                Tv * operator()(){
                        return u;
                }

                const Tv * operator()() const{
                        return u;
                }

                ~DeviceArray(void) {
                        if (u != NULL)
                        dfree((void*)u);
                }

                void dump(const char *msg="", FILE *stream=NULL) const {
                        HostArray<Tv> host(nx, ny, nz);
                        cudaErrCheck(cudaMemcpy(host(), u, num_bytes,
                                                cudaMemcpyDeviceToHost));
                        host.dump(msg, stream);
                }

                void dump_shape(const char *msg="", FILE *stream=NULL) const {
                        _dump_shape(nx, ny, nz, msg, stream);
                }

                Bounds bounds() const {
                        return fbounds(nx, ny, nz);
                }
};

