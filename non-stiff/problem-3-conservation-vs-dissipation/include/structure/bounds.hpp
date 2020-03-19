#pragma once
#include "function/copy_initializer_list.hpp"

class Bounds {
        public:
        size_t x[2];
        size_t y[2];
        size_t z[2];

        Bounds(
                std::initializer_list<size_t> x_={0,1}, 
                std::initializer_list<size_t> y_={0,1}, 
                std::initializer_list<size_t> z_={0,1})
               {
                copy_initializer_list(x, x_);
                assert(x[1] > x[0]);

                copy_initializer_list(y, y_);
                assert(y[1] > y[0]);

                copy_initializer_list(z, z_);
                assert(z[1] > z[0]);
        }

        void dump(const char *msg="", FILE *stream=NULL) const {
                        if (stream == NULL) stream = stdout;
                        fprintf(stream, "%s", msg);
                        fprintf(
                            stream,
                            "x = (%ld, %ld), y = (%ld, %ld), z = (%ld, %ld) \n",
                            x[0], x[1], y[0], y[1], z[0], z[1]);
        }

};
