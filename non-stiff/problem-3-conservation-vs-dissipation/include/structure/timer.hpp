#pragma once
 
#include <ctime>

std::clock_t c_start = std::clock();
// your_algorithm
std::clock_t c_end = std::clock();


class HostTimer {
        private:
                std::clock_t cstart, cstop;
        public:
                void start(void) {
                       cstart = std::clock();
                }

                void stop(void) {
                       cstop = std::clock();
                }

                double elapsed(void) {
                        double elapsedms = (double) 1000.0 * (cstop - cstart) / CLOCKS_PER_SEC;
                        return elapsedms;
                }

};

class Timer {
        private:
                cudaEvent_t cstart, cstop;

        public:
                void start(void) {
                        cudaEventCreate(&cstart);
                        cudaEventCreate(&cstop);
                }

                void stop(void) {
                        cudaEventRecord(cstop);
                        cudaEventSynchronize(cstop);
                        cudaDeviceSynchronize();
                }

                float elapsed(void) {
                        float elapsed = 0.0;
                        cudaEventElapsedTime(&elapsed, cstart, cstop);
                        return elapsed;
                }


};
