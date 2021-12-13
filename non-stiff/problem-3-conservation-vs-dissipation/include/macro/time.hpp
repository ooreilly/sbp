#pragma once
#include "structure/timer.hpp"
                                                     
#define time(x) {                                    \
        Timer timer;                                 \
        timer.start();                               \
        timer.stop();                                \
        printf("%s took %g ms \n", #x, timer.elapsed()); \
}
