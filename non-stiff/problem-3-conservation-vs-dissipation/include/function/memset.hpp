#pragma once
#include "function/host_err_check.hpp"

hostError_t hostMemset(void *out, int ch, size_t count)
{
        memset(out, ch, count);
        if (out == NULL) return hostMemsetFailure;
        return hostSuccess;
}

