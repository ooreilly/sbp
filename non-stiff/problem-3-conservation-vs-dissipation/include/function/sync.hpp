#pragma once

#define sync(buf) { sync_( (buf), __FILE__, __LINE__); }

template <typename Tv>
__inline__ void sync_(DeviceArray<Tv>& buf, const char *file, int line) 
{
        cudaErrCheck_(
            cudaMemcpy(buf.u, buf.v, buf.num_bytes, cudaMemcpyHostToDevice),
            file, line);
}
