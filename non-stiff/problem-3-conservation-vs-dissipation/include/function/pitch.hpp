#pragma once

// Determine pitch so that it is divisible by pitchByteSize
template <typename Tv>
size_t pitch(const size_t size, size_t pitchByteSize = 128) {
        size_t numBytes = size * sizeof(Tv);
        int v = pitchByteSize *
                (((int)numBytes - sizeof(Tv)) / pitchByteSize + 1) / sizeof(Tv);
        return v;
}

