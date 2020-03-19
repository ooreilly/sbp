cmake_minimum_required(VERSION 3.10)
project(ACOUSTIC VERSION 1.0 LANGUAGES CXX CUDA)
      
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(MPI REQUIRED )
set(MPI_DETERMINE_LIBRARY_VERSION ON)

include(CTest)
enable_testing()

include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
if (DEFINED ENV{ARCH})
        set(ARCH $ENV{ARCH})
else()
        set(ARCH sm_75)
endif()
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xptxas=-v -lineinfo -use_fast_math -arch=${ARCH}")

include_directories(src)
include_directories(include)

add_executable(acoustic src/acoustic.cu)    

