#pragma once
#include <mpi.h>

MPI_Datatype mpi_type(const char * val)
{
        return MPI_CHAR;
}

MPI_Datatype mpi_type(const float * val)
{
        return MPI_FLOAT;
}

MPI_Datatype mpi_type(const double * val)
{
        return MPI_DOUBLE;
}

MPI_Datatype mpi_type(const int * val)
{
        return MPI_INT;
}
