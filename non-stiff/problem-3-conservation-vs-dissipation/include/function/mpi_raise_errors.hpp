#pragma once
#include <mpi.h>

void mpi_raise_errors() {
        MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN); 
}


