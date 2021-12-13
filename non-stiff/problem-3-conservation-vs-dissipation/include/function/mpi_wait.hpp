#pragma once

int mpi_wait(MPI_Request& req, MPI_Status& stat) {
        return MPI_Wait(&req, &stat);
}
