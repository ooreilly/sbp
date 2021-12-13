#pragma once
#ifndef NDEBUG
#define mpi_err_check(err) mpi_err_check_((err), __FILE__, __LINE__, __func__)
#else
#define mpi_err_check(err) 
#endif

#ifndef NDEBUG
__inline__ void mpi_err_check_(int err, const char *file, int line, const char *func) {
        if (err != MPI_SUCCESS) {
                char error_string[2048];
                int length_of_error_string;
                MPI_Error_string((err), error_string, &length_of_error_string);
                fprintf(stderr, "MPI error: %s:%i %s(): %s\n", file, 
                        line, func, error_string);
                MPI_Abort(MPI_COMM_WORLD, err);
                fflush(stderr);
                exit(EXIT_FAILURE);
        }
}
#else
__inline__ void mpi_err_check_(int err, const char *file, int line,
                               const char *func) {}
#endif

