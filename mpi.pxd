# Pyrex module for MPI

cdef extern from "mpi.h":

    enum MPI_Error:
        MPI_SUCCESS

    ctypedef struct MPI_Comm_Imp:
        pass
    ctypedef MPI_Comm_Imp *MPI_Comm
    MPI_Comm MPI_COMM_WORLD

    int MPI_Init(int *, char ***)
    int MPI_Finalize()
    int MPI_Comm_rank(MPI_Comm, int *)

# end of file
