#ifndef MPI_COMM_H
#define MPI_COMM_H

#include <mpi.h>
#include "sor.h"

#ifdef	__cplusplus
extern "C" {
#endif
    //------------------------------
    // Communicator
    //------------------------------
    MPI_Comm CARTESIAN_COMM;

    //-------------------------------------------
    // Datatypes
    //-------------------------------------------
    MPI_Datatype mpi_column;
    MPI_Datatype mpi_row;
    MPI_Datatype mpi_block;

    //-------------------------------------------
    // Functions
    //-------------------------------------------
    int init_mpi();
    int init_mpi_cart(int* dims);
    void createDatatypes(int width, int height);
    void getNeighbors(sor* block);
    void mpi_sync(sor* block);

#ifdef	__cplusplus
 }
#endif

#endif //MPI_COMM_H
