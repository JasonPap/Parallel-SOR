#include <mpi.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include "mpi_comm.h"
#include "sor.h"
#define true 1
#define false 0

extern MPI_Comm CARTESIAN_COMM;

int main(void)
{
    ///parse arguments
    int matrix_width = 8;
    int matrix_height = 8;
    float threshold = 1;
    int h = 1;
    int w = 1;

    srand(time(NULL));
    int proc_num = init_mpi();
    //size of grid will be sqrt(proc_num) ^ 2
    int dims[2] = {(int) sqrt(proc_num), (int) sqrt(proc_num)};

    ///create cartesian topology/comunicator
    int rank_id = init_mpi_cart(dims);
    MPI_Barrier(CARTESIAN_COMM);

    ///initialize each process block
    sor* myblock = init_sor(rank_id, proc_num, matrix_width, matrix_height, h, w, threshold);
    int converged = 0;
    int round = 1;

    //start timer
    double start_time = MPI_Wtime();

    do
    {
        sync_ext(myblock);
        if ( round%20 == 0)
            converged = compute(myblock, true);
        else
            compute(myblock, false);


    }while(!converged);
    if ( rank_id == 0 )
    {
        printf("time: %f", MPI_Wtime() - start_time);
        MPI_Abort(CARTESIAN_COMM, 1);
    }
    return 0;
}
