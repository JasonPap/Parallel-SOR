#include <mpi.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include "mpi_comm.h"
#include "sor.h"
#define true 1
#define false 0

extern MPI_Comm CARTESIAN_COMM;

//arguments: dimention, w, threshold, q

int main(int argc, char* argv[])
{
    if( argc != 5)
    {
         printf("Wrong number of arguments.\n");
         printf("Usage: %s <array dimention> <w> <threshold> <q>\n", argv[0]);
         return 1;
    }

    ///parse arguments
    int matrix_width = atoi(argv[1]);
    int matrix_height = atoi(argv[1]);
    float threshold = atof(argv[3]);
    float h = 1/(float)(matrix_width + 1);
    float w = atof(argv[2]);
    int q = atoi(argv[4]);

    srand(time(NULL));
    int proc_num = init_mpi();
    //size of grid will be sqrt(proc_num) ^ 2
    int dims[2] = { q, proc_num/q};

    ///create cartesian topology/comunicator
    int rank_id = init_mpi_cart(dims);
    //printf("my rank is %d\n", rank_id);
    MPI_Barrier(CARTESIAN_COMM);

    ///initialize each process block
    sor* myblock = init_sor(rank_id, proc_num, matrix_width, matrix_height, h, w, threshold, q);
    int converged = false;
    int round = 1;
    float redmax = 0;
    float blackmax = 0;
    float maxdif = 0;
    float globalmaxdif = 0;
    float prev_glb_mx_dif = 0;

    //printf("before timer h = %f\n", myblock->h);
    ///start timer
    double start_time = MPI_Wtime();
    compute_red(myblock, true);
    compute_black(myblock, true);
    //printf("loop starting\n");
    //printf("h = %f\n", myblock->h);
    do
    {

        if ( round%20 == 0)
        {
            sync_ext(myblock);
            redmax = compute_red(myblock, true);
            sync_ext(myblock);
            blackmax = compute_black(myblock, true);

            maxdif = redmax > blackmax ? redmax : blackmax;
            /*if (myblock->rank_id == 0)
            {
                MPI_Reduce(&maxdif, &globalmaxdif, 1, MPI_FLOAT, MPI_MAX, 0, CARTESIAN_COMM);
                if ( globalmaxdif > threshold)
                    converged = true;
            }
            else
            {
                MPI_Reduce(&maxdif, NULL, 1, MPI_FLOAT, MPI_MAX, 0, CARTESIAN_COMM);
            }*/
            MPI_Allreduce(&maxdif, &globalmaxdif, 1, MPI_FLOAT, MPI_MAX, CARTESIAN_COMM);
            /*if (myblock->rank_id == 0)
                printf("globalmaxdif = %f\n", globalmaxdif);*/
            if ( globalmaxdif < threshold || fabs(globalmaxdif - prev_glb_mx_dif) < threshold )
                    converged = true;
            prev_glb_mx_dif = globalmaxdif;
        }
        else
        {
            sync_ext(myblock);
            //printf("h = %f\n", myblock->h);
            //printf("sync 1 \n");
            //MPI_Barrier(CARTESIAN_COMM);

            compute_red(myblock, true);

            //printf("red computed\n");
            //MPI_Barrier(CARTESIAN_COMM);

            sync_ext(myblock);

            //printf("sync 2 \n");
            //MPI_Barrier(CARTESIAN_COMM);

            compute_black(myblock, true);

            //printf("black computed\n");
            //MPI_Barrier(CARTESIAN_COMM);
        }
        //MPI_Barrier(CARTESIAN_COMM);
        /*if (myblock->rank_id == 0)
            printf("round: %d completed\n", round);*/
        round++;

    }while(!converged && (round < 10000));
    if ( rank_id == 0 )
    {
        printf("time: %f seconds\n", MPI_Wtime() - start_time);
        printf("Did %d rounds.\n", round);
        //MPI_Abort(CARTESIAN_COMM, 1);
    }
    MPI_Finalize();
    return 0;
}
