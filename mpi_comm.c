#include <stdio.h>
#include <stdlib.h>
#include "mpi_comm.h"

typedef enum direction {LEFT, RIGHT, UP, DOWN} direction;

//Initialize the MPI framework and return the
//number of active MPI processes.
int init_mpi()
{
    int numprocs;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    return numprocs;
}

//Create wrap-around cartesian MPI topology
//Parameter dims is a 2D array with the size of the grid of processes.
int init_mpi_cart(int* dims)
{
    int rank_id;

    //wrap-around
    int cyclic[2] = {1, 1};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, cyclic, 1, &CARTESIAN_COMM);

    //convert rank_id from COMM_WORLD to CARTESIAN
    int coords[2];
    memset(coords, '\0', sizeof(coords));
    if (CARTESIAN_COMM != MPI_COMM_NULL)
    {
        MPI_Cart_coords(CARTESIAN_COMM, rank_id, 2, coords);
        MPI_Cart_rank(CARTESIAN_COMM, coords, &rank_id);
    }
    else
    {
        printf("Could not properly set CARTESIAN_COMM_WORLD\n");
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    return rank_id;
}

//Create MPI datatypes needed for communication of processes
void createDatatypes(int width, int height)
{
    MPI_Type_vector(height, 1, width + 1, MPI_FLOAT, &mpi_column);
    MPI_Type_commit(&mpi_column);

    MPI_Type_vector(1, width, 0, MPI_FLOAT, &mpi_row);
    MPI_Type_commit(&mpi_row);

    MPI_Type_vector(1, width * height, 0, MPI_FLOAT, &mpi_block);
    MPI_Type_commit(&mpi_block);
}

// get neighbors process blocks and store them in the sor struct
void getNeighbors(sor* block)
{
    int current_row = block->coords[0];
    int current_col = block->coords[1];
    // Compute ranks
    int coords[2];
    // Upper
    int rank_upper;
    coords[0] = current_row - 1;
    coords[1] = current_col;
    MPI_Cart_rank(CARTESIAN_COMM, coords, &rank_upper);
    // Right
    int rank_right;
    coords[0] = current_row;
    coords[1] = current_col + 1;
    MPI_Cart_rank(CARTESIAN_COMM, coords, &rank_right);
    // Lower
    int rank_lower;
    coords[0] = current_row + 1;
    coords[1] = current_col;
    MPI_Cart_rank(CARTESIAN_COMM, coords, &rank_lower);
    // Left
    int rank_left;
    coords[0] = current_row;
    coords[1] = current_col - 1;
    MPI_Cart_rank(CARTESIAN_COMM, coords, &rank_left);

    block->rank_upper = rank_upper;
    block->rank_lower = rank_lower;
    block->rank_left = rank_left;
    block->rank_right = rank_right;
}

void mpi_sync(sor* block)
{
    //copy data to swap buffers so it wont conflict with received data
    memcpy(block->top_row, block->data, block->block_width * sizeof(float));
    memcpy(block->bottom_row, &(block->data[block->block_height * (block->block_width - 1) - 1]), block->block_width * sizeof(float));
    int i, j;
    for ( i = 0, j = 0; i < (block->block_width - 1) * block->block_height; i+=block->block_width , ++j )
        block->first_col[j] = block->data[i];
    for ( i = block->block_width - 1, j = 0; i < block->block_width * block->block_height; i+=block->block_width , ++j)
        block->last_col[j] = block->data[i];

    //sync with surrounding processes
    MPI_Status status;
    //left col
    MPI_Sendrecv(block->last_col, block->block_height, MPI_FLOAT, block->rank_right, LEFT,
                &(block->data[block->block_width - 1]), 1, mpi_column, block->rank_right, LEFT,
                CARTESIAN_COMM, &status );
    //right col
    MPI_Sendrecv(block->first_col, block->block_height, MPI_FLOAT, block->rank_left, RIGHT,
                block->data, 1, mpi_column, block->rank_left, RIGHT,
                CARTESIAN_COMM, &status );
    //first row
    MPI_Sendrecv(block->top_row, block->block_width, MPI_FLOAT, block->rank_upper, UP,
                block->data, 1, mpi_row, block->rank_upper, UP,
                CARTESIAN_COMM, &status );
    //last row
    MPI_Sendrecv(block->bottom_row, block->block_width, MPI_FLOAT, block->rank_lower, DOWN,
                &(block->data[block->block_height * (block->block_width - 1) - 1]), 1, mpi_row, block->rank_lower,
                DOWN, CARTESIAN_COMM, &status );

    block->generation += 1;
}
