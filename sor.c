#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "sor.h"
#include "mpi_comm.h"
#include "utils.h"

#define true 1
#define false 0

//extern MPI_Comm CARTESIAN_COMM;

//Initialize an sor struct
sor* init_sor(int rank, int num_proc, int m_width, int m_height, float p_h,
                float p_w, float p_threshold)
{
    sor* block = malloc(sizeof(sor));
    MPI_Cart_coords(CARTESIAN_COMM, rank, 2, block->coords);
    block->generation = 1;
    block->h = pow(p_h, 2);
    block->w = p_w;
    block->threshold = p_threshold;
    block->proc_num = num_proc;
    block->rank_id = rank;
    block->matrix_width = m_width;
    block->matrix_height = m_height;
    block->grid_size = sqrt(num_proc);
    block->block_width = m_width/block->grid_size + 2;
    block->block_height = m_height/block->grid_size + 2;
    block->data = malloc(block->block_width * block->block_height * sizeof(float));
    block->next_data = malloc(block->block_width * block->block_height * sizeof(float));

    createDatatypes(block->block_width, block->block_height);

    //initial values assigment
    int i;
    for ( i = 0; i < block->block_width * block->block_height; i++ )
        block->data[i] = 1;

    // set top row to zero
    if ( block->coords[0] != 0 )
    {
        for ( i = 0; i < block->block_width; i++ )
            block->data[i] = 0;
    }

    // set bottom row to zero
    if ( block->coords[0] != block->grid_size - 1)
    {
        for ( i = (block->block_width - 1) * block->block_height - 1; i < block->block_width * block->block_height; i++ )
            block->data[i] = 0;
    }

    // set first column to zero
    if ( block->coords[1] != 0 )
    {
        for ( i = 0; i < (block->block_width - 1) * block->block_height; i+=block->block_width )
            block->data[i] = 0;
    }

    // set last column to zero
    if ( block->coords[1] != block->grid_size - 1 )
    {
        for ( i = block->block_width - 1; i < block->block_width * block->block_height; i+=block->block_width )
            block->data[i] = 0;
    }

    getNeighbors(block);

    //allocate swap-buffers for send/receive
    block->top_row = malloc(block->block_width * sizeof(float));
    block->bottom_row = malloc(block->block_width * sizeof(float));
    block->first_col = malloc(block->block_height * sizeof(float));
    block->last_col = malloc(block->block_height * sizeof(float));

    return block;
}

//distribute values to all the MPI processes from the master (rank 0)
// !Not debugged!
void dispach_data(sor* block)
{
    if (block->rank_id == 0)
    {
        float* tmp_matrix = create1Darray(block->matrix_width * block->matrix_height);

        MPI_Datatype square;
        MPI_Type_vector(block->block_height, block->block_height,
            block->block_height + block->matrix_width - block->block_width, MPI_FLOAT, &square );
        MPI_Type_commit(&square);

        int row_offset = block->block_width * block->block_height*
                        (block->matrix_width/block->block_width);
        int offset_mult = 0;

        int i, j;
        for ( i = 1, j = 1; i < block->proc_num; i++ , j++)
        {
            MPI_Request request;
            MPI_Isend((float*)&(tmp_matrix[j*block->block_width + offset_mult*row_offset]), 1, square, i, i, CARTESIAN_COMM, &request);

            if ( (j+1)*block->block_width == block->matrix_width )
            {
                offset_mult++;
                j = -1;
            }
        }

        for ( i = 0; i < block->block_height; i++ )
        {
            memcpy(block->data + i*block->block_width, tmp_matrix + i*block->matrix_width,
                    block->block_width*sizeof(float));
        }

        free(tmp_matrix);
    }
    else
    {
        MPI_Request request;
        MPI_Irecv((float*)block->data, 1, mpi_block, 0, block->rank_id, CARTESIAN_COMM, &request);

        MPI_Status status;
        int retval = MPI_Wait(&request, &status);
        if ( retval != MPI_SUCCESS )
        {
            printf("ERROR on receive() function\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
}


// send and receive top/bottom rows and first/last columns
// to processes around current process.
void sync_ext(sor* block)
{
    mpi_sync(block);
}

int compute(sor* block, int cnv_check)
{
    int i, width, height;
    float *data, *next_data;
    width = block->block_width;
    height = block->block_height;
    data = block->data;
    next_data = block->next_data;
    int w = block->w;
    int h = block->h;
    if ( cnv_check )
    {
        float max = 0;
        float tmp;
        ///compute red cells
        for ( i = width + 1; i <= width*(height-1)-2; i+=2)
        {
            next_data[i] = (1 - w)*data[i] + (w/4)*( data[i-width] + data[i-1] + data[i+1] + data[i+width]);
            tmp = abs(next_data[i] - data[i]);
            if ( tmp > max )
                max = tmp;
        }

        ///compute black cells
        for ( i = width + 2; i <= width*(height-1)-2; i+=2)
        {
            next_data[i] = (1 - w)*data[i] + (w/4)*( next_data[i-width] + next_data[i-1] + next_data[i+1] + next_data[i+width]);
            tmp = abs(next_data[i] - data[i]);
            if ( tmp > max )
                max = tmp;
        }

        //comunicate max diff to master.

        if (block->rank_id == 0)
        {
            float maxOfmaxs;
            MPI_Reduce(&max, &maxOfmaxs, 1, MPI_FLOAT, MPI_MAX, 0, CARTESIAN_COMM);
            if ( maxOfmaxs > block->threshold)
                return true;
        }
        else
        {
            MPI_Reduce(&max, NULL, 1, MPI_FLOAT, MPI_MAX, 0, CARTESIAN_COMM);
        }

    }
    else
    {
        ///compute red cells
        for ( i = width + 1; i <= width*(height-1)-2; i+=2)
        {
            next_data[i] = (1 - w)*data[i] + (w/4)*( data[i-width] + data[i-1] + data[i+1] + data[i+width]);
        }

        ///compute black cells
        for ( i = width + 2; i <= width*(height-1)-2; i+=2)
        {
            next_data[i] = (1 - w)*data[i] + (w/4)*( next_data[i-width] + next_data[i-1] + next_data[i+1] + next_data[i+width]);
        }
    }
    return false;
}

float compute_red(sor* block, int cnv_check)
{
    int i, width, height;
    float *data, *next_data;
    width = block->block_width;
    height = block->block_height;
    data = block->data;
    next_data = block->next_data;
    int w = block->w;
    int h = block->h;
    if ( cnv_check )
    {
        float max = 0;
        float tmp;
        ///compute red cells
        for ( i = width + 1; i <= width*(height-1)-2; i+=2)
        {
            next_data[i] = (1 - w)*data[i] + (w/4)*( h + data[i-width] + data[i-1] + data[i+1] + data[i+width]);
            tmp = abs(next_data[i] - data[i]);
            if ( tmp > max )
                max = tmp;
        }
        return max;

    }
    else
    {
        ///compute red cells
        for ( i = width + 1; i <= width*(height-1)-2; i+=2)
        {
            next_data[i] = (1 - w)*data[i] + (w/4)*( h + data[i-width] + data[i-1] + data[i+1] + data[i+width]);
        }
        return 0;
    }

}

float compute_black(sor* block, int cnv_check)
{
    int i, width, height;
    float *data, *next_data;
    width = block->block_width;
    height = block->block_height;
    data = block->data;
    next_data = block->next_data;
    int w = block->w;
    int h = block->h;
    if ( cnv_check )
    {
        float max = 0;
        float tmp;
        ///compute black cells
        for ( i = width + 2; i <= width*(height-1)-2; i+=2)
        {
            next_data[i] = (1 - w)*data[i] + (w/4)*( h + next_data[i-width] + next_data[i-1] + next_data[i+1] + next_data[i+width]);
            tmp = abs(next_data[i] - data[i]);
            if ( tmp > max )
                max = tmp;
        }
        float* swap_tmp = block->data;
        block->data = block->next_data;
        block->next_data = swap_tmp;
        //printf("swaped\n");
        return max;

    }
    else
    {
        ///compute black cells
        for ( i = width + 2; i <= width*(height-1)-2; i+=2)
        {
            next_data[i] = (1 - w)*data[i] + (w/4)*( h + next_data[i-width] + next_data[i-1] + next_data[i+1] + next_data[i+width]);
        }
        float* swap_tmp = block->data;
        block->data = block->next_data;
        block->next_data = swap_tmp;
    }

    return false;
}

