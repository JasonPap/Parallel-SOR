#ifndef SOR_H
#define SOR_H

#ifdef	__cplusplus
extern "C" {
#endif

    typedef struct
    {
        float* data;
        float* next_data;
        float* top_row;
        float* bottom_row;
        float* first_col;
        float* last_col;

        float h;
        float w;
        float threshold;

        int coords[2];
        int proc_num;
        int rank_id;
        int grid_size;
        int generation;

        int rank_upper;
        int rank_lower;
        int rank_left;
        int rank_right;

        int matrix_width;
        int matrix_height;

        int block_width;
        int block_height;
    }sor;

    sor* init_sor(int rank, int num_proc, int matrix_width,
                int matrix_height, int p_h, int p_w, float p_threshold);
    int compute(sor* block, int cnv_check);
    float compute_red(sor* block, int cnv_check);
    float compute_black(sor* block, int cnv_check);
    void dispach_data(sor* block);

#ifdef	__cplusplus
 }
#endif

#endif //SOR_H
