#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
//#include <mkl.h>
#include <mpi.h>
#define ADDR(t, v) \
    &(t) { v }

const int DESC_LEN = 9;
typedef struct
{
    double *data;
    int *desc;
    int global_row;
    int global_col;
    int local_row;
    int local_col;
    int leading_dimension;
} Matrix;
int min(int a, int b)
{
    return (a < b ? a : b);
}

Matrix *create_matrix(int nproc_row, int nproc_col, int global_row, int global_col, int block_row, int block_col)
{
    int icontext;
    blacs_get_(ADDR(int, 0), ADDR(int, 0), &icontext);
    Matrix *matrix = malloc(sizeof(Matrix));
    matrix->desc = malloc(DESC_LEN * sizeof(int));
    matrix->global_col = global_col;
    matrix->global_row = global_row;
    matrix->local_col = global_col / nproc_col + block_col;
    matrix->local_row = global_row / nproc_row + block_row;
    matrix->leading_dimension = matrix->local_row + 8;
    matrix->data = malloc(matrix->leading_dimension * matrix->local_col * sizeof(double));
    int ierror;
    descinit_(matrix->desc, &(matrix->global_row), &(matrix->global_col), ADDR(int, block_row), ADDR(int, block_col), ADDR(int, 0), ADDR(int, 0), &icontext, &(matrix->leading_dimension), &ierror);
    printf("ierror %d\n", ierror);
    return matrix;
}

void mul(Matrix *a, Matrix *b, Matrix *c)
{
    assert(a->global_col == b->global_row);
    assert(a->global_row == c->global_row);
    assert(b->global_col == c->global_col);
    pdgemm_(ADDR(char, 'N'), ADDR(char, 'N'), &(a->global_row), &(b->global_col), &(a->global_col), ADDR(double, 1.0), &(a->data), ADDR(int, 1), ADDR(int, 1), &(a->desc), &(b->data), ADDR(int, 1), ADDR(int, 1), &(b->desc),
            ADDR(double, 0.0), &(c->data), ADDR(int, 1), ADDR(int, 1), &(c->desc));
}

void set(Matrix *a, int row, int col, double val)
{
    row++;
    col++;
    pdelset_(a->data, &row, &col, a->desc, &val);
}

double get(Matrix *a, int row, int col)
{
    double val;
    row++;
    col++;
    pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, a->data, &row, &col, a->desc);
    return val;
}

void free_matrix(Matrix *matrix)
{
    free(matrix->data);
    free(matrix->desc);
    free(matrix);
}

// #include <mpi.h>
// void blacs_pinfo_(int *, int *);
// void MPI_Dims_create(int, int, int[]);
// void blacs_get_(int *, int *, int *);
// void blacs_gridinit_(int *, char *, int *, int *);
// void blacs_gridinfo_(int *, int *, int *, int *, int *);
// void descinit_(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);

int main(int argc, char **argv)
{
    unsigned long const random_seed = 100;
    srand(random_seed);

    MPI_Init(&argc, &argv);

    if (argc != 2)
    {
        printf("(*'^')q < too few arguments!\n");
        return 0;
    }
    int m, n, k;
    m = n = k = atoi(argv[1]);
    int nproc, nproc_row, nproc_col, dims[2], ierror;
    int my_rank, my_row, my_col;
    int icontext;
    // int m = 1200, n = 800, k = 960;
    blacs_pinfo_(&my_rank, &nproc);
    dims[0] = dims[1] = 0;
    MPI_Dims_create(nproc, 2, dims);
    nproc_row = dims[0];
    nproc_col = dims[1];
    // sl_init_(&icontext, &nproc_row, &nproc_col);
    blacs_get_(ADDR(int, 0), ADDR(int, 0), &icontext);
    blacs_gridinit_(&icontext, ADDR(char, 'R'), &nproc_row, &nproc_col);
    blacs_gridinfo_(&icontext, &nproc_row, &nproc_col, &my_row, &my_col);
    if (my_rank == 0)
    {
        printf("pid %d\n", getpid());
        printf("dims: %d, %d\n", dims[0], dims[1]);
    }
    Matrix *mm = create_matrix(nproc_row, nproc_col, m, m, 40, 40);
    //blacs_barrier_(&icontext, ADDR(char, 'A'));

    printf("[%d] %d %d\n", my_rank, mm->local_row, mm->local_col);
    printf("[%d] %p %p\n", my_rank, mm->data, mm->desc);

    //pdelset_(mm->data, ADDR(int, 0), ADDR(int, 0), mm->desc, ADDR(double, 3.14));

    for (int i = 0; i < mm->global_row; i++)
    {
        for (int j = 0; j < mm->global_col; j++)
        {
            //    if (my_rank == 0)
            //     printf("%d %d\n", i, j);
            set(mm, i, j, (double)(i * 100 + j));
            //blacs_barrier_(&icontext, ADDR(char, 'A'));
        }
    }
    printf("[%d] %d %d\n", my_rank, mm->local_row, mm->local_col);
    blacs_barrier_(&icontext, ADDR(char, 'A'));

    for (int i = 0; i < mm->global_row; i++)
    {
        for (int j = 0; j < mm->global_col; j++)
        {
            double val = get(mm, i, j);
            if (my_rank == 0)
                printf("%03d ", (int)val);
        }
        if (my_rank == 0)
            printf("\n");
    }
    blacs_barrier_(&icontext, ADDR(char, 'A'));

    // printf("end");
    //free_matrix(mm);
    //free_matrix(mm);
    blacs_gridexit_(&icontext);
    blacs_exit_(ADDR(int, 0));
    //MPI_Finalize();

    return 0;
}