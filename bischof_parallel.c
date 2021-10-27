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

const int block_row = 1;
const int block_col = 1;

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
    return matrix;
}

void set(Matrix *a, int row, int col, double val)
{
    row++;
    col++;
    pdelset_(a->data, &row, &col, a->desc, &val);
}

double get(Matrix *a, int row, int col)
{
    assert(row < a->global_row);
    assert(col < a->global_col);
    assert(0 <= row);
    assert(0 <= col);
    double val;
    row++;
    col++;
    pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, a->data, &row, &col, a->desc);
    return val;
}

void print_matrix(char *msg, Matrix *matrix, int rank)
{
    if (rank == 0)
    {
        // printf("================================\n");
        printf("%s", msg);
        printf("[\n");
    }
    for (int i = 0; i < matrix->global_row; i++)
    {
        if (rank == 0)
        {
            printf("[");
        }
        for (int j = 0; j < matrix->global_col; j++)
        {
            double val = get(matrix, i, j);
            if (rank == 0)
                printf("%f,", val);
        }
        if (rank == 0)
            printf("],\n");
    }
    if (rank == 0)
        printf("]\n");
}

void pdgemv_wrap(char trans, int m, int n, double alpha, Matrix *a, int row_a, int col_a, Matrix *x, int row_x, int col_x, int inc_x, double beta, Matrix *y, int row_y, int col_y, int inc_y)
{
    row_a++;
    row_x++;
    row_y++;
    col_a++;
    col_x++;
    col_y++;
    pdgemv_(&trans, &m, &n, &alpha, a->data, &row_a, &col_a, a->desc, x->data, &row_x, &col_x, x->desc, &inc_x, &beta, y->data, &row_y, &col_y, y->desc, &inc_y);
}

void pdgemm_wrap(char trans_a, char trans_b, int m, int n, int k, double alpha, Matrix *a, int row_a, int col_a, Matrix *b, int row_b, int col_b, double beta, Matrix *c, int row_c, int col_c)
{
    // pdgemm_(ADDR(char, 'N'), ADDR(char, 'N'), &(a->global_row), &(b->global_col), &(a->global_col), ADDR(double, 1.0), &(a->data), ADDR(int, 1), ADDR(int, 1), &(a->desc), &(b->data), ADDR(int, 1), ADDR(int, 1), &(b->desc),
    //            ADDR(double, 0.0), &(c->data), ADDR(int, 1), ADDR(int, 1), &(c->desc));
    row_a++;
    row_b++;
    row_c++;
    col_a++;
    col_b++;
    col_c++;
    pdgemm_(&trans_a, &trans_b, &m, &n, &k, &alpha, a->data, &row_a, &col_a, a->desc, b->data, &row_b, &col_b, b->desc, &beta, c->data, &row_c, &col_c, c->desc);
}

void pdsymm_wrap(char side, char uplo, int m, int n, double alpha, Matrix *a, int row_a, int col_a, Matrix *b, int row_b, int col_b, double beta, Matrix *c, int row_c, int col_c)
{
    row_a++;
    row_b++;
    row_c++;
    col_a++;
    col_b++;
    col_c++;
    pdsymm_(&side, &uplo, &m, &n, &alpha, a->data, &row_a, &col_a, a->desc, b->data, &row_b, &col_b, b->desc, &beta, c->data, &row_c, &col_c, c->desc);
}

void pdsyr2k_wrap(char uplo, char trans, int m, int k, double alpha, Matrix *a, int row_a, int col_a, Matrix *b, int row_b, int col_b, double beta, Matrix *c, int row_c, int col_c)
{
    row_a++;
    row_b++;
    row_c++;
    col_a++;
    col_b++;
    col_c++;
    pdsyr2k_(&uplo, &trans, &m, &k, &alpha, a->data, &row_a, &col_a, a->desc, b->data, &row_b, &col_b, b->desc, &beta, c->data, &row_c, &col_c, c->desc);
}

void pdgeqrf_wrap(int m, int n, Matrix *matrix, int row, int col, double *tau)
{
    assert(m >= n);
    row++;
    col++;
    int info;
    double *work = malloc(1 * sizeof(double));

    // calculate length of work
    int lwork = -1;
    pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info);
    lwork = (int)(work[0]);

    // reallocate and compute QR
    free(work);
    work = malloc(lwork * sizeof(double));
    pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info);

    //    printf("lwork %f\n", work[0]);
    //  printf("info %d\n", info);
    free(work);
}

void pdgeqrt_wrap(int rank, int proc_row, int proc_col, int m, int n, Matrix *matrix, int row, int col, Matrix *T)
{
    int icontext;
    blacs_get_(ADDR(int, 0), ADDR(int, 0), &icontext);

    assert(m >= n);
    assert(T->global_row >= n);
    assert(T->global_col >= n);
    double *tau = malloc(n * sizeof(double));
    pdgeqrf_wrap(m, n, matrix, row, col, tau);
    set(T, 0, 0, tau[0]);
    Matrix *Y = create_matrix(proc_row, proc_col, m, n, block_row, block_col);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            set(T, i, j, 0.0);
        }
    }

    for (int i = 0; i < m; i++)
    {
        if (rank == 0)
        {
            printf("i %d\n", i);
        }

        if (i == 0)
        {
            set(Y, i, 0, 1.0);
        }
        else
        {
            set(Y, i, 0, get(matrix, row + i, col));
        }
    }

    Matrix *y = create_matrix(proc_row, proc_col, m, 1, block_row, block_col);
    for (int j = 1; j < n; j++)
    {
        if (rank == 0)
        {
            printf("j %d\n", j);
            printf("pid %d\n", getpid());
        }

        for (int i = 0; i < m; i++)
        {
            if (i < j)
            {
                set(y, i, 0, 0.0);
            }
            else if (j == i)
            {
                set(y, i, 0, 1.0);
            }
            else
            {
                set(y, i, 0, get(matrix, row + i, col + j));
            }
        }
        Matrix *z = create_matrix(proc_row, proc_col, j, 1, block_row, block_col);
        Matrix *tmp = create_matrix(proc_row, proc_col, j, 1, block_row, block_col);

        pdgemm_wrap('T', 'N', j, 1, m, 1.0, Y, 0, 0, y, 0, 0, 0.0, tmp, 0, 0);

        pdgemm_wrap('N', 'N', j, 1, j, -tau[j], T, 0, 0, tmp, 0, 0, 0.0, z, 0, 0);

        for (int i = 0; i < m; i++)
        {
            if (rank == 0)
            {
                printf("loop %d/%d\n", i, m - 1);
                printf("%d\n", i);
            }
            // printf("[%d] start %d\n", rank, i);
            double ttt = get(y, i, 0);
            if (rank == 0)
            {
                printf("ok %f\n", ttt);
            }
            // printf("[%d] done %d\n", rank, i);
            set(Y, i, j, ttt);
            // printf("%d\n", icontext);
        }
        for (int i = 0; i < j; i++)
        {
            set(T, i, j, get(z, i, 0));
        }
        set(T, j, j, tau[j]);

        free_matrix(tmp);
        free_matrix(z);
    }
    free_matrix(y);
    free_matrix(Y);
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
    int rank, my_row, my_col;
    int icontext;
    // int m = 1200, n = 800, k = 960;
    blacs_pinfo_(&rank, &nproc);
    dims[0] = dims[1] = 0;
    MPI_Dims_create(nproc, 2, dims);
    nproc_row = dims[0];
    nproc_col = dims[1];
    // sl_init_(&icontext, &nproc_row, &nproc_col);
    blacs_get_(ADDR(int, 0), ADDR(int, 0), &icontext);
    blacs_gridinit_(&icontext, ADDR(char, 'R'), &nproc_row, &nproc_col);
    blacs_gridinfo_(&icontext, &nproc_row, &nproc_col, &my_row, &my_col);
    if (rank == 0)
    {
        printf("pid %d\n", getpid());
        printf("dims: %d, %d\n", dims[0], dims[1]);
    }
    Matrix *a = create_matrix(nproc_row, nproc_col, m, m, 1, 1);
    Matrix *b = create_matrix(nproc_row, nproc_col, m, m, 1, 1);
    Matrix *c = create_matrix(nproc_row, nproc_col, m, m, 1, 1);
    // blacs_barrier_(&icontext, ADDR(char, 'A'));

    // pdelset_(mm->data, ADDR(int, 0), ADDR(int, 0), mm->desc, ADDR(double, 3.14));

    for (int i = 0; i < a->global_row; i++)
    {
        for (int j = 0; j < a->global_col; j++)
        {
            //    if (my_rank == 0)
            //     printf("%d %d\n", i, j);
            set(a, i, j, (double)rand() / (double)(RAND_MAX));
            set(b, i, j, (double)rand() / (double)(RAND_MAX));
            // blacs_barrier_(&icontext, ADDR(char, 'A'));
        }
    }

    // print_matrix("a = ", a, rank);
    // print_matrix("b = ", b, rank);
    // printf("[%d] %d %d\n", my_rank, mm->local_row, mm->local_col);
    blacs_barrier_(&icontext, ADDR(char, 'A'));
    // pdgemm_wrap('N', 'N', m, m, m, 1.0, a, 0, 0, b, 0, 0, 0.0, c, 0, 0);

    Matrix *T = create_matrix(nproc_row, nproc_col, m, 12, block_row, block_col);
    pdgeqrt_wrap(rank, nproc_row, nproc_col, m, 12, a, 0, 0, T);

    //  print_matrix("c = ", c, rank);
    blacs_barrier_(&icontext, ADDR(char, 'A'));

    // printf("end");
    free_matrix(a);
    free_matrix(b);
    free_matrix(c);

    // free_matrix(mm);
    blacs_gridexit_(&icontext);
    blacs_exit_(ADDR(int, 0));
    // MPI_Finalize();

    return 0;
}