#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
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

// #include <mpi.h>
// void blacs_pinfo_(int *, int *);
// void MPI_Dims_create(int, int, int[]);
// void blacs_get_(int *, int *, int *);
// void blacs_gridinit_(int *, char *, int *, int *);
// void blacs_gridinfo_(int *, int *, int *, int *, int *);
// void descinit_(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);

int main(int argc, char **argv)
{

    MPI_Init(&argc, &argv);
    if (argc != 2)
    {
        printf("(*'^')q < too few arguments!\n");
        return 0;
    }
    int m, n, k;
    m = n = k = atoi(argv[1]);
    int nproc, nprow, npcol, dims[2], ierror;
    int my_rank, my_row, my_col;
    int icontext;
    int dlen_ = 9;
    char transa = 'N';
    // int m = 1200, n = 800, k = 960;
    int mb = 40, nb = 40;
    // mb = nb = atoi(argv[2]);
    double *a, *b, *c, *dummy;
    int llda, lldb, lldc, ndummy = 2000;
    int lmrow, lncol, lkrow, lkcol;
    int *a_desc = (int *)malloc(dlen_ * sizeof(int));
    int *b_desc = (int *)malloc(dlen_ * sizeof(int));
    int *c_desc = (int *)malloc(dlen_ * sizeof(int));
    int niter = 3;
    int i, ii, jj;
    double t1, t2, etime = 0.0;
    blacs_pinfo_(&my_rank, &nproc);
    dims[0] = dims[1] = 0;
    MPI_Dims_create(nproc, 2, dims);
    if (my_rank == 0)
        printf("dims: %d, %d\n", dims[0], dims[1]);
    nprow = dims[0];
    npcol = dims[1];

    blacs_get_(ADDR(int, 0), ADDR(int, 0), &icontext);
    blacs_gridinit_(&icontext, ADDR(char, 'R'), &nprow, &npcol);
    blacs_gridinfo_(&icontext, &nprow, &npcol, &my_row, &my_col);
    if (my_rank == 2)
    {
        printf("%d %d %d\n", my_rank, my_row, my_col);
    }
    lmrow = min(m, m / nprow) + mb;
    lncol = min(n, n / npcol) + nb;
    lkrow = min(k, k / nprow) + mb;
    lkcol = min(k, k / npcol) + nb;

    if (my_rank == 2)
    {
        printf("%d %d %d %d\n", lmrow, lncol, lkrow, lkcol);
    }
    llda = lmrow + 8;
    lldb = lkrow + 8;
    lldc = lmrow + 8;

    a = (double *)malloc(llda * lkcol * sizeof(double));
    b = (double *)malloc(lldb * lncol * sizeof(double));
    c = (double *)malloc(lldc * lncol * sizeof(double));
    dummy = (double *)malloc(ndummy * ndummy * sizeof(double));

    printf("[%d] %p\n", my_rank, a);

    // dummy = malloc(ndummy * ndummy * sizeof(double));

    descinit_(a_desc, &m, &k, &mb, &nb, ADDR(int, 0), ADDR(int, 0), &icontext, &llda, &ierror);
    descinit_(b_desc, &k, &n, &mb, &nb, ADDR(int, 0), ADDR(int, 0), &icontext, &lldb, &ierror);
    descinit_(c_desc, &m, &n, &mb, &nb, ADDR(int, 0), ADDR(int, 0), &icontext, &lldc, &ierror);
    printf("ierror = %d\n", ierror);

    // initialize A
    for (jj = 1; jj <= k; jj++)
    {
        for (ii = 1; ii <= m; ii++)
        {
            // if (my_rank == 0){
            pdelset_(a, &ii, &jj, a_desc, ADDR(double, (double)(min(ii, jj))));
            //      printf("%d %d\n", jj, ii);
            //       }
        }
    }

    // initialize B
    for (jj = 1; jj <= n; jj++)
    {
        for (ii = 1; ii <= k; ii++)
        {
            //      if (my_rank == 0){
            pdelset_(b, &ii, &jj, b_desc, ADDR(double, 1.0));
            //      printf("%d %d\n", jj, ii);
            //     }
        }
    }

    // initialize C
    for (jj = 1; jj <= n; jj++)
    {
        for (ii = 1; ii <= m; ii++)
        {
            //	    if (my_rank==0){
            pdelset_(c, &ii, &jj, c_desc, ADDR(double, 1.0));
            //      printf("%d %d\n", jj, ii);
            //          }
        }
    }

    for (i = 1; i <= niter; i++)
    {
        if (my_rank == 0)
            printf("iter : %d\n", i);
        t1 = MPI_Wtime();
        // pdgemv_(ADDR(char, 'N'), &m, &n, ADDR(double, 1.0), a, ADDR(int, 1), ADDR(int, 1), a_desc, x, ADDR(int, 1), ADDR(int, 1), x_desc, ADDR(int, 1), ADDR(int, 1), y, ADDR(int, 1), ADDR(int, 1), y_desc, ADDR(int, 1));
        pdgemm_(ADDR(char, 'N'), ADDR(char, 'N'), &m, &n, &k, ADDR(double, 1.0), a, ADDR(int, 1), ADDR(int, 1), a_desc, b, ADDR(int, 1), ADDR(int, 1), b_desc, ADDR(double, 2.0), c, ADDR(int, 1), ADDR(int, 1), c_desc);
        t2 = MPI_Wtime();

        etime = etime + (t2 - t1);
        if (my_rank == 0)
        {
            printf("t1 = %f\n", t1);
            printf("t2 = %f\n", t2);
        }
    }

    ////END/////////
    printf("(*'w')b < node %d complete computation \n", my_rank);
    free(a);
    free(b);
    free(c);

    free(a_desc);
    free(b_desc);
    free(c_desc);
    free(dummy);
    if (my_rank == 0)
    {
        printf("GFLOPS : %f\n", (double)(m) * (double)(n) * (double)(k)*2.0 * (double)(niter) / ((double)(etime)*1.0e9));
    }
    Matrix *mm = create_matrix(nprow, npcol, 100, 100, 20, 20);
    MPI_Finalize();

    return 0;
}