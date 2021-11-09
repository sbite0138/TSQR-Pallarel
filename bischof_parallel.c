#include <assert.h>
// #include <mkl.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h> // for a timing routine.
#include <unistd.h>
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define ADDR(t, v) \
    &(t) { v }
#define measure_time(x)                 \
    do                                  \
    {                                   \
        double start = omp_get_wtime(); \
        {                               \
            x;                          \
        }                               \
        double end = omp_get_wtime();   \
    } while (0)
#define rprintf(...)                                \
    do                                              \
    {                                               \
        if (rank == 0)                              \
        {                                           \
            printf(__VA_ARGS__);                    \
        }                                           \
        blacs_barrier_(&icontext, ADDR(char, 'A')); \
                                                    \
    } while (0)

const int DESC_LEN = 9;
int block_row;
int block_col;

int my_row;
int my_col;
int icontext;

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
    // int icontext;
    // blacs_get_(ADDR(int, 0), ADDR(int, 0), &icontext);
    Matrix *matrix = malloc(sizeof(Matrix));
    matrix->desc = malloc(DESC_LEN * sizeof(int));

    matrix->global_col = global_col;
    matrix->global_row = global_row;
    matrix->local_col = global_col / nproc_col + block_col;
    matrix->local_row = global_row / nproc_row + block_row;
    matrix->leading_dimension = matrix->local_row + 1;
    matrix->data = malloc(matrix->leading_dimension * matrix->local_col * sizeof(double));
    // matrix->data = malloc(1024 * 1024 * sizeof(double));
    int ierror;
    descinit_(matrix->desc, &(matrix->global_row), &(matrix->global_col), ADDR(int, block_row), ADDR(int, block_col), ADDR(int, 0), ADDR(int, 0), &icontext, &(matrix->leading_dimension), &ierror);
    assert(ierror == 0);
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
                printf("%.18f,", val);
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
    printf("aaa\n");
    double *work = malloc(32 * sizeof(double));
    printf("bbb\n");

    // calculate length of work
    int lwork = -1;
    pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info);
    printf("ccc\n");
    lwork = (int)(work[0] + 2);
    printf("lwork %f\n", work[0]);
    // reallocate and compute QR
    free(work);
    work = malloc(lwork * sizeof(double));
    pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info);
    printf("info %d\n", info);
    free(work);
}

void pdgeqrt_wrap(int rank, int proc_row, int proc_col, int m, int n, Matrix *matrix, int row, int col, Matrix *T)
{

    // int icontext;
    // blacs_get_(ADDR(int, 0), ADDR(int, 0), &icontext);
    rprintf("%d %d\n", row, col);
    rprintf("%d %d\n", m, n);

    assert(m >= n);
    assert(T->global_row >= n);
    assert(m + row <= matrix->global_row);
    assert(n + col <= matrix->global_col);

    assert(T->global_col >= n);
    double *tau = malloc((n + col) * sizeof(double));
    rprintf("start\n");
    pdgeqrf_wrap(m, n, matrix, row, col, tau);

    rprintf("ok\n");
    int desc[9];
    desc[0] = 1;
    desc[1] = matrix->desc[1];
    desc[2] = 1;
    desc[3] = n;
    desc[4] = matrix->desc[4];
    desc[5] = matrix->desc[5];
    desc[6] = matrix->desc[6];
    desc[7] = matrix->desc[7];
    desc[8] = 1;
    double val;
    // for (int i = 0; i < n; i++)
    // {
    //     pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, col + 1 + i), desc);
    //     printf("%lf ", val);
    // }
    // printf("\n");
    // printf("%d (%d, %d) tau :", rank, my_row, my_col);
    // for (int i = 0; i < n; i++)
    // {
    //     double val;
    //     pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + i), desc);
    //     printf("%f ", val);

    //     //        printf("%f ", tau[i]);
    // }
    // printf("\n");

    // {
    //     for (int i = 0; i < n; i++)
    //     {
    //         printf("%f ", tau[i]);
    //     }
    //     printf("\n");
    // }

    Matrix *Y = create_matrix(proc_row, proc_col, m, n, block_row, block_col);
    rprintf("start\n");
    for (int i = 0; i < T->global_row; i++)
    {
        for (int j = 0; j < T->global_col; j++)
        {
            set(T, i, j, 0.0);
        }
    }
    rprintf("done\n");

    pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1), desc);
    // set(T, 0, 0, tau[0]);
    set(T, 0, 0, val);
    // print_matrix("T=", T, rank);

    set(Y, 0, 0, 1.0);
    for (int i = 1; i < m; i++)
    {
        set(Y, i, 0, get(matrix, row + i, col));
    }

    Matrix *y = create_matrix(proc_row, proc_col, m, 1, block_row, block_col);

    for (int j = 1; j < n; j++)
    {
        // print_matrix("T=", T, rank);

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

        pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc);
        //  printf("%lf ", val);
        pdgemm_wrap('N', 'N', j, 1, j, -val, T, 0, 0, tmp, 0, 0, 0.0, z, 0, 0);
        // print_matrix("tmp=", tmp, rank);
        // print_matrix("z=", z, rank);

        for (int i = 0; i < m; i++)
        {
            // if (rank == 0)
            // {
            //     printf("loop %d/%d\n", i, m - 1);
            //     printf("%d\n", i);
            // }
            // printf("[%d] start %d\n", rank, i);
            // printf("[%d] done %d\n", rank, i);
            set(Y, i, j, get(y, i, 0));
            // printf("%d\n", icontext);
        }
        for (int i = 0; i < j; i++)
        {
            set(T, i, j, get(z, i, 0));
        }

        pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + j), desc);
        set(T, j, j, val);

        free_matrix(tmp);
        free_matrix(z);
    }
    free_matrix(y);
    free_matrix(Y);
    free(tau);
}

void free_matrix(Matrix *matrix)
{
    free(matrix->data);
    free(matrix->desc);
    free(matrix);
}

void bischof(int rank, int nproc_row, int nproc_col, int N, int L, Matrix *A, Matrix *T, Matrix *Y)
{
    int nb = L;
    int ldt_iter = nb;
    assert(MIN(N, L) >= nb && nb >= 1);
    assert(ldt_iter >= nb);
    assert(N % L == 0);
    //    double *const T_iter = calloc(ldt_iter * L, sizeof(double));
    Matrix *T_iter = create_matrix(nproc_row, nproc_col, L, L, block_col, block_row);
    for (int k = 0; k < N / L - 1; k++)
    {
        // print_matrix("A = ", N, N, A, lda);
        rprintf("iteration %d/%d\n", k + 1, N / L - 1);
        int Nk = N - L - k * L;
        // Aの(k+1,k)ブロック以下をQR分解する
        rprintf("check 1\n");

        pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter);

        ////////////////////////////////!!!!!!!!!!!!!
        blacs_barrier_(&icontext, ADDR(char, 'A'));

        // LAPACK_COL_MAJOR, Nk, L, L, &A[(k + 1) * L + lda * k * L], lda, T_iter, ldt_iter);

        // measure_time(LAPACKE_dgeqrt(LAPACK_COL_MAJOR, Nk, L, L, &A[(k + 1) * L + lda * k * L], lda, T_iter, ldt_iter));
        // print_matrix("T iter = ", L, L, T_iter, ldt_iter);
        // TにT_iterを代入
        // print_matrix("T iter = ", L, L, T_iter, ldt_iter);
        rprintf("check 2\n");
        for (int i = 0; i < L; i++)
        {
            for (int j = 0; j < L; j++)
            {
                if (j >= i)
                { // T[i + (j + L * k) * ldt] = T_iter[i + j * ldt_iter];
                    set(T, i, j + L * k, get(T_iter, i, j));
                }
                else
                // T[i + (j + L * k) * ldt] = T_iter[i + j * ldt_iter] = 0.0;
                {
                    set(T, i, j + L * k, 0.0);
                    set(T_iter, i, j, 0.0);
                }
            }
        };
        //        print_matrix("T = ", L, N, T, ldt);

        // Aを更新する
        rprintf("check 3\n");

        //        measure_time(
        Matrix *V = create_matrix(nproc_row, nproc_col, Nk, L, block_row, block_col);
        // double *V = malloc(Nk * L * sizeof(double));
        Matrix *update_tmp = create_matrix(nproc_row, nproc_col, Nk, Nk, block_row, block_col);
        // double *update_tmp = malloc(Nk * Nk * sizeof(double));
        Matrix *update_P = create_matrix(nproc_row, nproc_col, Nk, L, block_row, block_col);
        // double *update_P = malloc(Nk * L * sizeof(double));
        Matrix *update_beta = create_matrix(nproc_row, nproc_col, L, L, block_row, block_col);
        // double *update_beta = malloc(L * L * sizeof(double));
        Matrix *update_Q = create_matrix(nproc_row, nproc_col, Nk, L, block_row, block_col);
        // double *update_Q = malloc(Nk * L * sizeof(double));
        //	);
        rprintf("check 4\n");

        // Aの(k+1,k)ブロックの先頭アドレス
        // sdouble *a_part = &A[(k + 1) * L + lda * k * L];
        int pad_row = (k + 1) * L;
        int pad_col = k * L;
        // print_matrix("A part = ", Nk, L, a_part, lda);

        // Vにa_partの下三角部分を代入する(LAPACKE_dgeqrtを実行しているので，a_partの下三角部分にはQのcompact-WY表現の一部が入っている)
        for (int i = 0; i < Nk; i++)
        {
            for (int j = 0; j < L; j++)
            {
                if (i == j)
                {
                    set(V, i, j, 1.0);
                    // V[i + j * Nk] = 1.0;
                }
                else if (i > j)
                {
                    set(V, i, j, get(A, pad_row + i, pad_col + j));
                    // V[i + j * Nk] = a_part[i + j * lda];
                    set(A, pad_row + i, pad_col + j, 0.0);
                    // a_part[i + j * lda] = 0.0;
                }
                else
                {
                    set(V, i, j, 0.0);
                    // V[i + j * Nk] = 0.0;
                }
            }
        };
        rprintf("check 5\n");

        // print_matrix("V = ", Nk, L, V, Nk);

        // YにVを代入
        for (int j = 0; j < L; j++)
        {
            for (int i = 0; i < Nk; i++)
            {
                // Y[i + (j + L * k) * ldy] = V[i + j * Nk];
                //  rprintf("%d %d\n", i, j);
                blacs_barrier_(&icontext, ADDR(char, 'A'));
                set(Y, i, j + L * k, get(V, i, j));
                blacs_barrier_(&icontext, ADDR(char, 'A'));
            }
        }
        rprintf("check 6\n");

        // a_partの下三角部分を0クリアし，上三角部分の要素の値を，その要素の位置と対称な位置へ代入する（Aは対称行列なので）
        for (int i = L * k + L; i < N; i++)
        {
            for (int j = L * k; j < L * (k + 1); j++)
            {
                // A[j + i * lda] = A[i + j * lda];
                set(A, j, i, get(A, i, j));
            }
        };
        rprintf("check 7\n");

        // 山本有作先生の『キャッシュマシン向け三重対角化アルゴリズムの性能予測方式』で説明されているBischofのアルゴリズム中のPを構築する
        // measure_time(
        // cblas_dsymm(CblasColMajor, CblasLeft, CblasLower, Nk, L, 1.0, &A[(k + 1) * L + lda * (k + 1) * L], lda, V, Nk, 0.0, update_tmp, Nk);
        pdsymm_wrap('L', 'L', Nk, L, 1.0, A, (k + 1) * L, (k + 1) * L, V, 0, 0, 0.0, update_tmp, 0, 0);
        // cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Nk, L, L, 1.0, update_tmp, Nk, T_iter, ldt_iter, 0.0, update_P, Nk);
        pdgemm_wrap('N', 'N', Nk, L, L, 1.0, update_tmp, 0, 0, T_iter, 0, 0, 0.0, update_P, 0, 0);
        //);
        rprintf("check 8\n");

        // betaを構築する
        // measure_time(
        pdgemm_wrap('T', 'N', L, L, Nk, 0.5, V, 0, 0, update_P, 0, 0, 0.0, update_tmp, 0, 0);
        // cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, L, L, Nk, 0.5, V, Nk, update_P, Nk, 0.0, update_tmp, Nk);
        pdgemm_wrap('T', 'N', L, L, L, 1.0, T_iter, 0, 0, update_tmp, 0, 0, 0.0, update_beta, 0, 0);
        rprintf("check 9\n");
        // cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, L, L, L, 1.0, T_iter, ldt_iter, update_tmp, Nk, 0.0, update_beta, L);
        //);
        // Qを構築する
        // measure_time(
        for (int i = 0; i < Nk; i++)
        {
            for (int j = 0; j < L; j++)
            {
                set(update_Q, i, j, get(update_P, i, j));
                // update_Q[i + j * Nk] = update_P[i + j * Nk];
            }
        }
        //);
        rprintf("check 10\n");

        // measure_time(
        pdgemm_wrap('N', 'N', Nk, L, L, -1.0, V, 0, 0, update_beta, 0, 0, 1.0, update_Q, 0, 0);
        // cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Nk, L, L, -1.0, V, Nk, update_beta, L, 1.0, update_Q, Nk);
        //);
        //  print_matrix("update_Q", Nk, L, update_Q, Nk);
        rprintf("check 11\n");

        // Aをrank-2k更新する
        // measure_time(
        pdsyr2k_wrap('L', 'N', Nk, L, -1.0, V, 0, 0, update_Q, 0, 0, 1.0, A, (k + 1) * L, (k + 1) * L);
        // cblas_dsyr2k(CblasColMajor, CblasLower, CblasNoTrans, Nk, L, -1.0, V, Nk, update_Q, Nk, 1.0, &A[(k + 1) * L + lda * (k + 1) * L], lda);
        //);
        rprintf("check 12\n");
        blacs_barrier_(&icontext, ADDR(char, 'A'));

        free_matrix(V);
        free_matrix(update_P);
        free_matrix(update_beta);
        free_matrix(update_Q);
        free_matrix(update_tmp);
        rprintf("check 13\n");
    }
    free_matrix(T_iter);
    for (int i = 0; i < N; i++)
    {
        for (int j = i; j < N; j++)
        {
            set(A, i, j, get(A, j, i));
            // A[i + j * lda] = A[j + i * lda];
        }
    }
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
        printf("Usage %s matrix_size\n", argv[0]);
        return 0;
    }

    int m, n, k;
    m = atoi(argv[1]);
    //  n = atoi(argv[2]);
    block_row = 64; // atoi(argv[3]);
    block_col = 64; // atoi(argv[4]);

    int nproc, nproc_row, nproc_col, dims[2], ierror;
    int L = 2;
    int rank;

    // int icontext;
    // int m = 1200, n = 800, k = 960;
    blacs_pinfo_(&rank, &nproc);
    dims[0] = dims[1] = 0;

    MPI_Dims_create(nproc, 2, dims);
    nproc_row = dims[0];
    nproc_col = dims[1];

    printf("nproc_row: %d nproc_col: %d\n", nproc_row, nproc_col);
    // sl_init_(&icontext, &nproc_row, &nproc_col);
    blacs_get_(ADDR(int, 0), ADDR(int, 0), &icontext);
    blacs_gridinit_(&icontext, ADDR(char, 'R'), &nproc_row, &nproc_col);
    blacs_gridinfo_(&icontext, &nproc_row, &nproc_col, &my_row, &my_col);

    Matrix *A = create_matrix(nproc_row, nproc_col, m, m, block_row, block_col);
    Matrix *T = create_matrix(nproc_row, nproc_col, L, m, block_row, block_col);
    Matrix *Y = create_matrix(nproc_row, nproc_col, m, m, block_row, block_col);
    for (size_t i = 0; i < m; ++i)
    {
        for (size_t j = 0; j < m; ++j)
        {
            double r = (double)(rand()) / RAND_MAX;
            set(A, i, j, r);
            set(A, j, i, r);
        }
    }

    bischof(rank, nproc_row, nproc_col, m, L, A, T, Y);
    rprintf("free\n");
    print_matrix("A=,", A, rank);
    free_matrix(A);
    free_matrix(T);
    free_matrix(Y);
    MPI_Finalize();
    // MPI_Finalize();
    // if (rank == 0)
    //     printf("a = np.matrix(a)\nR = np.matrix(R)\nY = np.matrix(Y)\nT = np.matrix(T)\ntmp = -Y*T*Y.T\nfor i in range(len(tmp)):\n    tmp[i, i] += 1.0\ntmp = tmp*R\ntmp = tmp-a\nprint(\"error:\", np.linalg.norm(tmp, 2)/np.linalg.norm(a, 2))");

    return 0;
}