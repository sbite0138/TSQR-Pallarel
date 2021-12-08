#include <assert.h>
// #include <mkl.h>
#include <cblas.h>
#include <lapacke.h>
#include <mpi.h>
#include <omp.h> // for a timing routine.
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define ADDR(t, v) \
    &(t) { v }
#define rprintf(...)                                \
    do                                              \
    {                                               \
        if (rank == 0)                              \
        {                                           \
            sbprintf(__VA_ARGS__);                  \
        }                                           \
        blacs_barrier_(&icontext, ADDR(char, 'A')); \
                                                    \
    } while (0)

#define measure_time(x)                                                                                                 \
    do                                                                                                                  \
    {                                                                                                                   \
        double start = omp_get_wtime();                                                                                 \
        {                                                                                                               \
            x;                                                                                                          \
        }                                                                                                               \
        double end = omp_get_wtime();                                                                                   \
        if (print_checkcode == false)                                                                                   \
            rprintf("@ {\"rank\":%d, \"line\":%d, \"cmd\":\"%s\", \"time\":%.18f}\n", rank, __LINE__, #x, end - start); \
    } while (0)

typedef struct StringBuffer
{
    char *buf;
    size_t size;
    size_t next;
} StringBuffer;
StringBuffer sBuffer;

void initBuffer()
{
    StringBuffer *sb = &sBuffer;
    sb->size = 65536;
    sb->buf = (char *)malloc(sb->size * sizeof(char));
    sb->next = 0;
}

void sbprintf(const char *fmt, ...)
{
    va_list args;
    va_start(args, fmt);

    sBuffer.next += vsprintf(&sBuffer.buf[sBuffer.next], fmt, args);
    va_end(args);
    if (sBuffer.next > 8 * sBuffer.size / 10)
    {
        sBuffer.size *= 2;
        sBuffer.buf = realloc(sBuffer.buf, sBuffer.size * sizeof(char));
        fprintf(stderr, "realloced!\n");
    }
}
void dumpBuffer()
{
    sBuffer.buf[sBuffer.next++] = '\0';
    puts(sBuffer.buf);
}

const size_t DESC_LEN = 9;
int block_row;
int block_col;

int my_row;
int my_col;
int proc_num;
int proc_row_num;
int proc_col_num;
int icontext;
int rank;
bool print_checkcode = true;

typedef struct
{
    double *data;
    int *desc;
    size_t global_row;
    size_t global_col;
    size_t local_row;
    size_t local_col;
    size_t leading_dimension;
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

void free_matrix(Matrix *matrix)
{
    free(matrix->data);
    free(matrix->desc);
    free(matrix);
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
                printf("%12.8lf,", val);
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

void pdgemr2d_wrap(int m, int n, Matrix *A, int row_a, int col_a, Matrix *B, int row_b, int col_b)
{
    row_a++;
    col_a++;
    row_b++;
    col_b++;
    pdgemr2d_(&m, &n, A->data, &row_a, &col_a, A->desc, B->data, &row_b, &col_b, B->desc, &icontext);
}

void pdlaset_wrap(char uplo, int m, int n, double alpha, double beta, Matrix *A, int row, int col)
{
    row++;
    col++;
    pdlaset_(&uplo, &m, &n, &alpha, &beta, A->data, &row, &col, A->desc);
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
    double *work = malloc(32 * sizeof(double));

    // calculate length of work
    int lwork = -1;
    (pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info));
    lwork = (int)(work[0] + 2);
    // reallocate and compute QR
    free(work);
    work = malloc(lwork * sizeof(double));
    (pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info));
    free(work);
}

void pdgeqrt_wrap(int rank, int proc_row, int proc_col, int m, int n, Matrix *matrix, int row, int col, Matrix *T)
{

    // int icontext;
    // blacs_get_(ADDR(int, 0), ADDR(int, 0), &icontext);
    assert(m >= n);
    assert(T->global_row >= n);
    assert(m + row <= matrix->global_row);
    assert(n + col <= matrix->global_col);

    assert(T->global_col >= n);
    double *tau;
    measure_time(
        tau = malloc((n + col) * sizeof(double)));
    measure_time(pdgeqrf_wrap(m, n, matrix, row, col, tau));
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
    Matrix *Y;
    measure_time(
        Y = create_matrix(proc_row, proc_col, m, n, block_row, block_col););
    measure_time(
        pdlaset_wrap('A', T->global_row, T->global_col, 0.0, 0.0, T, 0, 0));
    measure_time(pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col), desc));

    // set(T, 0, 0, tau[0]);
    set(T, 0, 0, val);
    // print_matrix("T=", T, rank);
    measure_time(pdgemr2d_wrap(m, n, matrix, row, col, Y, 0, 0));
    measure_time(pdlaset_wrap('U', Y->global_row, Y->global_col, 0.0, 1.0, Y, 0, 0));
    Matrix *y;
    Matrix *z;
    Matrix *tmp;

    measure_time(y = create_matrix(proc_row, proc_col, m, n, block_row, block_col);
                 z = create_matrix(proc_row, proc_col, n, 1, block_row, block_col);
                 tmp = create_matrix(proc_row, proc_col, n, 1, block_row, block_col));

    measure_time(pdgemr2d_wrap(m, n, matrix, row, col, y, 0, 0));
    measure_time(pdlaset_wrap('U', y->global_row, y->global_col, 0.0, 1.0, y, 0, 0));
    for (int j = 1; j < n; j++)
    {
        // (pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0));
        // (set(y, j, 0, 1.0));
        // (pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0));
        // // measure_time(pdgemm_wrap('T', 'N', j, 1, m, 1.0, Y, 0, 0, y, 0, 0, 0.0, tmp, 0, 0));
        measure_time(pdgemv_wrap('T', m, j, 1.0, Y, 0, 0, y, 0, j, 1, 0.0, tmp, 0, 0, 1));
        measure_time(pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc));
        // printf("%lf ", val);
        measure_time(pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, T, 0, j, 1));
        // measure_time(pdgemm_wrap('N', 'N', j, 1, j, -val, T, 0, 0, tmp, 0, 0, 0.0, z, 0, 0));

        // print_matrix("tmp=", tmp, rank);
        // print_matrix("z=", z, rank);
        // (pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j));
        // measure_time(pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j));
        measure_time(pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc));
        // val=3.14;
        measure_time(set(T, j, j, val));
    }
    measure_time(
        free_matrix(tmp);
        free_matrix(z);
        free_matrix(y);
        free_matrix(Y);
        free(tau););
}

void bischof(int rank, int nproc_row, int nproc_col, int N, int L, Matrix *A, Matrix *T, Matrix *Y)
{
    int nb = L;
    int ldt_iter = nb;
    assert(MIN(N, L) >= nb && nb >= 1);
    assert(ldt_iter >= nb);
    assert(N % L == 0);
    Matrix *T_iter = create_matrix(nproc_row, nproc_col, L, L, block_col, block_row);
    Matrix *V = create_matrix(nproc_row, nproc_col, N - L, L, block_row, block_col);
    Matrix *update_tmp = create_matrix(nproc_row, nproc_col, N - L, N - L, block_row, block_col);
    Matrix *update_P = create_matrix(nproc_row, nproc_col, N - L, L, block_row, block_col);
    Matrix *update_beta = create_matrix(nproc_row, nproc_col, L, L, block_row, block_col);
    Matrix *update_Q = create_matrix(nproc_row, nproc_col, N - L, L, block_row, block_col);

    for (int k = 0; k < N / L - 1; k++)
    {
        //  rprintf("iteration %d/%d\n", k + 1, N / L - 1);
        int Nk = N - L - k * L;
        // Aの(k+1,k)ブロック以下をQR分解する
        measure_time(pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter));
        measure_time(pdlaset_wrap('L', L, L, 0.0, 0.0, T, 0, L * k));

        // print_matrix("A part = ", Nk, L, a_part, lda);
        int pad_row = (k + 1) * L;
        int pad_col = k * L;

        // VにA_partの下三角部分を代入する(LAPACKE_dgeqrtを実行しているので，a_partの下三角部分にはQのcompact-WY表現の一部が入っている)
        measure_time(pdgemr2d_wrap(Nk, L, A, pad_row, pad_col, V, 0, 0));
        measure_time(pdlaset_wrap('U', Nk, L, 0.0, 1.0, V, 0, 0));
        // print_matrix("V = ", Nk, L, V, Nk);

        // YにVを代入
        measure_time(pdgemr2d_wrap(Nk, L, V, 0, 0, Y, 0, L * k));

        // print_matrix("Y=", Y, rank);

        // 山本有作先生の『キャッシュマシン向け三重対角化アルゴリズムの性能予測方式』で説明されているBischofのアルゴリズム中のPを構築する
        measure_time(pdsymm_wrap('L', 'L', Nk, L, 1.0, A, (k + 1) * L, (k + 1) * L, V, 0, 0, 0.0, update_tmp, 0, 0));
        // print_matrix("T_iter", T_iter, rank);

        measure_time(pdgemm_wrap('N', 'N', Nk, L, L, 1.0, update_tmp, 0, 0, T_iter, 0, 0, 0.0, update_P, 0, 0));
        // print_matrix("update_P", update_P, rank);

        // betaを構築する
        measure_time(pdgemm_wrap('T', 'N', L, L, Nk, 0.5, V, 0, 0, update_P, 0, 0, 0.0, update_tmp, 0, 0));
        measure_time(pdgemm_wrap('T', 'N', L, L, L, 1.0, T_iter, 0, 0, update_tmp, 0, 0, 0.0, update_beta, 0, 0));
        // print_matrix("update_P", update_beta, rank);

        //  Qを構築する
        measure_time(pdgemr2d_wrap(Nk, L, update_P, 0, 0, update_Q, 0, 0));
        // print_matrix("update_Q=", update_Q, rank);

        measure_time(pdgemm_wrap('N', 'N', Nk, L, L, -1.0, V, 0, 0, update_beta, 0, 0, 1.0, update_Q, 0, 0));
        // print_matrix("update_Q", Nk, L, update_Q, Nk);

        // Aをrank-2k更新する
        // print_matrix("V=", V, rank);
        // print_matrix("update_Q=", update_Q, rank);

        measure_time(pdsyr2k_wrap('L', 'N', Nk, L, -1.0, V, 0, 0, update_Q, 0, 0, 1.0, A, (k + 1) * L, (k + 1) * L));
    }
    free_matrix(V);
    free_matrix(update_P);
    free_matrix(update_beta);
    free_matrix(update_Q);
    free_matrix(update_tmp);
    free_matrix(T_iter);
}

bool is_power_of_2(int n)
{
    while (n > 1)
    {
        n /= 2;
    }
    return n == 1;
}

void TSQR_init(int proc_row_id, int proc_col_id, int m, int n, Matrix *matrix, int row, int col, double **data, int *id)
{
    assert(is_power_of_2(proc_num));
    assert(is_power_of_2(m));
    assert(m % proc_num == 0);
    assert(m >= proc_num);

    int m_part = m / proc_num;
    int n_part = n;
    printf("n = %d", n);
    Matrix matrixes[proc_num];
    // printf("\n%d\n", proc_num);

    // Ex: for 2x3 grid, mapping [row_id, col_id] => pos is like:
    // 0 1 2
    // 3 4 5
    int my_grid_pos = proc_row_id * proc_col_num + proc_col_id;
    *id = my_grid_pos;
    int cnt = 0;
    for (int i = 0; i < proc_row_num; i++)
    {
        // for (int j = 0; j < proc_col_num; j++)
        for (int j = 0; j < proc_col_num; j++)
        {
            // printf("%d %d\n", i, j);

            // Ex: for 2x3 grid, mapping [row_id, col_id] => pos is like:
            // 0 1 2
            // 3 4 5
            int current_grid_pos = i * proc_col_num + j;

            // init each proc's matrix
            matrixes[current_grid_pos].global_row = m_part;
            matrixes[current_grid_pos].global_col = n_part;
            // scalapackの仕様的にlocal_col/row>=block_col/rowである必要があるか？
            // なければ以下の最後の + block_col/row は不要
            matrixes[current_grid_pos].local_row = m_part; //+ block_row;
            matrixes[current_grid_pos].local_col = n_part; // + block_col;

            matrixes[current_grid_pos].leading_dimension = matrixes[current_grid_pos].local_row;
            if (my_grid_pos == current_grid_pos)
            {
                // printf("[%d] %d ok\n", rank, icontext);
                matrixes[current_grid_pos].data = malloc(matrixes[current_grid_pos].leading_dimension * matrixes[current_grid_pos].local_col * sizeof(double));
                *data = matrixes[current_grid_pos].data;
            }
            matrixes[current_grid_pos].desc = malloc(DESC_LEN * sizeof(int));

            int ierror;
            descinit_(matrixes[current_grid_pos].desc, &(matrixes[current_grid_pos].global_row),
                      &(matrixes[current_grid_pos].global_col), ADDR(int, matrixes[current_grid_pos].global_row), ADDR(int, matrixes[current_grid_pos].global_col), ADDR(int, i), ADDR(int, j), &icontext, &(matrixes[current_grid_pos].leading_dimension), &ierror);

            assert(ierror == 0);
            pdgemr2d_wrap(matrixes[current_grid_pos].global_row, matrixes[current_grid_pos].global_col, matrix, row + matrixes[current_grid_pos].global_row * cnt, col, &(matrixes[current_grid_pos]), 0, 0);
            if (my_grid_pos == current_grid_pos)
            {
                printf("desc %d\n", matrixes[current_grid_pos].desc[0]);
                printf("desc %d\n", matrixes[current_grid_pos].desc[1]);
                printf("desc %d\n", matrixes[current_grid_pos].desc[2]);
                printf("desc %d\n", matrixes[current_grid_pos].desc[3]);
                printf("desc %d\n", matrixes[current_grid_pos].desc[4]);
                printf("desc %d\n", matrixes[current_grid_pos].desc[5]);
                printf("desc %d\n", matrixes[current_grid_pos].desc[6]);
                printf("desc %d\n", matrixes[current_grid_pos].desc[7]);
                printf("desc %d\n", matrixes[current_grid_pos].desc[8]);

                printf("%lf ", (*data)[0]);
                printf("%lf ", (*data)[1]);
                printf("%lf\n", (*data)[2]);
            }
            cnt++;
        }
    }
}

void TSQR(int rank, int proc_row_id, int proc_col_id, int m, int n, Matrix *matrix, int row, int col, Matrix *T)
{
    int id;
    double *data;
    TSQR_init(proc_row_id, proc_col_id, m, n, matrix, row, col, &data, &id);
    blacs_barrier_(&icontext, ADDR(char, 'A'));
    for (int k = 0; k < proc_num; k++)
    {
        if (k == id)
        {
            for (int i = 0; i < m / proc_num; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    printf("%lf ", data[j * (m / proc_num) + i]);
                    // printf("%d ", j * m / proc_num + i);
                }
                printf("\n");
            }
            printf("\n");
        }
        blacs_barrier_(&icontext, ADDR(char, 'A'));
    }
    int m_part = m / proc_num;
    int n_part = n;
    double *tau = malloc(n * sizeof(double));
    double *R = calloc(n_part * n_part, sizeof(double));
    int current_blocknum = 1;
    LAPACKE_dgeqrf(LAPACK_COL_MAJOR, m_part, n_part, data, m_part, tau);
    LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'U', n_part, n_part, data, m_part, R, n_part);
    int k = 1;
    int i = id;
    int p = proc_num;
    while ((1 << k) < proc_num)
    {
        printf("%d\n", k);
        if (i % (1 << k) == 0 && (i + (1 << (k - 1))) < p)
        {
            int j = i + (1 << (k - 1));
            printf("[%d] send to %d\n", i, j);
            double *R_tmp = calloc(2 * n_part * n_part, sizeof(double));
            double *R_rcv = calloc(n_part * n_part, sizeof(double));
            MPI_Status st;
            int ret = MPI_Recv(R_rcv, n_part * n_part, MPI_DOUBLE, j, 0, MPI_COMM_WORLD, &st);
            assert(ret == MPI_SUCCESS);
            LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'U', n_part, n_part, R, n_part, R_tmp, 2 * n_part);
            LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'U', n_part, n_part, R_rcv, n_part, &R_tmp[n_part], 2 * n_part);
            LAPACKE_dgeqrf(LAPACK_COL_MAJOR, 2 * n_part, n_part, R_tmp, 2 * n_part, tau);
            LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'U', n_part, n_part, R_tmp, 2 * n_part, R, n_part);
            free(R_tmp);
            free(R_rcv);
        }
        else if (i % (1 << k) == (1 << (k - 1)))
        {
            printf("[%d] send to %d\n", i, (1 << (k - 1)));
            int ret = MPI_Send(R, n_part * n_part, MPI_DOUBLE, i - (1 << (k - 1)), 0, MPI_COMM_WORLD);
            assert(ret == MPI_SUCCESS);
        }

        k++;
        blacs_barrier_(&icontext, ADDR(char, 'A'));
    }
    if (id == 0)
    {
        for (int i = 0; i < n_part; i++)
        {
            for (int j = 0; j < n_part; j++)
            {
                printf("%lf ", R[i * n_part + j]);
            }
            printf("\n");
        }
    }
}

int main(int argc, char **argv)
{
    // init
    MPI_Init(&argc, &argv);
    if (argc != 2)
    {
        printf("Usage %s matrix_size\n", argv[0]);
        return 0;
    }
    printf("do nothing");
    initBuffer();
    printf("\n");
    int m, n, k;
    m = atoi(argv[1]);
    // n = atoi(argv[2]);
    block_row = 60;        // atoi(argv[3]);
    block_col = block_row; // atoi(argv[4]);
    unsigned long const random_seed = 100;
    srand(random_seed);

    int nproc, dims[2], ierror;
    int L = 32;

    blacs_pinfo_(&rank, &nproc);
    proc_num = nproc;
    dims[0] = dims[1] = 0;

    MPI_Dims_create(nproc, 2, dims);
    proc_row_num = dims[0];
    proc_col_num = dims[1];

    blacs_get_(ADDR(int, 0), ADDR(int, 0), &icontext);
    blacs_gridinit_(&icontext, ADDR(char, 'R'), &proc_row_num, &proc_col_num);
    blacs_gridinfo_(&icontext, &proc_row_num, &proc_col_num, &my_row, &my_col);
    n = 10;
    // main calc
    printf("%d %d", m, n);

    Matrix *A = create_matrix(proc_row_num, proc_col_num, m, n, block_row, block_col);
    Matrix *T = create_matrix(proc_row_num, proc_col_num, n, n, block_row, block_col);

    measure_time(for (size_t i = 0; i < A->global_row; ++i) {
        for (size_t j = 0; j < A->global_col; ++j)
        {
            // double r = i * A->global_col + j;
            double r = (double)(rand()) / RAND_MAX;
            set(A, i, j, r);
            // set(A, j, i, r);
        }
    });
    blacs_barrier_(&icontext, ADDR(char, 'A'));
    print_matrix("A=", A, rank);
    printf("start\n");
    TSQR(rank, my_row, my_col, m, n, A, 0, 0, T);
    free_matrix(A);

    MPI_Finalize();
    if (rank == 0)
        dumpBuffer();
    return 0;
}
