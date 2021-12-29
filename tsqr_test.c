#include <assert.h>
//#include <mkl.h>
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
#define rprintf(...)                                   \
    do                                                 \
    {                                                  \
        if (rank == 0)                                 \
        {                                              \
            sbprintf(__VA_ARGS__);                     \
        }                                              \
        blacs_barrier_(&icontext_2d, ADDR(char, 'A')); \
                                                       \
    } while (0)

#define measure_time(x)                                                                                                                     \
    do                                                                                                                                      \
    {                                                                                                                                       \
        level++;                                                                                                                            \
        double start = omp_get_wtime();                                                                                                     \
        {                                                                                                                                   \
            x;                                                                                                                              \
        }                                                                                                                                   \
        double end = omp_get_wtime();                                                                                                       \
        level--;                                                                                                                            \
        if (print_checkcode == false)                                                                                                       \
            rprintf("@ {\"rank\":%d,\"level\":%d, \"line\":%d, \"cmd\":\"%s\", \"time\":%.18f}\n", rank, level, __LINE__, #x, end - start); \
    } while (0)
#define append(vec, vec_size, val)                           \
    do                                                       \
    {                                                        \
        vec = realloc(vec, (vec_size + 1) * sizeof(size_t)); \
        (vec)[vec_size] = val;                               \
        vec_size++;                                          \
    } while (false)
int level = 0;
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
int icontext_2d;
int icontext_1d;
bool print_checkcode = false;

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

typedef struct
{
    double *data;
    size_t row;
    size_t col;
    size_t leading_dimension;
} LocalMatrix;

LocalMatrix *create_local_matrix(size_t row, size_t col)
{
    LocalMatrix *lmat = malloc(sizeof(LocalMatrix));
    lmat->row = row;
    lmat->col = col;
    lmat->leading_dimension = row;
    lmat->data = malloc(lmat->leading_dimension * lmat->col * sizeof(double));
    return lmat;
}

void free_local_matrix(LocalMatrix *lmat)
{
    free(lmat->data);
    free(lmat);
}

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
    matrix->data = malloc(matrix->leading_dimension * matrix->local_col * (size_t)sizeof(double));
    // matrix->data = malloc(1024 * 1024 * sizeof(double));
    int ierror;
    descinit_(matrix->desc, &(matrix->global_row), &(matrix->global_col), ADDR(int, block_row), ADDR(int, block_col), ADDR(int, 0), ADDR(int, 0), &icontext_2d, &(matrix->leading_dimension), &ierror);
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
                printf("%.18lf,", val);
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
    pdgemr2d_(&m, &n, A->data, &row_a, &col_a, A->desc, B->data, &row_b, &col_b, B->desc, &icontext_2d);
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
    size_t lwork = -1;
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
        tau = malloc(((size_t)n + (size_t)col) * sizeof(double)));
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

bool can_TSQR(int m, int n)
{
    int block_size = m / proc_num;
    printf("block_size %d, m %d, proc_num %d, m mod proc_num %d\n", block_size, m, proc_num, m % proc_num);
    // rprintf("m %d n %d proc_num %d\n", m, n, proc_num)

    return (m >= proc_num * n) &&
           (m % proc_num <= block_size);
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
        if (can_TSQR(Nk, L))
        {
            rprintf("# TSQR!\n");
            measure_time(TSQR_HR(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter));
        }
        else
        {
            measure_time(pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter));
        }
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
        if (n % 2 != 0)
            return false;
        n /= 2;
    }
    return n == 1;
}

void exchange(int m_total, int n, int m_part, int m_last, int block_num, int have1, int have2, int need1, int need2, Matrix *mat, int rank)
{
    int m_head = m_part;
    int m_tail = m_total - m_head;
    assert(need1 != -1);
    assert(have1 != -1);
    int t = proc_num - (block_num % proc_num);
    printf("%d\n", m_last);
    int recv1 = -1;
    int recv2 = -1;
    int send1 = -1;
    int send2 = -1;
    // set recv
    recv1 = need1 % proc_num;
    if (need2 != -1)
    {
        recv2 = need2 % proc_num;
    }
    send1 = have1;
    if (send1 >= t)
        send1 = (have1 - t) / 2 + t;
    if (have2 != -1)
    {
        send2 = have2;
        if (send2 >= t)
            send2 = (have2 - t) / 2 + t;
    }
    printf("status [%d] %d %d %d %d\n", rank, send1, send2, recv1, recv2);
    printf("send   [%d] %d %d\n", rank, m_head, m_tail);

    MPI_Request reqs[4];

    int reqs_idx = 0;
    // send 1
    double *data_send1 = malloc(m_head * n * sizeof(double));
    LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', m_head, n, mat->data, m_total, data_send1, m_head);

    MPI_Isend(data_send1, m_head * n, MPI_DOUBLE, send1, have1, MPI_COMM_WORLD, &reqs[reqs_idx++]);
    //  recv 1
    int m_recv1 = (need1 == block_num - 1) ? m_last : m_part;
    double *data_recv1 = malloc(m_recv1 * n * sizeof(double));
    MPI_Irecv(data_recv1, m_recv1 * n, MPI_DOUBLE, recv1, need1, MPI_COMM_WORLD, &reqs[reqs_idx++]);

    // send 2
    if (send2 != -1)
    {
        double *data_send2 = malloc(m_tail * n * sizeof(double));
        LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', m_tail, n, mat->data, m_total, data_send2, m_tail);
        MPI_Isend(data_send2, m_tail * n, MPI_DOUBLE, send2, have2, MPI_COMM_WORLD, &reqs[reqs_idx++]);
    }

    // recv 2
    if (recv2 != -1)
    {
        int m_recv2 = (need2 == block_num - 1) ? m_last : m_part;
        double *data_recv2 = malloc(m_recv2 * n * sizeof(double));
        MPI_Irecv(data_recv2, m_recv2 * n, MPI_DOUBLE, recv2, need2, MPI_COMM_WORLD, &reqs[reqs_idx++]);
        printf("recv   [%d] %d \n", rank, m_recv2);
        printf("mlast   [%d] %d\n", rank, (need2 == block_num - 1) ? m_last : m_part);
    }

    for (int i = 0; i < reqs_idx; i++)
    {

        MPI_Status st;
        MPI_Wait(&reqs[i], &st);
    }
}
void TSQR_init(int rank, int m, int n, Matrix *matrix, int row, int col, double **data)
{

    assert(is_power_of_2(proc_num));
    // assert(is_power_of_2(m));
    // assert(m % proc_num == 0);
    assert(m >= proc_num * n);

    size_t m_part = m / proc_num;
    size_t block_num = (m + m_part - 1) / m_part;
    int have1 = -1;
    int have2 = -1;
    int need1 = -1;
    int need2 = -1;

    if (rank < block_num % proc_num)
    {
        have1 = rank;
        have2 = proc_num + rank;
    }
    else
    {
        have1 = rank;
    }
    int t = proc_num - (block_num % proc_num);
    if (rank < t)
    {
        need1 = rank;
    }
    else
    {
        need1 = 2 * (rank - t) + t;
        need2 = need1 + 1;
    }
    int n_part = n;
    int desc[DESC_LEN];
    int ierror;
    int pad = 0;
    if (m % proc_num != 0)
    {
        if (rank == 0)
        {
            pad = m % proc_num;
        }
    }
    int m_total = numroc_(&m, &m_part, &rank, ADDR(int, 0), &proc_num);
    // printf("[%d] descinit 3 %d %d %d %d %d %d\n", rank, m, n, m_part, n_part, icontext_1d, m_total);
    descinit_(desc, ADDR(int, m), ADDR(int, n), &m_part, &n, ADDR(int, 0), ADDR(int, 0), &icontext_1d, &m_total, &ierror);
    printf("[%d] done 3 %d %d %d %d %d %d\n", rank, m, n, m_part, n_part, icontext_1d, m_total);
    Matrix mat;
    mat.desc = desc;
    mat.data = malloc(((size_t)m_total) * (size_t)n_part * sizeof(double));
    printf("[%d] %d %d %d %d\n", rank, m, n, row, col);
    pdgemr2d_wrap(m, n, matrix, row, col, &mat, 0, 0);
    *data = mat.data;
    int m_last;
    if (m % m_part == 0)
    {
        m_last = m_part;
    }
    else
    {
        m_last = m % m_part;
    }
    exchange(m_total, n, m_part, m_last, block_num, have1, have2, need1, need2, &mat, rank);
    MPI_Barrier(MPI_COMM_WORLD);
}

// TODO:冗長な二重ポインタを取り除く
// TODO:引数にm,nではなくm_part,n_partを与えるようにする（見通しが良くなりそうなので）
void TSQR(int id, int m_part, int n_part, double **data, double **Y, size_t **Y_heads, double **R, double **tau)
{
    MPI_Barrier(MPI_COMM_WORLD);
    // for (int k = 0; k < proc_num; k++)
    // {
    //     if (k == id)
    //     {
    //         for (int i = 0; i < m / proc_num; i++)
    //         {
    //             for (int j = 0; j < n; j++)
    //             {
    //                 printf("%lf ", (*data)[j * (m / proc_num) + i]);
    //                 // printf("%d ", j * m / proc_num + i);
    //             }
    //             printf("\n");
    //         }
    //         printf("\n");
    //     }
    //     blacs_barrier_(&icontext, ADDR(char, 'A'));
    // }
    //    printf("m %d n %d proc_num %d\n", m, n, proc_num);
    //   printf("m_part %d n_part %d\n", m_part, n_part);
    assert(m_part >= n_part);
    *tau = malloc((size_t)n_part * sizeof(double));
    size_t tau_size = n_part;
    *R = calloc(n_part * n_part, sizeof(double));
    *Y = calloc(m_part * n_part, sizeof(double));
    size_t Y_heads_size = 0;
    *Y_heads = NULL;
    append(*Y_heads, Y_heads_size, 0);
    size_t Y_size = m_part * n_part;

    int current_blocknum = 1;
    LAPACKE_dgeqrf(LAPACK_COL_MAJOR, m_part, n_part, *data, m_part, *tau);
    LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'U', n_part, n_part, *data, m_part, *R, n_part);
    LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'L', m_part, n_part, *data, m_part, *Y, m_part);
    LAPACKE_dlaset(LAPACK_COL_MAJOR, 'U', m_part, n_part, 0.0, 1.0, *Y, m_part);
    int k = 1;
    int i = id;
    int p = proc_num;
    // printf("[%d] initial R: %lf\n", i, (*R)[0]);
    while ((1 << k) <= proc_num)
    {
        // printf("%d\n", k);
        if (i % (1 << k) == 0 && (i + (1 << (k - 1))) < p)
        {
            int j = i + (1 << (k - 1));
            double *R_tmp = calloc(2 * n_part * n_part, sizeof(double));
            double *R_rsv = calloc(n_part * n_part, sizeof(double));
            MPI_Status st;
            int ret = MPI_Recv(R_rsv, n_part * n_part, MPI_DOUBLE, j, 0, MPI_COMM_WORLD, &st);
            assert(ret == MPI_SUCCESS);
            // printf("[%d] recv from %d(%lf)\n", i, j, R_rsv[0]);
            //  write to Y

            *Y = realloc(*Y, (Y_size + n_part * n_part + n_part * n_part) * sizeof(double));
            LAPACKE_dlaset(LAPACK_COL_MAJOR, 'A', 2 * n_part, n_part, 0.0, 0.0, &((*Y)[Y_size]), 2 * n_part);

            append(*Y_heads, Y_heads_size, Y_size);

            *tau = realloc(*tau, (tau_size + n_part) * sizeof(double));
            // write to R
            LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'U', n_part, n_part, *R, n_part, &((*Y)[Y_size]), 2 * n_part);
            LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'U', n_part, n_part, R_rsv, n_part, &((*Y)[Y_size + n_part]), 2 * n_part);

            ret = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, 2 * n_part, n_part, &((*Y)[Y_size]), 2 * n_part, &((*tau)[tau_size]));
            // printf("[%d] ret = %d\n", i, ret);
            LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'U', n_part, n_part, &((*Y)[Y_size]), 2 * n_part, *R, n_part);
            // printf("[%d] R %lf\n", i, (*R)[0]);
            Y_size += n_part * n_part * 2;
            tau_size += n_part;
            free(R_tmp);
            free(R_rsv);
        }
        else if (i % (1 << k) == (1 << (k - 1)))
        {
            int ret = MPI_Send(*R, n_part * n_part, MPI_DOUBLE, i - (1 << (k - 1)), 0, MPI_COMM_WORLD);
            assert(ret == MPI_SUCCESS);
            // printf("[%d] send to %d (%lf)\n", i, i - (1 << (k - 1)), (*R)[0]);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        k++;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    append(*Y_heads, Y_heads_size, Y_size);

    // if (id == -1)
    // {
    //     for (int i = 0; i < n_part; i++)
    //     {
    //         for (int j = 0; j < n_part; j++)
    //         {
    //             printf("%lf ", (*R)[i * n_part + j]);
    //         }
    //         printf("\n");
    //     }
    // }
}

void construct_TSQR_Q(int id, int m_part, int n_part, double *Y, size_t *Y_heads, double *tau, double **Q_ret)
{
    int k = 0;
    int p = proc_num;
    int i = id;
    double *Q = NULL;
    int Q_dim[2];
    Q_dim[0] = Q_dim[1] = 0;
    if (i == 0)
    {
        Q_dim[0] = n_part;
        Q_dim[1] = n_part;

        Q = malloc(Q_dim[0] * Q_dim[1] * sizeof(double));
        LAPACKE_dlaset(LAPACK_COL_MAJOR, 'A', Q_dim[0], Q_dim[1], 0.0, 1.0, Q, Q_dim[0]);
    }
    while ((1 << k) < p)
    {
        k++;
    }
    while (k >= 1)
    {

        if (i % (1 << k) == 0 && i + (1 << (k - 1)) < p)
        {
            size_t m_Y = Y_heads[k + 1] - Y_heads[k];
            assert(m_Y % n_part == 0);
            m_Y /= n_part;
            assert(Q_dim[0] <= m_Y);
            double *Q_tmp = calloc(m_Y * n_part, sizeof(double));
            LAPACKE_dlaset(LAPACK_COL_MAJOR, 'A', m_Y, n_part, 0.0, 0.0, Q_tmp, m_Y);
            LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', Q_dim[0], Q_dim[1], Q, Q_dim[0], Q_tmp, m_Y);
            LAPACKE_dormqr(LAPACK_COL_MAJOR, 'L', 'N', m_Y, n_part, n_part, &Y[Y_heads[k]], m_Y, &tau[n_part * (k)], Q_tmp, m_Y);
            LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', Q_dim[0], Q_dim[1], Q_tmp, m_Y, Q, Q_dim[0]);
            int Q_send_dim[2];
            Q_send_dim[0] = m_Y - Q_dim[0];
            Q_send_dim[1] = Q_dim[1];
            int ret = MPI_Send(Q_send_dim, 2, MPI_INT, i + (1 << (k - 1)), 0, MPI_COMM_WORLD);
            assert(ret == MPI_SUCCESS);
            double *Q_send = calloc(Q_send_dim[0] * Q_send_dim[1], sizeof(double));
            LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', Q_send_dim[0], Q_send_dim[1], &Q_tmp[Q_dim[0]], m_Y, Q_send, Q_send_dim[0]);

            ret = MPI_Send(Q_send, Q_send_dim[0] * Q_send_dim[1], MPI_DOUBLE, i + (1 << (k - 1)), 0, MPI_COMM_WORLD);

            assert(ret == MPI_SUCCESS);
            free(Q_tmp);
            free(Q_send);
        }
        else if (i % (1 << k) == (1 << (k - 1)))
        {
            MPI_Status st;
            int ret = MPI_Recv(Q_dim, 2, MPI_INT, i - (1 << (k - 1)), 0, MPI_COMM_WORLD, &st);
            assert(ret == MPI_SUCCESS);

            Q = realloc(Q, Q_dim[0] * Q_dim[1] * sizeof(double));

            ret = MPI_Recv(Q, Q_dim[0] * Q_dim[1], MPI_DOUBLE, i - (1 << (k - 1)), 0, MPI_COMM_WORLD, &st);
            assert(ret == MPI_SUCCESS);
        }
        k--;
    }
    size_t m_Y = Y_heads[1] - Y_heads[0];
    assert(m_Y % n_part == 0);
    m_Y /= n_part;

    double *Q_tmp = calloc(m_Y * n_part, sizeof(double));
    LAPACKE_dlaset(LAPACK_COL_MAJOR, 'A', Q_dim[0], Q_dim[1], 0.0, 0.0, Q_tmp, m_Y);
    // printf("------------------------------------\n");
    // printf("[\n");
    // for (int row = 0; row < Q_dim[0]; row++)
    // {
    //     printf("[");
    //     for (int col = 0; col < Q_dim[1]; col++)
    //     {
    //         printf("%lf, ", Q[col * Q_dim[0] + row]);
    //     }
    //     printf("],\n");
    // }
    // printf("]");

    LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', Q_dim[0], Q_dim[1], Q, Q_dim[0], Q_tmp, m_Y);
    // printf("------------------------------------\n");
    // printf("[\n");
    // for (int row = 0; row < m_Y; row++)
    // {
    //     printf("[");
    //     for (int col = 0; col < n_part; col++)
    //     {
    //         printf("%lf, ", Q_tmp[col * m_Y + row]);
    //     }
    //     printf("],\n");
    // }
    // printf("]");
    LAPACKE_dormqr(LAPACK_COL_MAJOR, 'L', 'N', m_Y, n_part, n_part, &Y[Y_heads[0]], m_Y, &tau[0], Q_tmp, m_Y);
    MPI_Barrier(MPI_COMM_WORLD);

    free(Q);
    *Q_ret = Q_tmp;
    return;
    // LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', Q_dim[0], Q_dim[1], Q_tmp, m_Y, Q, Q_dim[0]);
}

void modified_LU_decomposition(int id, int m_part, int n_part, double *Y, double **S_ret)
{
    // m is a number of row of ***local*** matrix!
    double *S = calloc(n_part, sizeof(double));
    if (id == 0)
    {
        for (int i = 0; i < n_part; i++)
        {
            if (Y[i + i * m_part] < 0)
            {
                S[i] = 1.0;
            }
            else
            {
                S[i] = -1.0;
            }
            Y[i + i * m_part] -= S[i];
            cblas_dscal(m_part - 1 - i, 1.0 / Y[i + i * m_part], &Y[i + 1 + i * m_part], 1);
            cblas_dger(CblasColMajor, m_part - 1 - i, n_part - 1 - i, -1.0, &Y[i + 1 + i * m_part], 1, &Y[i + (i + 1) * m_part], m_part, &Y[i + 1 + (i + 1) * m_part], m_part);
        }
        *S_ret = S;
        double *U = malloc(n_part * n_part * sizeof(double));
        LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'U', n_part, n_part, Y, m_part, U, n_part);
        MPI_Bcast(U, n_part * n_part, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(S, n_part, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        free(U);
    }
    else
    {
        double *U = malloc(n_part * n_part * sizeof(double));
        MPI_Bcast(U, n_part * n_part, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(S, n_part, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // printf("[%d] U %lf\n", rank, U[0]);
        // printf("[%d] S %lf\n", rank, S[0]);
        // printf("[%d] Y %lf\n", rank, Y[0]);
        // printf("[%d] m %d\n", rank, m);
        // printf("[%d] n %d\n", rank, n);
        *S_ret = S;

        cblas_dtrsm(CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m_part, n_part, 1.0, U, n_part, Y, m_part);
        free(U);
    }
    // free(S);
}

void TSQR_HR(int rank_2d, int proc_row_id, int proc_col_id, int m, int n, Matrix *A, int row, int col, Matrix *T_ret)
{
    int id;
    double *data;
    int rank;
    int nrow, ncol, tmp;
    blacs_gridinfo_(&icontext_1d, &nrow, &ncol, &rank, &tmp);

    TSQR_init(rank, m, n, A, row, col, &data);
    return;
    id = rank;
    MPI_Barrier(MPI_COMM_WORLD);
    double *Y, *R, *tau;
    size_t *Y_heads;
    int m_part = m / proc_num;
    if (rank == proc_num - 1)
    {
        m_part += m % proc_num;
    }

    int n_part = n;
    TSQR(id, m_part, n_part, &data, &Y, &Y_heads, &R, &tau);

    // assert(m % proc_num == 0);
    MPI_Barrier(MPI_COMM_WORLD);
    // if (id == 0)
    // {
    //     printf("R=[\n");
    //     for (int row = 0; row < n_part; row++)
    //     {
    //         printf("[");
    //         for (int col = 0; col < n_part; col++)
    //         {
    //             printf("%.18lf, ", R[col * n_part + row]);
    //             fflush(stdout);
    //         }
    //         printf("],\n");
    //     }
    //     printf("]\n");
    // }

    double *Q;
    construct_TSQR_Q(id, m_part, n_part, Y, Y_heads, tau, &Q);
    // if (id == 0)
    //     printf("Q=[\n");
    // MPI_Barrier(MPI_COMM_WORLD);
    // for (int p = 0; p < proc_num; p++)
    // {
    //     if (p == id)
    //     {
    //         // printf("#p = %d\n", id);
    //         for (int i = 0; i < m / proc_num; i++)
    //         {
    //             printf("[");
    //             fflush(stdout);
    //             for (int j = 0; j < n; j++)
    //             {
    //                 printf("%.18lf ,", Q[j * m / proc_num + i]);
    //                 fflush(stdout);
    //             }
    //             printf("],\n");
    //             fflush(stdout);
    //         }
    //     }
    //     else
    //     {
    //         // MPI_Barrierはfflushではprintfの出力順を指定できないので，sleep(1)で無理やり出力を安定させる
    //         sleep(1);
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     // blacs_barrier_(&icontext, ADDR(char, 'A'));
    // }
    // if (id == 0)
    //     printf("]\n");
    MPI_Barrier(MPI_COMM_WORLD);

    double *S; // = malloc(n * sizeof(double));
    modified_LU_decomposition(id, m_part, n, Q, &S);

    double *T;
    if (id == 0)
    {

        double *SY1 = calloc(n * n, sizeof(double));
        T = calloc(n * n, sizeof(double));
        LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'L', n, n, Q, m_part, SY1, n);
        LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'U', n, n, Q, m_part, T, n);

        LAPACKE_dtrtri(LAPACK_COL_MAJOR, 'L', 'U', n, SY1, n);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                T[i + j * n] *= S[j];
                R[i + j * n_part] *= S[i];
            }
        }

        cblas_dtrmm(CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, n, n, -1.0, SY1, n, T, n);
        free(SY1);
        LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'U', n, n, R, n, Q, m_part);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // if (id == 0)
    // {
    //     printf("T=[\n");
    //     for (int row = 0; row < n_part; row++)
    //     {
    //         printf("[");
    //         for (int col = 0; col < n_part; col++)
    //         {
    //             printf("%.18lf, ", T[col * n_part + row]);
    //             fflush(stdout);
    //         }
    //         printf("],\n");
    //     }
    //     printf("]\n");
    // printf("R=[\n");
    // for (int row = 0; row < n_part; row++)
    // {
    //     printf("[");
    //     for (int col = 0; col < n_part; col++)
    //     {
    //         printf("%.18lf, ", R[col * n_part + row]);
    //         fflush(stdout);
    //     }
    //     printf("],\n");
    // }
    // for (int row = n_part; row < m; row++)
    // {
    //     printf("[");
    //     for (int col = 0; col < n_part; col++)
    //     {
    //         printf("0.0, ");
    //         fflush(stdout);
    //     }
    //     printf("],\n");
    // }

    //     printf("]\n");
    // }
    int pad = 0;
    if (m % proc_num != 0)
    {

        if (id == proc_num - 1)
        {
            pad = m % proc_num;

            double *send = malloc(pad * n_part * sizeof(double));
            double *shirnk = malloc((m_part - pad) * n_part * sizeof(double));
            LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', m_part - pad, n_part, Q, m_part, shirnk, m_part - pad);
            LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', pad, n_part, &(Q[m_part - pad]), m_part, send, pad);
            int ret = MPI_Send(send, pad * n_part, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            assert(ret == MPI_SUCCESS);
            free(send);
            free(Q);
            Q = shirnk;
            m_part -= pad;
            pad = 0;
        }

        if (id == 0)
        {
            pad = m % proc_num;

            double *recv = malloc(pad * n_part * sizeof(double));
            double *expand = malloc((m_part + pad) * n_part * sizeof(double));
            MPI_Status st;
            int ret = MPI_Recv(recv, pad * n_part, MPI_DOUBLE, proc_num - 1, 0, MPI_COMM_WORLD, &st);
            assert(ret == MPI_SUCCESS);
            LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', m_part, n_part, Q, m_part, expand, (m_part + pad));
            LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', pad, n_part, recv, pad, &(expand[m_part]), (m_part + pad));
            free(Q);
            free(recv);
            Q = expand;
            //  m_part += pad;
        }
    }

    // if (id == 0)
    // {
    //     printf("S=[\n");
    //     for (int i = 0; i < n; i++)
    //     {
    //         printf("[");
    //         for (int j = 0; j < n; j++)
    //         {
    //             if (i == j)
    //             {
    //                 printf("%.18lf, ", S[i]);
    //             }
    //             else
    //             {
    //                 printf("0, ");
    //             }
    //         }
    //         printf("],\n");
    //     }
    //     printf("]\n");

    //     printf("U=[\n");
    //     for (int i = 0; i < n; i++)
    //     {
    //         printf("[");
    //         for (int j = 0; j < n; j++)
    //         {
    //             if (i <= j)
    //             {
    //                 printf("%.18lf, ", Q[i + j * m_part]);
    //             }
    //             else
    //             {
    //                 printf("0, ");
    //             }
    //         }
    //         printf("],\n");
    //     }
    //     printf("]\n");

    //     printf("Y=[\n");
    //     for (int i = 0; i < m_part; i++)
    //     {
    //         printf("[");
    //         for (int j = 0; j < n; j++)
    //         {
    //             if (i > j)
    //             {
    //                 printf("%.18lf, ", Q[i + j * m_part]);
    //             }
    //             else if (i == j)
    //             {
    //                 printf("1, ");
    //             }
    //             else
    //             {
    //                 printf("0, ");
    //             }
    //         }
    //         printf("],\n");
    //     }
    // }
    // MPI_Barrier(MPI_COMM_WORLD);
    // sleep(1);

    // for (int p = 1; p < proc_num; p++)
    // {
    //     if (p == id)
    //     {
    //         // printf("#p = %d\n", id);
    //         for (int i = 0; i < m_part; i++)
    //         {
    //             printf("[");
    //             fflush(stdout);
    //             for (int j = 0; j < n; j++)
    //             {
    //                 printf("%.18lf ,", Q[j * m_part + i]);
    //                 fflush(stdout);
    //             }
    //             printf("],\n");
    //             fflush(stdout);
    //         }
    //     }
    //     else
    //     {
    //         // MPI_Barrierはfflushではprintfの出力順を指定できないので，sleep(1)で無理やり出力を安定させる
    //         sleep(1);
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     // blacs_barrier_(&icontext, ADDR(char, 'A'));
    // }
    // if (id == 0)
    //     printf("]\n");
    // MPI_Barrier(MPI_COMM_WORLD);

    int desc[DESC_LEN];
    int ierror;
    // printf("[%d] m_part %d\n", id, m_part);
    // printf("[%d] pad %d\n", id, pad);
    printf("[%d] descinit 1 %d %d %d %d %d %d\n", rank, m, n, m_part, n, icontext_1d, m_part + pad);
    descinit_(&desc, ADDR(int, m),
              ADDR(int, n), ADDR(int, m_part), ADDR(int, n), ADDR(int, 0), ADDR(int, 0), &icontext_1d, ADDR(int, m_part + pad), &ierror);
    printf("[%d] done 1 %d %d %d %d %d %d\n", rank, m, n, m_part, n, icontext_1d, m_part + pad);
    Matrix mat;
    mat.data = Q;
    mat.desc = desc;
    pdgemr2d_wrap(m, n, &mat, 0, 0, A, row, col);

    MPI_Barrier(MPI_COMM_WORLD);

    printf("[%d] descinit 2 %d %d %d %d %d %d\n", rank, n, n, n, n, icontext_1d, n);
    descinit_(&desc, ADDR(int, n),
              ADDR(int, n), ADDR(int, n), ADDR(int, n), ADDR(int, 0), ADDR(int, 0), &icontext_1d, ADDR(int, n), &ierror);
    printf("[%d] done 2 %d %d %d %d %d %d\n", rank, n, n, n, n, icontext_1d, n);
    mat.data = T;
    mat.desc = desc;
    pdgemr2d_wrap(n, n, &mat, 0, 0, T_ret, 0, 0);
    MPI_Barrier(MPI_COMM_WORLD);

    return;
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
    int rank;
    initBuffer();
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

    blacs_get_(ADDR(int, 0), ADDR(int, 0), &icontext_2d);
    blacs_gridinit_(&icontext_2d, ADDR(char, 'R'), &proc_row_num, &proc_col_num);
    blacs_gridinfo_(&icontext_2d, &proc_row_num, &proc_col_num, &my_row, &my_col);

    blacs_get_(ADDR(int, 0), ADDR(int, 0), &icontext_1d);
    blacs_gridinit_(&icontext_1d, ADDR(char, 'R'), ADDR(int, proc_row_num *proc_col_num), ADDR(int, 1));

    // main calc
    // printf("%d %d", m, n);
    if (rank == 0 && print_checkcode)
    {
        printf("import numpy as np\n");
    }
    n = 2;
    Matrix *A = create_matrix(proc_row_num, proc_col_num, m, n, block_row, block_col);
    Matrix *T = create_matrix(proc_row_num, proc_col_num, L, m, block_row, block_col);
    Matrix *Y = create_matrix(proc_row_num, proc_col_num, m, m, block_row, block_col);
    measure_time(for (size_t i = 0; i < A->global_row; ++i) {
        for (size_t j = i; j < A->global_col; ++j)
        {
            double r = (double)(rand()) / RAND_MAX;
            set(A, i, j, r);
            set(A, j, i, r);
        }
    });
    blacs_barrier_(&icontext_2d, ADDR(char, 'A'));
    if (print_checkcode)
    {
        rprintf("import numpy as np\n");

        print_matrix("A=", A, rank);
    }

    // measure_time(bischof(rank, proc_row_num, proc_col_num, m, L, A, T, Y));
    // if (print_checkcode)
    // {
    //     for (int i = 0; i < m; i++)
    //     {
    //         for (int j = i; j < m; j++)
    //         {
    //             set(A, i, j, get(A, j, i));
    //         }
    //     }

    //     for (int i = 0; i < m; i++)
    //     {
    //         for (int j = 0; j < m; j++)
    //         {
    //             if (abs(i - j) > L)
    //             {
    //                 set(A, i, j, 0.0);
    //             }
    //         }
    //     }
    //     print_matrix("B=", A, rank);
    //     rprintf("\nA = np.matrix(A)\nB = np.matrix(B)\ne=0\nfor i, j in zip(sorted(np.linalg.eigvals(A)), sorted(np.linalg.eigvals(B))):\n    print(i, j)\n    e+=(i-j)**2\n\nprint('error=',e**0.5)\n");
    //     // print_matrix("T=", T, rank);
    //     // print_matrix("Y=", Y, rank);
    // }

    MPI_Barrier(MPI_COMM_WORLD);
    // print_matrix("A=", A, rank);
    MPI_Barrier(MPI_COMM_WORLD);

    TSQR_HR(rank, my_row, my_col, m, n, A, 0, 0, T);
    MPI_Barrier(MPI_COMM_WORLD);
    // print_matrix("ret=", A, rank);
    // print_matrix("T=", T, rank);

    free_matrix(A);

    MPI_Finalize();
    if (rank == 0)
        dumpBuffer();
    return 0;
}
