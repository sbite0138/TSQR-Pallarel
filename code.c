#include <assert.h>
#ifdef DEBUG
#include <cblas.h>
#include <lapacke.h>
#else
#include <mkl.h>
#endif
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
#define rprintf(...)               \
    do                             \
    {                              \
        if (rank_global == 0)      \
        {                          \
            sbprintf(__VA_ARGS__); \
        }                          \
                                   \
    } while (0)

#define measure_time(x)                                                                                                                            \
    do                                                                                                                                             \
    {                                                                                                                                              \
        if (print_checkcode == false)                                                                                                              \
            rprintf("# enter %d\n", __LINE__);                                                                                                     \
        level++;                                                                                                                                   \
        double start = omp_get_wtime();                                                                                                            \
        {                                                                                                                                          \
            x;                                                                                                                                     \
        }                                                                                                                                          \
        double end = omp_get_wtime();                                                                                                              \
        level--;                                                                                                                                   \
        if (print_checkcode == false)                                                                                                              \
            rprintf("@ {\"rank\":%d,\"level\":%d, \"line\":%d, \"cmd\":\"%s\", \"time\":%.18f}\n", rank_global, level, __LINE__, #x, end - start); \
        if (print_checkcode == false)                                                                                                              \
            rprintf("# exit %d\n", __LINE__);                                                                                                      \
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
        fprintf(stderr, "realloced! %d\n", sBuffer.size);
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
bool use_tsqr = true;
int rank_global;

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

double pdlange_(char *norm, int *m, int *n, double *a, int *ia, int *ja, int *desca, double *work);

double pdlange_wrap(Matrix *mat)
{
    double *work;
    return pdlange_(ADDR(char, 'F'), &(mat->global_row), &(mat->global_col), mat->data, ADDR(int, 1), ADDR(int, 1), mat->desc, work);
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

    measure_time(tau = malloc(((size_t)n + (size_t)col) * sizeof(double)));
    measure_time(pdgeqrf_wrap(m, n, matrix, row, col, tau));
    // double *work = malloc(n * (n - 1) / 2);

    // pdlarft_(ADDR(char, 'F'), ADDR(char, 'C'), ADDR(int, 1),
    //          ADDR(int, n), matrix->data, ADDR(int, row + 1), ADDR(int, col + 1), matrix->desc, tau, T->data, work);
    // free(work);
    // return;
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

    measure_time(Y = create_matrix(proc_row, proc_col, m, n, block_row, block_col));
    measure_time(pdlaset_wrap('A', T->global_row, T->global_col, 0.0, 0.0, T, 0, 0));
    measure_time(pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col), desc));

    // set(T, 0, 0, tau[0]);
    measure_time(set(T, 0, 0, val));
    // print_matrix("T=", T, rank);
    measure_time(pdgemr2d_wrap(m, n, matrix, row, col, Y, 0, 0));
    measure_time(pdlaset_wrap('U', Y->global_row, Y->global_col, 0.0, 1.0, Y, 0, 0));
    Matrix *y;
    Matrix *z;
    Matrix *tmp;

    measure_time(y = create_matrix(proc_row, proc_col, m, n, block_row, block_col));
    measure_time(z = create_matrix(proc_row, proc_col, n, 1, block_row, block_col));
    measure_time(tmp = create_matrix(proc_row, proc_col, n, 1, block_row, block_col));

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
    measure_time(free_matrix(tmp));
    measure_time(free_matrix(z));
    measure_time(free_matrix(y));
    measure_time(free_matrix(Y));
    measure_time(free(tau));
}

bool can_TSQR(int m, int n)
{
    int block_size = m / proc_num;
    // printf("block_size %d, m %d, proc_num %d, m mod proc_num %d\n", block_size, m, proc_num, m % proc_num);
    // rprintf("m %d n %d proc_num %d\n", m, n, proc_num)

    return use_tsqr && (m >= proc_num * n);
}

void bischof(int rank, int nproc_row, int nproc_col, int N, int L, Matrix *A, Matrix *T, Matrix *Y)
{
    int nb = L;
    int ldt_iter = nb;
    assert(MIN(N, L) >= nb && nb >= 1);
    assert(ldt_iter >= nb);
    assert(N % L == 0);

    Matrix *T_iter;
    Matrix *V;
    Matrix *update_tmp;
    Matrix *update_P;
    Matrix *update_beta;
    Matrix *update_Q;
    measure_time(T_iter = create_matrix(nproc_row, nproc_col, L, L, block_col, block_row));
    measure_time(V = create_matrix(nproc_row, nproc_col, N - L, L, block_row, block_col));
    measure_time(update_tmp = create_matrix(nproc_row, nproc_col, N - L, N - L, block_row, block_col));
    measure_time(update_P = create_matrix(nproc_row, nproc_col, N - L, L, block_row, block_col));
    measure_time(update_beta = create_matrix(nproc_row, nproc_col, L, L, block_row, block_col));
    measure_time(update_Q = create_matrix(nproc_row, nproc_col, N - L, L, block_row, block_col));

    for (int k = 0; k < N / L - 1; k++)
    {
        //  rprintf("iteration %d/%d\n", k + 1, N / L - 1);
        int Nk = N - L - k * L;
        // A???(k+1,k)?????????????????????QR????????????
        if (can_TSQR(Nk, L))
        {
            rprintf("# TSQR\n");
            measure_time(TSQR_HR(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter));
        }
        else
        {
            rprintf("# Householder\n");
            measure_time(pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter));
        }
        measure_time(pdlaset_wrap('L', L, L, 0.0, 0.0, T, 0, L * k));
        measure_time(pdgemr2d_wrap(L, L, T_iter, 0, 0, T, 0, L * k));
        // print_matrix("A part = ", Nk, L, a_part, lda);
        int pad_row = (k + 1) * L;
        int pad_col = k * L;

        // V???A_part?????????????????????????????????(LAPACKE_dgeqrt??????????????????????????????a_part????????????????????????Q???compact-WY?????????????????????????????????)
        measure_time(pdgemr2d_wrap(Nk, L, A, pad_row, pad_col, V, 0, 0));
        measure_time(pdlaset_wrap('U', Nk, L, 0.0, 1.0, V, 0, 0));
        // print_matrix("V = ", Nk, L, V, Nk);

        // Y???V?????????
        measure_time(pdgemr2d_wrap(Nk, L, V, 0, 0, Y, 0, L * k));

        // print_matrix("Y=", Y, rank);

        // ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????Bischof???????????????????????????P???????????????
        measure_time(pdsymm_wrap('L', 'L', Nk, L, 1.0, A, (k + 1) * L, (k + 1) * L, V, 0, 0, 0.0, update_tmp, 0, 0));
        // print_matrix("T_iter", T_iter, rank);

        measure_time(pdgemm_wrap('N', 'N', Nk, L, L, 1.0, update_tmp, 0, 0, T_iter, 0, 0, 0.0, update_P, 0, 0));
        // print_matrix("update_P", update_P, rank);

        // beta???????????????
        measure_time(pdgemm_wrap('T', 'N', L, L, Nk, 0.5, V, 0, 0, update_P, 0, 0, 0.0, update_tmp, 0, 0));
        measure_time(pdgemm_wrap('T', 'N', L, L, L, 1.0, T_iter, 0, 0, update_tmp, 0, 0, 0.0, update_beta, 0, 0));
        // print_matrix("update_P", update_beta, rank);

        //  Q???????????????
        measure_time(pdgemr2d_wrap(Nk, L, update_P, 0, 0, update_Q, 0, 0));
        // print_matrix("update_Q=", update_Q, rank);

        measure_time(pdgemm_wrap('N', 'N', Nk, L, L, -1.0, V, 0, 0, update_beta, 0, 0, 1.0, update_Q, 0, 0));
        // print_matrix("update_Q", Nk, L, update_Q, Nk);

        // A???rank-2k????????????
        // print_matrix("V=", V, rank);
        // print_matrix("update_Q=", update_Q, rank);

        measure_time(pdsyr2k_wrap('L', 'N', Nk, L, -1.0, V, 0, 0, update_Q, 0, 0, 1.0, A, (k + 1) * L, (k + 1) * L));
    }
    measure_time(free_matrix(V));
    measure_time(free_matrix(update_P));
    measure_time(free_matrix(update_beta));
    measure_time(free_matrix(update_Q));
    measure_time(free_matrix(update_tmp));
    measure_time(free_matrix(T_iter));
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
    int recv1 = -1;
    int recv2 = -1;
    int send1 = -1;
    int send2 = -1;
    int m_recv1 = 0;
    int m_recv2 = 0;
    double *data_recv1 = NULL;
    double *data_recv2 = NULL;
    double *data_send1 = NULL;
    double *data_send2 = NULL;
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
    // printf("status [%d] %d %d %d %d\n", rank, send1, send2, recv1, recv2);
    // printf("have/need [%d] %d %d %d %d\n", rank, have1, have2, need1, need2);
    // printf("send   [%d] %d %d\n", rank, m_head, m_tail);

    MPI_Request reqs[4];

    int reqs_idx = 0;
    // send 1
    measure_time(data_send1 = malloc(m_head * n * sizeof(double)));
    measure_time(LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', m_head, n, mat->data, m_total, data_send1, m_head));

    measure_time(MPI_Isend(data_send1, m_head * n, MPI_DOUBLE, send1, have1, MPI_COMM_WORLD, &reqs[reqs_idx++]));
    //  recv 1
    m_recv1 = (need1 == block_num - 1) ? m_last : m_part;
    measure_time(data_recv1 = malloc(m_recv1 * n * sizeof(double)));
    measure_time(MPI_Irecv(data_recv1, m_recv1 * n, MPI_DOUBLE, recv1, need1, MPI_COMM_WORLD, &reqs[reqs_idx++]));

    // send 2
    if (send2 != -1)
    {
        measure_time(data_send2 = malloc(m_tail * n * sizeof(double)));
        measure_time(LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', m_tail, n, &(mat->data[m_head]), m_total, data_send2, m_tail));
        measure_time(MPI_Isend(data_send2, m_tail * n, MPI_DOUBLE, send2, have2, MPI_COMM_WORLD, &reqs[reqs_idx++]));
    }

    // recv 2
    if (recv2 != -1)
    {
        m_recv2 = (need2 == block_num - 1) ? m_last : m_part;
        measure_time(data_recv2 = malloc(m_recv2 * n * sizeof(double)));
        measure_time(MPI_Irecv(data_recv2, m_recv2 * n, MPI_DOUBLE, recv2, need2, MPI_COMM_WORLD, &reqs[reqs_idx++]));
        // printf("recv   [%d] %d \n", rank, m_recv2);
        // printf("mlast   [%d] %d\n", rank, (need2 == block_num - 1) ? m_last : m_part);
    }

    for (int i = 0; i < reqs_idx; i++)
    {
        MPI_Status st;
        measure_time(MPI_Wait(&reqs[i], &st));
    }

    if (recv2 == -1)
    {
        // printf("[%d] %p\n", rank, mat->data);
        measure_time(mat->data = realloc(mat->data, (m_recv1)*n * sizeof(double)));
        // printf("[%d] %p\n", rank, mat->data);
        measure_time(LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', m_recv1, n, data_recv1, m_recv1, mat->data, m_recv1));
        mat->global_row = mat->local_row = mat->leading_dimension = m_recv1;
    }
    else
    {
        measure_time(mat->data = realloc(mat->data, (m_recv1 + m_recv2) * n * sizeof(double)));

        measure_time(LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', m_recv1, n, data_recv1, m_recv1, mat->data, m_recv1 + m_recv2));
        measure_time(LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', m_recv2, n, data_recv2, m_recv2, &((mat->data)[m_recv1]), m_recv1 + m_recv2));
        mat->global_row = mat->local_row = mat->leading_dimension = m_recv1 + m_recv2;
    }
    if (data_recv1 != NULL)
    {
        measure_time(free(data_recv1));
    }
    if (data_recv2 != NULL)
    {
        measure_time(free(data_recv2));
    }
    if (data_send1 != NULL)
    {
        measure_time(free(data_send1));
    }
    if (data_send2 != NULL)
    {
        measure_time(free(data_send2));
    }
}

// TODO: exchange??????1??????????????????
void exchange2(int m_global, int n, int m_local, Matrix *mat, int rank)
{
    size_t m_part = m_global / proc_num;
    size_t block_num = (m_global + m_part - 1) / m_part;
    int m_head = m_global / proc_num;
    assert(m_head <= m_local);
    int m_tail = m_local - m_head;
    int m_last;
    if (m_global % m_part == 0)
    {
        m_last = m_part;
    }
    else
    {
        m_last = m_global % m_part;
    }
    int have1 = -1;
    int have2 = -1;
    int need1 = -1;
    int need2 = -1;

    int t = proc_num - (block_num % proc_num);

    if (rank < t)
    {
        have1 = rank;
    }
    else
    {
        have1 = 2 * (rank - t) + t;
        have2 = have1 + 1;
    }

    if (rank < block_num % proc_num)
    {
        need1 = rank;
        need2 = proc_num + rank;
    }
    else
    {
        need1 = rank;
    }

    int recv1 = -1;
    int recv2 = -1;
    int send1 = -1;
    int send2 = -1;
    int m_recv1 = 0;
    int m_recv2 = 0;
    double *data_recv1 = NULL;
    double *data_recv2 = NULL;
    double *data_send1 = NULL;
    double *data_send2 = NULL;

    // set recv
    if (need1 < t)
    {
        recv1 = need1;
    }
    else
    {
        recv1 = (need1 - t) / 2 + t;
    }
    if (need2 != -1)
    {
        if (need2 < t)
        {
            recv2 = need2;
        }
        else
        {
            recv2 = (need2 - t) / 2 + t;
        }
    }

    send1 = have1 % proc_num;
    if (have2 != -1)
    {
        send2 = have2 % proc_num;
    }
    // printf("status [%d] %d %d %d %d\n", rank, send1, send2, recv1, recv2);
    // printf("have/need [%d] %d %d %d %d\n", rank, have1, have2, need1, need2);
    // printf("send   [%d] %d %d\n", rank, m_head, m_tail);

    MPI_Request reqs[4];

    int reqs_idx = 0;
    // send 1
    measure_time(data_send1 = malloc(m_head * n * sizeof(double)));
    measure_time(LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', m_head, n, mat->data, m_local, data_send1, m_head));

    measure_time(MPI_Isend(data_send1, m_head * n, MPI_DOUBLE, send1, have1, MPI_COMM_WORLD, &reqs[reqs_idx++]));
    //  recv 1
    m_recv1 = (need1 == block_num - 1) ? m_last : m_part;
    measure_time(data_recv1 = malloc(m_recv1 * n * sizeof(double)));
    measure_time(MPI_Irecv(data_recv1, m_recv1 * n, MPI_DOUBLE, recv1, need1, MPI_COMM_WORLD, &reqs[reqs_idx++]));

    // send 2
    if (send2 != -1)
    {
        measure_time(data_send2 = malloc(m_tail * n * sizeof(double)));
        measure_time(LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', m_tail, n, &(mat->data[m_head]), m_local, data_send2, m_tail));
        measure_time(MPI_Isend(data_send2, m_tail * n, MPI_DOUBLE, send2, have2, MPI_COMM_WORLD, &reqs[reqs_idx++]));
    }

    // recv 2
    if (recv2 != -1)
    {
        m_recv2 = (need2 == block_num - 1) ? m_last : m_part;
        measure_time(data_recv2 = malloc(m_recv2 * n * sizeof(double)));
        measure_time(MPI_Irecv(data_recv2, m_recv2 * n, MPI_DOUBLE, recv2, need2, MPI_COMM_WORLD, &reqs[reqs_idx++]));
        // printf("recv   [%d] %d \n", rank, m_recv2);
        // printf("mlast   [%d] %d\n", rank, (need2 == block_num - 1) ? m_last : m_part);
    }

    for (int i = 0; i < reqs_idx; i++)
    {
        MPI_Status st;
        measure_time(MPI_Wait(&reqs[i], &st));
    }

    if (recv2 == -1)
    {
        // printf("[%d] %p\n", rank, mat->data);
        measure_time(mat->data = realloc(mat->data, (m_recv1)*n * sizeof(double)));
        // printf("[%d] %p\n", rank, mat->data);
        measure_time(LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', m_recv1, n, data_recv1, m_recv1, mat->data, m_recv1));
        mat->global_row = mat->local_row = mat->leading_dimension = m_recv1;
    }
    else
    {
        measure_time(mat->data = realloc(mat->data, (m_recv1 + m_recv2) * n * sizeof(double)));

        measure_time(LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', m_recv1, n, data_recv1, m_recv1, mat->data, m_recv1 + m_recv2));
        measure_time(LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', m_recv2, n, data_recv2, m_recv2, &((mat->data)[m_recv1]), m_recv1 + m_recv2));
        mat->global_row = mat->local_row = mat->leading_dimension = m_recv1 + m_recv2;
    }
    if (data_recv1 != NULL)
    {
        measure_time(free(data_recv1));
    }
    if (data_recv2 != NULL)
    {
        measure_time(free(data_recv2));
    }
    if (data_send1 != NULL)
    {
        measure_time(free(data_send1));
    }
    if (data_send2 != NULL)
    {
        measure_time(free(data_send2));
    }
}

void TSQR_init(int rank, int m, int n, Matrix *matrix, int row, int col, double **data, int *m_ret)
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
    int m_total;
    measure_time(m_total = numroc_(&m, &m_part, &rank, ADDR(int, 0), &proc_num));
    // printf("[%d] descinit 3 %d %d %d %d %d %d\n", rank, m, n, m_part, n_part, icontext_1d, m_total);
    measure_time(descinit_(desc, ADDR(int, m), ADDR(int, n), &m_part, &n, ADDR(int, 0), ADDR(int, 0), &icontext_1d, &m_total, &ierror));
    // printf("[%d] done 3 %d %d %d %d %d %d\n", rank, m, n, m_part, n_part, icontext_1d, m_total);
    Matrix mat;
    mat.desc = desc;
    measure_time(mat.data = malloc(((size_t)m_total) * (size_t)n_part * sizeof(double)));
    // printf("[%d] %d %d %d %d\n", rank, m, n, row, col);
    measure_time(pdgemr2d_wrap(m, n, matrix, row, col, &mat, 0, 0));
    int m_last;
    if (m % m_part == 0)
    {
        m_last = m_part;
    }
    else
    {
        m_last = m % m_part;
    }
    measure_time(exchange(m_total, n, m_part, m_last, block_num, have1, have2, need1, need2, &mat, rank));
    measure_time(MPI_Barrier(MPI_COMM_WORLD));
    *data = mat.data;
    *m_ret = mat.global_row;
}

// TODO:??????????????????????????????????????????
void TSQR(int id, int m_part, int n_part, double **data, double **Y, size_t **Y_heads, double **R, double **tau)
{
    measure_time(MPI_Barrier(MPI_COMM_WORLD));
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
    measure_time(*tau = malloc((size_t)n_part * sizeof(double)));
    size_t tau_size = n_part;
    measure_time(*R = calloc(n_part * n_part, sizeof(double)));
    measure_time(*Y = calloc(m_part * n_part, sizeof(double)));
    size_t Y_heads_size = 0;
    *Y_heads = NULL;
    measure_time(append(*Y_heads, Y_heads_size, 0));
    size_t Y_size = m_part * n_part;

    int current_blocknum = 1;
    measure_time(LAPACKE_dgeqrf(LAPACK_COL_MAJOR, m_part, n_part, *data, m_part, *tau));
    measure_time(LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'U', n_part, n_part, *data, m_part, *R, n_part));
    measure_time(LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'L', m_part, n_part, *data, m_part, *Y, m_part));
    measure_time(LAPACKE_dlaset(LAPACK_COL_MAJOR, 'U', m_part, n_part, 0.0, 1.0, *Y, m_part));
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
            double *R_tmp;
            measure_time(R_tmp = calloc(2 * n_part * n_part, sizeof(double)));
            double *R_rsv;
            measure_time(R_rsv = calloc(n_part * n_part, sizeof(double)));
            MPI_Status st;
            int ret;
            measure_time(ret = MPI_Recv(R_rsv, n_part * n_part, MPI_DOUBLE, j, 0, MPI_COMM_WORLD, &st));
            assert(ret == MPI_SUCCESS);
            // printf("[%d] recv from %d(%lf)\n", i, j, R_rsv[0]);
            //  write to Y

            measure_time(*Y = realloc(*Y, (Y_size + n_part * n_part + n_part * n_part) * sizeof(double)));
            measure_time(LAPACKE_dlaset(LAPACK_COL_MAJOR, 'A', 2 * n_part, n_part, 0.0, 0.0, &((*Y)[Y_size]), 2 * n_part));

            measure_time(append(*Y_heads, Y_heads_size, Y_size));

            measure_time(*tau = realloc(*tau, (tau_size + n_part) * sizeof(double)));
            // write to R
            measure_time(LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'U', n_part, n_part, *R, n_part, &((*Y)[Y_size]), 2 * n_part));
            measure_time(LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'U', n_part, n_part, R_rsv, n_part, &((*Y)[Y_size + n_part]), 2 * n_part));

            measure_time(ret = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, 2 * n_part, n_part, &((*Y)[Y_size]), 2 * n_part, &((*tau)[tau_size])));
            // printf("[%d] ret = %d\n", i, ret);
            measure_time(LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'U', n_part, n_part, &((*Y)[Y_size]), 2 * n_part, *R, n_part));
            // printf("[%d] R %lf\n", i, (*R)[0]);
            Y_size += n_part * n_part * 2;
            tau_size += n_part;
            measure_time(free(R_tmp));
            measure_time(free(R_rsv));
        }
        else if (i % (1 << k) == (1 << (k - 1)))
        {
            int ret;
            measure_time(ret = MPI_Send(*R, n_part * n_part, MPI_DOUBLE, i - (1 << (k - 1)), 0, MPI_COMM_WORLD));
            assert(ret == MPI_SUCCESS);
            // printf("[%d] send to %d (%lf)\n", i, i - (1 << (k - 1)), (*R)[0]);
        }
        measure_time(MPI_Barrier(MPI_COMM_WORLD));

        k++;
    }
    measure_time(MPI_Barrier(MPI_COMM_WORLD));

    measure_time(append(*Y_heads, Y_heads_size, Y_size));

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

        measure_time(Q = malloc(Q_dim[0] * Q_dim[1] * sizeof(double)));
        measure_time(LAPACKE_dlaset(LAPACK_COL_MAJOR, 'A', Q_dim[0], Q_dim[1], 0.0, 1.0, Q, Q_dim[0]));
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
            measure_time(LAPACKE_dlaset(LAPACK_COL_MAJOR, 'A', m_Y, n_part, 0.0, 0.0, Q_tmp, m_Y));
            measure_time(LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', Q_dim[0], Q_dim[1], Q, Q_dim[0], Q_tmp, m_Y));
            measure_time(LAPACKE_dormqr(LAPACK_COL_MAJOR, 'L', 'N', m_Y, n_part, n_part, &Y[Y_heads[k]], m_Y, &tau[n_part * (k)], Q_tmp, m_Y));
            measure_time(LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', Q_dim[0], Q_dim[1], Q_tmp, m_Y, Q, Q_dim[0]));
            int Q_send_dim[2];
            Q_send_dim[0] = m_Y - Q_dim[0];
            Q_send_dim[1] = Q_dim[1];
            int ret;
            measure_time(ret = MPI_Send(Q_send_dim, 2, MPI_INT, i + (1 << (k - 1)), 0, MPI_COMM_WORLD));
            assert(ret == MPI_SUCCESS);
            double *Q_send;
            measure_time(Q_send = calloc(Q_send_dim[0] * Q_send_dim[1], sizeof(double)));
            measure_time(LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', Q_send_dim[0], Q_send_dim[1], &Q_tmp[Q_dim[0]], m_Y, Q_send, Q_send_dim[0]));

            measure_time(ret = MPI_Send(Q_send, Q_send_dim[0] * Q_send_dim[1], MPI_DOUBLE, i + (1 << (k - 1)), 0, MPI_COMM_WORLD));

            assert(ret == MPI_SUCCESS);
            measure_time(free(Q_tmp));
            measure_time(free(Q_send));
        }
        else if (i % (1 << k) == (1 << (k - 1)))
        {
            MPI_Status st;
            int ret;
            measure_time(ret = MPI_Recv(Q_dim, 2, MPI_INT, i - (1 << (k - 1)), 0, MPI_COMM_WORLD, &st));
            assert(ret == MPI_SUCCESS);

            measure_time(Q = realloc(Q, Q_dim[0] * Q_dim[1] * sizeof(double)));
            measure_time(ret = MPI_Recv(Q, Q_dim[0] * Q_dim[1], MPI_DOUBLE, i - (1 << (k - 1)), 0, MPI_COMM_WORLD, &st));
            assert(ret == MPI_SUCCESS);
        }
        k--;
    }
    size_t m_Y = Y_heads[1] - Y_heads[0];
    assert(m_Y % n_part == 0);
    m_Y /= n_part;

    double *Q_tmp;
    measure_time(Q_tmp = calloc(m_Y * n_part, sizeof(double)));
    measure_time(LAPACKE_dlaset(LAPACK_COL_MAJOR, 'A', Q_dim[0], Q_dim[1], 0.0, 0.0, Q_tmp, m_Y));
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

    measure_time(LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', Q_dim[0], Q_dim[1], Q, Q_dim[0], Q_tmp, m_Y));
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
    measure_time(LAPACKE_dormqr(LAPACK_COL_MAJOR, 'L', 'N', m_Y, n_part, n_part, &Y[Y_heads[0]], m_Y, &tau[0], Q_tmp, m_Y));
    measure_time(MPI_Barrier(MPI_COMM_WORLD));

    measure_time(free(Q));
    *Q_ret = Q_tmp;
    return;
    // LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', Q_dim[0], Q_dim[1], Q_tmp, m_Y, Q, Q_dim[0]);
}

void modified_LU_decomposition(int id, int m_part, int n_part, double *Y, double **S_ret)
{
    // m is a number of row of ***local*** matrix!
    double *S;
    measure_time(S = calloc(n_part, sizeof(double)));
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
            measure_time(cblas_dscal(m_part - 1 - i, 1.0 / Y[i + i * m_part], &Y[i + 1 + i * m_part], 1));
            measure_time(cblas_dger(CblasColMajor, m_part - 1 - i, n_part - 1 - i, -1.0, &Y[i + 1 + i * m_part], 1, &Y[i + (i + 1) * m_part], m_part, &Y[i + 1 + (i + 1) * m_part], m_part));
        }
        *S_ret = S;
        double *U;
        measure_time(U = malloc(n_part * n_part * sizeof(double)));
        measure_time(LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'U', n_part, n_part, Y, m_part, U, n_part));
        measure_time(MPI_Bcast(U, n_part * n_part, MPI_DOUBLE, 0, MPI_COMM_WORLD));
        measure_time(MPI_Bcast(S, n_part, MPI_DOUBLE, 0, MPI_COMM_WORLD));
        measure_time(free(U));
    }
    else
    {
        double *U;
        measure_time(U = malloc(n_part * n_part * sizeof(double)));
        measure_time(MPI_Bcast(U, n_part * n_part, MPI_DOUBLE, 0, MPI_COMM_WORLD));
        measure_time(MPI_Bcast(S, n_part, MPI_DOUBLE, 0, MPI_COMM_WORLD));
        // printf("[%d] U %lf\n", rank, U[0]);
        // printf("[%d] S %lf\n", rank, S[0]);
        // printf("[%d] Y %lf\n", rank, Y[0]);
        // printf("[%d] m %d\n", rank, m);
        // printf("[%d] n %d\n", rank, n);
        *S_ret = S;

        measure_time(cblas_dtrsm(CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m_part, n_part, 1.0, U, n_part, Y, m_part));
        measure_time(free(U));
    }
    // free(S);
}

void TSQR_HR(int rank_2d, int proc_row_id, int proc_col_id, int m, int n, Matrix *A, int row, int col, Matrix *T_ret)
{
    //  print_matrix("A=", A, rank_2d);
    int id;
    double *data;
    int rank;
    int nrow, ncol, tmp;
    measure_time(blacs_gridinfo_(&icontext_1d, &nrow, &ncol, &rank, &tmp));
    int m_local;
    measure_time(TSQR_init(rank, m, n, A, row, col, &data, &m_local));

    measure_time(MPI_Barrier(MPI_COMM_WORLD));
    // for (int p = 0; p < proc_num; p++)
    // {
    //     if (p == rank)
    //     {
    //         // printf("#p = %d\n", id);
    //         for (int i = 0; i < m_part; i++)
    //         {
    //             printf("[");
    //             fflush(stdout);
    //             for (int j = 0; j < n; j++)
    //             {
    //                 printf("%lf ,", data[j * m_part + i]);
    //                 fflush(stdout);
    //             }
    //             printf("],\n");
    //             fflush(stdout);
    //         }
    //     }
    //     else
    //     {
    //         // MPI_Barrier???fflush??????printf??????????????????????????????????????????sleep???????????????????????????????????????
    //         usleep(1000 * 100);
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     // blacs_barrier_(&icontext, ADDR(char, 'A'));
    // }
    id = rank;
    measure_time(MPI_Barrier(MPI_COMM_WORLD));
    double *Y, *R, *tau;
    size_t *Y_heads;
    // int m_part = m / proc_num;
    // if (rank == proc_num - 1)
    // {
    //     m_part += m % proc_num;
    //    }

    int n_part = n;
    measure_time(TSQR(id, m_local, n_part, &data, &Y, &Y_heads, &R, &tau));

    // assert(m % proc_num == 0);
    measure_time(MPI_Barrier(MPI_COMM_WORLD));
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
    measure_time(construct_TSQR_Q(id, m_local, n_part, Y, Y_heads, tau, &Q));

    measure_time(MPI_Barrier(MPI_COMM_WORLD));
    double *S; // = malloc(n * sizeof(double));
    measure_time(modified_LU_decomposition(id, m_local, n, Q, &S));

    double *T;
    if (id == 0)
    {

        double *SY1;
        measure_time(SY1 = calloc(n * n, sizeof(double)));
        measure_time(T = calloc(n * n, sizeof(double)));
        measure_time(LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'L', n, n, Q, m_local, SY1, n));
        measure_time(LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'U', n, n, Q, m_local, T, n));

        measure_time(LAPACKE_dtrtri(LAPACK_COL_MAJOR, 'L', 'U', n, SY1, n));
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                T[i + j * n] *= S[j];
                R[i + j * n_part] *= S[i];
            }
        }

        measure_time(cblas_dtrmm(CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, n, n, -1.0, SY1, n, T, n));
        measure_time(free(SY1));
        measure_time(LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'U', n, n, R, n, Q, m_local));
    }
    measure_time(MPI_Barrier(MPI_COMM_WORLD));

    int desc[DESC_LEN];
    int ierror;
    Matrix mat;
    mat.data = Q;
    measure_time(exchange2(m, n, m_local, &mat, id));
    int ld;
    measure_time(ld = numroc_(&m, ADDR(int, m / proc_num), &rank, ADDR(int, 0), &proc_num));

    measure_time(descinit_(&desc, ADDR(int, m), ADDR(int, n), ADDR(int, m / proc_num), ADDR(int, n), ADDR(int, 0), ADDR(int, 0), &icontext_1d, ADDR(int, ld), &ierror));
    mat.desc = desc;

    measure_time(pdgemr2d_wrap(m, n, &mat, 0, 0, A, row, col));

    measure_time(MPI_Barrier(MPI_COMM_WORLD));

    // printf("[%d] descinit 2 %d %d %d %d %d %d\n", rank, n, n, n, n, icontext_1d, n);
    measure_time(descinit_(&desc, ADDR(int, n), ADDR(int, n), ADDR(int, n), ADDR(int, n), ADDR(int, 0), ADDR(int, 0), &icontext_1d, ADDR(int, n), &ierror));
    // printf("[%d] done 2 %d %d %d %d %d %d\n", rank, n, n, n, n, icontext_1d, n);
    mat.data = T;
    mat.desc = desc;

    measure_time(pdgemr2d_wrap(n, n, &mat, 0, 0, T_ret, 0, 0));
    measure_time(MPI_Barrier(MPI_COMM_WORLD));

    return;
}

int main(int argc, char **argv)
{

    // init
    MPI_Init(&argc, &argv);
    if (argc != 4)
    {
        printf("Usage %s matrix_size band_size qr_type('house' or 'tsqr')\n", argv[0]);
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
    int L = atoi(argv[2]);

    blacs_pinfo_(&rank, &nproc);
    rank_global = rank;
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
    if (argv[3][0] == 't')
    {
        rprintf("# use_tsqr\n");
        use_tsqr = true;
    }
    else
    {
        rprintf("# use householder\n");
        use_tsqr = false;
    }

    // main calc
    // printf("%d %d", m, n);
    if (rank == 0 && print_checkcode)
    {
        printf("import numpy as np\n");
    }
    n = 4;
    Matrix *A = create_matrix(proc_row_num, proc_col_num, m, m, block_row, block_col);
    Matrix *T = create_matrix(proc_row_num, proc_col_num, L, m, block_row, block_col);
    Matrix *Y = create_matrix(proc_row_num, proc_col_num, m, m, block_row, block_col);
    measure_time(for (size_t i = 0; i < A->global_row; ++i) {
        for (size_t j = 0; j < A->global_col; ++j)
        {
            double r = (double)(rand()) / RAND_MAX;
            set(A, i, j, r);
            set(A, j, i, r);
            //  set(A, i, j, 1.0);
        }
    });
    double norm_a = pdlange_wrap(A);

    blacs_barrier_(&icontext_2d, ADDR(char, 'A'));
    if (print_checkcode)
    {
        rprintf("import numpy as np\n");

        print_matrix("A=", A, rank);
    }
    measure_time(bischof(rank, proc_row_num, proc_col_num, m, L, A, T, Y));
    if (true)
    {
        for (int i = 0; i < m; i++)
        {
            for (int j = i; j < m; j++)
            {
                set(A, i, j, get(A, j, i));
            }
        }

        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < m; j++)
            {
                if (abs(i - j) > L)
                {
                    set(A, i, j, 0.0);
                }
            }
        }
    }
    if (print_checkcode)
    {
        print_matrix("B=", A, rank);
        rprintf("\nA = np.matrix(A)\nB = np.matrix(B)\ne=0\nfor i, j in zip(sorted(np.linalg.eigvals(A)), sorted(np.linalg.eigvals(B))):\n    print(i, j)\n    e+=(i-j)**2\n\nprint('error=',e**0.5)\n");
        print_matrix("T=", T, rank);
        print_matrix("Y=", Y, rank);
        rprintf("L = %d\n", L);
        rprintf("N = %d\n", m);
        rprintf("T=np.matrix(T)\n");
        rprintf("Y=np.matrix(Y)\n");
        rprintf("assert(N%%L==0)\n");
        rprintf("for i in range(N//L-1):\n");
        rprintf("    t=T[:L,L*i:L*(i+1)]\n");
        rprintf("    y=Y[:N-L*(i+1),L*i:L*(i+1)]\n");
        rprintf("    q=y*t*y.T\n");
        rprintf("    q=np.eye(q.shape[0])-q\n");
        rprintf("    print(np.linalg.norm(np.eye(q.shape[0])-q*q.T))\n");
    }
    double norm_b = pdlange_wrap(A);

    rprintf("# norm of A = %.18lf\n", norm_a);
    rprintf("# norm of B = %.18lf\n", norm_b);

    MPI_Barrier(MPI_COMM_WORLD);
    // print_matrix("A=", A, rank);
    MPI_Barrier(MPI_COMM_WORLD);

    // TSQR_HR(rank, my_row, my_col, m, n, A, 0, 0, T);
    MPI_Barrier(MPI_COMM_WORLD);
    // print_matrix("ret=", A, rank);
    // print_matrix("T=", T, rank);

    free_matrix(A);

    MPI_Finalize();
    if (rank == 0)
        dumpBuffer();
    return 0;
}
