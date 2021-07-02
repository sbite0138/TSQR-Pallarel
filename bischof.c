#include <assert.h>
#include <cblas.h>
#include <lapacke.h>
#include <math.h>
#include <omp.h> // for a timing routine.
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define MIN(a, b) ((a) < (b) ? (a) : (b))
void print_matrix(char *msg, int m, int n, double *a, int lda)
{
    printf("================================\n");
    printf("%s\n", msg);
    printf("[\n");
    for (int i = 0; i < m; i++)
    {

        printf("[");
        for (int j = 0; j < n; j++)
        {
            printf("%12.8lf ,", a[j + lda * i]);
        }
        printf("],");
        printf("\n");
    }
    printf("]\n");
}

void construct_Q_from_compact_WY(int m, int n, double *Y, double *T, int ldy, int ldt)
{
    double *YT = calloc(m * n, sizeof(double));
    double *Q = calloc(m * m, sizeof(double));

    for (int i = 0; i < m; i++)
    {
        for (int j = i; j < n; j++)
        {
            if (i == j)
                Y[j + i * ldy] = 1.0;
            else
                Y[j + i * ldy] = 0.0;
        }
    }
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < i; j++)
        {
            T[j + i * ldt] = 0.0;
        }
    }
    print_matrix("Y=", m, n, Y, ldy);
    print_matrix("T=", n, n, T, ldt);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, n, 1.0, Y, ldy, T, ldt, 0.0, YT, n);
    print_matrix("YT=", m, n, YT, n);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, m, n, 1.0, YT, n, Y, ldy, 0.0, Q, m);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            Q[j + i * m] = (i == j ? 1.0 : 0.0) - Q[j + i * m];
        }
    }
    print_matrix("Q=", m, m, Q, m);
    free(YT);
    free(Q);
}

double calc_Frobenius_norm(int m, int n, double *a, int lda)
{
    double norm = 0.0;
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            norm += a[j + lda * i] * a[j + lda * i];
        }
    }
    return sqrt(norm);
}
double *gen_matrix(int m, int n, int lda)
{
    assert(lda >= n);
    double *a = malloc(sizeof(double) * m * lda);
    for (size_t i = 0; i < m; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            a[j + lda * i] = 1.0 - 2.0 * (double)(rand()) / RAND_MAX;
            // printf("%ld %ld\n", i, j);
        }
    }
    return a;
}
double *construct_Q(int m, int n, double *A, int lda, double *T, int ldt)
{
    double *V = malloc(sizeof(double) * m * n);
    double *VT = malloc(sizeof(double) * m * n);
    double *Q = malloc(sizeof(double) * m * m);

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
                V[j + n * i] = 1.0;
            else if (i > j)
            {

                V[j + n * i] = A[j + lda * i];
            }
            else
            {
                V[j + n * i] = 0.0;
            }
        }
    }
    print_matrix("A=", m, n, A, lda);
    print_matrix("V=", m, n, V, n);

    // print_matrix("t= ", L, L, t, L + 1);
    // print_matrix("v=", n, n, v, n);
    print_matrix("t=", n, n, T, ldt);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, n, 1.0, V, n, T, ldt, 0.0, VT, n);

    print_matrix("VT =", m, n, VT, ldt);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, m, n, -1.0, VT, n, V, n, 0.0, Q, m);
    for (int i = 0; i < m; i++)
    {
        Q[i + m * i] += 1.0;
    }
    free(V);
    free(VT);
    return Q;
}

void test()
{
    int m = 100;
    int n = 5;
    int lda = n + 1;
    int ldt = n + 3;
    double *a = gen_matrix(m, n, lda);
    double *b = malloc(sizeof(double) * m * n);
    double *v = malloc(sizeof(double) * m * n);
    double *vt = malloc(sizeof(double) * m * n);

    double *q; //= malloc(sizeof(double) * n * n);
    double norm_a = calc_Frobenius_norm(m, n, a, lda);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            b[j + n * i] = a[j + lda * i];
            v[j + n * i] = 0.0;
        }
    }

    // print_matrix("A = ", n, n, a, n);

    double *const t = calloc(ldt * n, sizeof(double));
    LAPACKE_dgeqrt3(LAPACK_ROW_MAJOR, m, n, a, lda, t, ldt);
    q = construct_Q(m, n, a, lda, t, ldt);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i > j)
            {

                a[j + lda * i] = 0;
            }
        }
    }
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, m, 1.0, q, m, a, lda, 0.0, vt, n);

    // print_matrix("t= ", L, L, t, L + 1);
    // print_matrix("v=", n, n, v, n);
    // print_matrix("t=", n, n, t, n);
    print_matrix("q=", m, m, q, m);
    print_matrix("A = ", m, n, b, n);
    print_matrix("qr=", m, n, vt, n);

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            b[j + n * i] -= vt[j + n * i];
        }
    }
    double norm_qr = calc_Frobenius_norm(m, n, b, n);
    printf("error = %e\n", norm_qr / norm_a);
    free(t);
}

// QR分解はLAPACKを使って大丈夫

//https://www.hpc.nec/documents/sdk/SDK_NLC/UsersGuide/man/dsyr2k.html
void bischof(int matrix_layout, int N, double *a, int lda)
{
    int L = 10;
    int nb = L;
    assert(L % nb == 0);
    assert(MIN(N, L) >= nb && nb >= 1);
    int ldt = L;
    assert(ldt >= nb);

    assert(N % L == 0);
    for (int k = 0; k < N / L - 1; k++)
    {

        double *const t = calloc(ldt * nb, sizeof(double));
        int Nk = (N / L - k + 1) * L;
        printf("%d\n", ldt);
        LAPACKE_dgeqrt(LAPACK_ROW_MAJOR, Nk, L, nb, &a[k * L + lda * k * L], lda, t, ldt);
        // print_matrix("t= ", L, L, t, L + 1);
        printf("ok %p\n", t);

        print_matrix("t=", nb, L, t, ldt);
        free(t);
    }
    printf("ok");
}

int main(void)
{
    unsigned long const random_seed = 10;
    printf("%ld\n", time(NULL));
    sranddev(); //srand(time(NULL));だと最初のrand()の返り値が偏る cf:https://stackoverflow.com/questions/32489058/trouble-generating-random-numbers-in-c-on-a-mac-using-xcode
    printf("%lf\n", (double)(rand()) / RAND_MAX);
    test();
    return 0;
    size_t const m = 30;

    size_t const lda = m + 4;

    double *const a = malloc(sizeof(double) * m * lda);
    for (size_t i = 0; i < m; ++i)
    {
        for (size_t j = 0; j < m; ++j)
        {
            a[j + lda * i] = (double)(rand()) / RAND_MAX;
            // printf("%ld %ld\n", i, j);
        }
    }
    print_matrix("A= ", m, m, a, lda);
    double const t1 = omp_get_wtime();
    /* ----------------- START ----------------- */
    bischof(LAPACK_ROW_MAJOR, m, a, lda);
    /* ----------------- END ----------------- */
    double const t2 = omp_get_wtime();
    print_matrix("res=", m, m, a, lda);
    //    print_matrix("t=", n, n, t, ldt);
    return EXIT_SUCCESS;
}
