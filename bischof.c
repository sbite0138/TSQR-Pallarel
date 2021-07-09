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
    return;
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
    // print_matrix("Y=", m, n, Y, ldy);
    // print_matrix("T=", n, n, T, ldt);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, n, 1.0, Y, ldy, T, ldt, 0.0, YT, n);
    // print_matrix("YT=", m, n, YT, n);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, m, n, 1.0, YT, n, Y, ldy, 0.0, Q, m);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            Q[j + i * m] = (i == j ? 1.0 : 0.0) - Q[j + i * m];
        }
    }
    // print_matrix("Q=", m, m, Q, m);
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
    printf("%d\n", sizeof(double) * m * lda);
    double *a = malloc(sizeof(double) * m * lda);

    for (size_t i = 0; i < m; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            a[j + lda * i] = 10.0 - 20.0 * (double)(rand()) / RAND_MAX;
            // printf("%ld %ld\n", i, j);
        }
    }
    return a;
}

double *construct_Q(int m, int n, double *A, int lda, double *T, int ldt)
{
    printf("in construct Q\n");
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
    // print_matrix("A=", m, n, A, lda);
    // print_matrix("V=", m, n, V, n);
    // print_matrix("t=", n, n, T, ldt);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, n, 1.0, V, n, T, ldt, 0.0, VT, n);

    // print_matrix("VT =", m, n, VT, n);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, m, n, -1.0, VT, n, V, n, 0.0, Q, m);
    for (int i = 0; i < m; i++)
    {
        Q[i + m * i] += 1.0;
    }
    free(V);
    free(VT);
    printf("out construct Q\n");

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
    LAPACKE_dgeqrt(LAPACK_ROW_MAJOR, m, n, n, a, lda, t, ldt);
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
    // print_matrix("q=", m, m, q, m);
    // print_matrix("A = ", m, n, b, n);
    // print_matrix("qr=", m, n, vt, n);

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
    free(a);
    free(b);
    free(q);
    free(v);
    free(vt);
}

// QR分解はLAPACKを使って大丈夫

//https://www.hpc.nec/documents/sdk/SDK_NLC/UsersGuide/man/dsyr2k.html
void bischof(int matrix_layout, int N, double *a, int lda, double *Q)
{
    int L = 500;
    int nb = L;
    assert(L % nb == 0);
    assert(MIN(N, L) >= nb && nb >= 1);
    int ldt = L;
    assert(ldt >= nb);
    double *Qnext = malloc(N * N * sizeof(double));
    double *Qtmp = malloc(N * N * sizeof(double));
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            Q[j + i * N] = (i == j ? 1.0 : 0.0);
        }
    }

    assert(N % L == 0);
    for (int k = 0; k < N / L - 1; k++)
    //for (int k = 0; k < 1; k++)
    {

        double *const t = calloc(ldt * nb, sizeof(double));
        int Nk = N - L - k * L;

        LAPACKE_dgeqrt(LAPACK_ROW_MAJOR, Nk, L, L, &a[k * L + lda * (k + 1) * L], lda, t, ldt);
        for (int i = 0; i < L; i++)
        {
            for (int j = 0; j < i; j++)
            {
                t[j + i * ldt] = 0.0;
            }
        } // print_matrix("qr a =", N, N, a, lda);

        /// construct Q
        double *tmp = construct_Q(Nk, L, &a[k * L + lda * (k + 1) * L], lda, t, ldt);
        //print_matrix("tmp = ", Nk, Nk, tmp, Nk);

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (i < L * (k + 1) || j < (L * (k + 1)))
                    Qnext[j + i * N] = (i == j ? 1.0 : 0.0);
                else
                    Qnext[j + i * N] = tmp[(j - L * (k + 1)) + (i - L * (k + 1)) * Nk];
            }
        }
        // print_matrix("Qnext = ", N, N, Qnext, N);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, Q, N, Qnext, N, 0.0, Qtmp, N);
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                Q[j + i * N] = Qtmp[j + i * N];
            }
        }
        // print_matrix("Q = ", N, N, Q, N);
        free(tmp);
        // construct Q done

        // update A
        double *V = malloc(Nk * L * sizeof(double));
        double *update_tmp = malloc(Nk * Nk * sizeof(double));
        double *update_P = malloc(Nk * L * sizeof(double));
        double *update_beta = malloc(L * L * sizeof(double));
        double *update_Q = malloc(Nk * L * sizeof(double));

        double *a_part = &a[k * L + lda * (k + 1) * L];
        for (int i = 0; i < Nk; i++)
        {
            for (int j = 0; j < L; j++)
            {
                if (i == j)
                {
                    V[j + i * L] = 1.0;
                }
                else if (i > j)
                {
                    V[j + i * L] = a_part[j + i * lda];
                }
                else
                {
                    V[j + i * L] = 0.0;
                }
            }
        }
        //print_matrix("V= ", Nk, L, V, L);

        for (int i = L * k + L; i < N; i++)
        {
            for (int j = L * k; j < L * (k + 1); j++)
            {
                if (i > j + L)
                {
                    a[j + i * lda] = a[i + j * lda] = 0.0;
                }
                else
                {
                    a[i + j * lda] = a[j + i * lda];
                }
            }
        }

        // print_matrix("sym a =", N, N, a, lda);
        // construct P
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nk, L, Nk, 1.0, &a[(k + 1) * L + lda * (k + 1) * L], lda, V, L, 0.0, update_tmp, Nk);
        // print_matrix("before V =", Nk, L, V, L);
        // print_matrix("before a[(k + 1) * L + lda * (k + 1) * L] =", Nk, Nk, &a[(k + 1) * L + lda * (k + 1) * L], lda);
        // print_matrix("before update_tmp =", Nk, L, update_tmp, Nk);
        // print_matrix("before t =", L, L, t, ldt);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nk, L, L, 1.0, update_tmp, Nk, t, ldt, 0.0, update_P, L);
        // print_matrix("before P =", Nk, L, update_P, L);

        // construct beta
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, L, L, Nk, 0.5, V, L, update_P, L, 0.0, update_tmp, Nk);
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, L, L, L, 1.0, t, ldt, update_tmp, Nk, 0.0, update_beta, L);
        // print_matrix("before beta =", L, L, update_beta, L);

        // construct Q
        for (int i = 0; i < Nk; i++)
        {
            for (int j = 0; j < L; j++)
            {
                update_Q[j + i * L] = update_P[j + i * L];
            }
        }
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nk, L, L, -1.0, V, L, update_beta, L, 1.0, update_Q, L);

        // print_matrix("before a =", N, N, a, lda);
        // print_matrix("before V =", Nk, L, V, L);
        // print_matrix("before Q =", Nk, L, update_Q, L);

        // rank-2k update
        cblas_dsyr2k(CblasRowMajor, CblasUpper, CblasNoTrans, Nk, L, -1.0, V, L, update_Q, L, 1.0, &a[(k + 1) * L + lda * (k + 1) * L], lda);

        // DGEMM update
        // cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Nk, Nk, Nk, 1.0, tmp, Nk, &a[(k + 1) * L + lda * (k + 1) * L], lda, 0.0, tmp2, Nk);
        //cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nk, Nk, Nk, 1.0, tmp2, Nk, tmp, Nk, 0.0, &a[(k + 1) * L + lda * (k + 1) * L], lda);

        for (int i = 0; i < N; i++)
        {
            for (int j = i; j < N; j++)
            {
                a[i + j * lda] = a[j + i * lda];
            }
        }

        // print_matrix("after a =", N, N, a, lda);
        free(update_P);
        free(update_beta);
        free(update_Q);
        free(update_tmp);
        //print_matrix("@A = ", N, N, a, lda);

        free(V);

        free(t);
    }
    free(Qtmp);
    free(Qnext);
}

int main(void)
{
    unsigned long const random_seed = 10;
    //sranddev(); //srand(time(NULL));だと最初のrand()の返り値が偏る cf:https://stackoverflow.com/questions/32489058/trouble-generating-random-numbers-in-c-on-a-mac-using-xcode
    srand(random_seed);
    printf("%lf\n", (double)(rand()) / RAND_MAX);
    //  test();
    //return 0;
    size_t const m = 3000;

    size_t const lda = m + 1;

    double *const a = malloc(sizeof(double) * m * lda);
    double *const b = malloc(sizeof(double) * m * lda);
    double *const tmp = malloc(sizeof(double) * m * lda);

    for (size_t i = 0; i < m; ++i)
    {
        for (size_t j = 0; j < m; ++j)
        {
            a[j + lda * i] = a[i + j * lda] = 1.0 * (double)(rand()) / RAND_MAX;
            // printf("%ld %ld\n", i, j);
        }
    }
    for (size_t i = 0; i < m; ++i)
    {
        for (size_t j = 0; j < m; ++j)
        {
            b[j + lda * i] = a[j + i * lda];
            // printf("%ld %ld\n", i, j);
        }
    }
    double *Q = malloc(m * m * sizeof(double));

    print_matrix("A= ", m, m, a, lda);
    double norm_a = calc_Frobenius_norm(m, m, a, lda);
    double const t1 = omp_get_wtime();
    /* ----------------- START ----------------- */
    bischof(LAPACK_ROW_MAJOR, m, a, lda, Q);
    /* ----------------- END ----------------- */
    double const t2 = omp_get_wtime();
    print_matrix("res=", m, m, a, lda);
    print_matrix("Q=", m, m, Q, m);
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m, m, m, 1.0, Q, m, Q, m, 0.0, tmp, lda);
    print_matrix("QQ^T =", m, m, tmp, lda);
    printf("norm of QQ^T = %e\n", calc_Frobenius_norm(m, m, tmp, lda) / sqrt(m));
    // cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m, m, m, 1.0, Q, m, b, lda, 0.0, a, lda);
    // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, m, m, 1.0, a, lda, Q, m, 0.0, b, lda);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, m, m, 1.0, Q, m, a, lda, 0.0, tmp, lda);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, m, m, 1.0, tmp, lda, Q, m, 0.0, a, lda);
    print_matrix("a =", m, m, a, lda);
    for (size_t i = 0; i < m; ++i)
    {
        for (size_t j = 0; j < m; ++j)
        {
            b[j + lda * i] -= a[j + i * lda];
        }
    }

    print_matrix("b =", m, m, b, lda);

    printf("norm  = %e\n", calc_Frobenius_norm(m, m, b, lda) / norm_a);
    //    print_matrix("t=", n, n, t, ldt);
    printf("END\n");

    free(a);
    free(b);
    free(tmp);
    free(Q);
    return EXIT_SUCCESS;
}
