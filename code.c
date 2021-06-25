#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <omp.h> // for a timing routine.
#include <cblas.h>
#include <lapacke.h>
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
            printf("%lf ,", a[j + lda * i]);
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

int tsqr(int matrix_layout, int m, int n, double *a, int lda, double *t, int ldt)
{
    if (m <= n)
    {
        LAPACKE_dgeqrt3(matrix_layout, m, n, a, lda, t, ldt);
        return 0;
    }
    int m_top = m / 2;
    int m_bottom = m - m_top;

    double *t_top = calloc(n * ldt, sizeof(double));
    double *t_bottom = calloc(n * ldt, sizeof(double));
    double *a_org = calloc(m * lda, sizeof(double));
    LAPACKE_dlacpy(matrix_layout, 'A', m, n, a, lda, a_org, lda);

    double *a_top = a;
    double *a_bottom = &a[m_top * lda];
    // print_matrix(m_top, n, a_top, lda);
    // print_matrix(m_bottom, n, a_bottom, lda);

    LAPACKE_dgeqrt3(matrix_layout, m_bottom, n, a_bottom, lda, t_bottom, ldt);
    LAPACKE_dgeqrt3(matrix_layout, m_top, n, a_top, lda, t_top, ldt);
    print_matrix("a=", m, n, a, lda);

    // TODO: forの範囲を変えて高速化
    for (int i = 0; i < m_top; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i > j)
            {
                a_top[j + i * lda] = 0.0;
            }
        }
    }

    // TODO: forの範囲を変えて高速化
    for (int i = 0; i < m_bottom; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i > j)
            {
                a_bottom[j + i * lda] = 0.0;
            }
        }
    }
    // (R1;R2) -> Q'R
    LAPACKE_dgeqrt3(matrix_layout, m, n, a, lda, t, ldt);
    // a_org=a_org-R
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i <= j)
            {
                a_org[j + i * lda] -= a[j + i * lda];
            }
        }
    }
    // aは上三角行列の形になっている
    double *R = a;
    // t = R_1()
    // TODO: LAPACKEのルーチンを使ってできないか考える
    for (int i = 0; i < n; i++)
    {
        for (int j = i; j < n; j++)
        {
            t[j + i * lda] = R[j + i * lda];
        }
    }

    // a_org=LU
    int *ipiv = calloc(n, sizeof(int));
    assert(LAPACKE_dgetrf(matrix_layout, m, n, a_org, lda, ipiv) == 0);

    for (int i = 0; i < n; i++)
    {
        // TODO:行列の置換が起きているときも動作するようにする
        assert(ipiv[i] == i + 1);
    }
    print_matrix("a_org= ", m, n, a_org, lda);
    // print_matrix(n, n, t, ldt);
    // t = t^-1 = R_1^-1

    LAPACKE_dtrtri(matrix_layout, 'U', 'N', n, t, ldt);

    print_matrix("Line 107", n, n, t, ldt);
    double *U = calloc(n * ldt, sizeof(double));

    // TODO: LAPACKEのルーチンを使ってできないか考える
    for (int i = 0; i < n; i++)
    {
        for (int j = i; j < n; j++)
        {
            U[j + i * lda] = a_org[j + i * lda];
        }
    }
    // LAPACKE_dlacpy(matrix_layout, 'U', n, n, a_org, lda, U, ldt);

    double *Yinv = calloc(n * ldt, sizeof(double));
    // Yinv = L
    print_matrix("Yinv 0", n, n, Yinv, ldt);

    // Yinv = a_orgの下三角部分
    // TODO: LAPACKEのルーチンを使ってできないか考える
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            Yinv[j + i * lda] = a_org[j + i * lda];
        }
    }
    print_matrix("Yinv 1", n, n, Yinv, ldt);

    for (int i = 0; i < n; i++)
    {
        Yinv[i + i * ldt] = 1.0;
    }
    print_matrix("Yinv 2", n, n, Yinv, ldt);

    // Yinv = Yinv^-1 = L^-1
    LAPACKE_dtrtri(matrix_layout, 'L', 'U', n, Yinv, ldt);
    print_matrix("Yinv 3", n, n, Yinv, ldt);

    // t2 =  R_1^-1 * Yinv^T
    double *t2 = calloc(n * ldt, sizeof(double));

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, n, 1.0, t, ldt, Yinv, ldt, 0.0, t2, ldt);
    print_matrix("Line 123", n, n, t2, ldt);

    // t = -Ut2
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, -1.0, U, ldt, t2, ldt, 0.0, t, ldt);
    print_matrix("Line 127", n, n, t, ldt);

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i > j)
            {
                a[j + i * lda] = a_org[j + i * lda];
            }
        }
    }
    free(U);
    free(t2);
    free(Yinv);
    free(ipiv);
    free(t_top);
    free(t_bottom);
    free(a_org);

    return 0;
}

int main(void)
{
    size_t const m = 8;
    size_t const n = 3;

    unsigned long const random_seed = 10;

    size_t const lda = n + 4;
    size_t const ldt = n + 4;

    double *const a = malloc(sizeof(double) * m * lda);
    double *const t = malloc(sizeof(double) * n * ldt);
    srand(random_seed);
    for (size_t i = 0; i < m; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            a[j + lda * i] = (double)(rand()) / RAND_MAX;
            // printf("%ld %ld\n", i, j);
        }
    }
    printf("test\naaaaaa");
    print_matrix("Aasasas= ", m, n, a, lda);
    double const t1 = omp_get_wtime();
    /* ----------------- START ----------------- */
    tsqr(LAPACK_ROW_MAJOR, m, n, a, lda, t, ldt);
    /* ----------------- END ----------------- */
    double const t2 = omp_get_wtime();
    print_matrix("res=", m, n, a, lda);
    //    print_matrix("t=", n, n, t, ldt);
    construct_Q_from_compact_WY(m, n, a, t, lda, ldt);

    srand(random_seed);
    for (size_t i = 0; i < m; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            a[j + lda * i] = (double)(rand()) / RAND_MAX;
        }
    }
    printf("start next\n");
    LAPACKE_dgeqrt3(LAPACK_ROW_MAJOR, m, n, a, lda, t, ldt);
    print_matrix("res=", m, n, a, lda);
    print_matrix("t=", n, n, t, ldt);
    construct_Q_from_compact_WY(m, n, a, t, lda, ldt);

    printf("Time [sec.]: %e\n", t2 - t1);
    printf("Performance [GFLOPS]: %e\n", (2.0 / 3.0) * n * n * n / (t2 - t1) * 1.0e-9);

    return EXIT_SUCCESS;
}
