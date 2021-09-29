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
            printf("%lf ,", a[i + lda * j]);
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
//cf: http://slis.tsukuba.ac.jp/~fujisawa.makoto.fu/cgi-bin/wiki/index.php?LU%CA%AC%B2%F2
void LU_decomp(int M, int N, double *A, int lda)
{

    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j <= i && j < N; j++)
        {
            double lu = A[i + j * lda];
            for (int k = 0; k < j; k++)
            {
                lu -= A[i + k * lda] * A[k + j * lda];
            }
            A[i + j * lda] = lu;
        }
        for (int j = i + 1; j < N; j++)
        {
            double lu = A[i + j * lda];
            for (int k = 0; k < i; k++)
            {
                lu -= A[i + k * lda] * A[k + j * lda];
            }
            A[i + j * lda] = lu / A[i + i * lda];
        }
    }
}
int main(void)
{
    size_t const m = 8;
    size_t const n = 4;

    unsigned long const random_seed = 10;

    size_t const lda = m + 4;
    size_t const ldt = n + 4;

    double *const a = malloc(sizeof(double) * n * lda);
    srand(random_seed);
    for (size_t i = 0; i < m; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            a[i + j * lda] = (double)(rand()) / RAND_MAX;
            // printf("%ld %ld\n", i, j);
        }
    }
    print_matrix("A= ", m, n, a, lda);
    double const t1 = omp_get_wtime();
    /* ----------------- START ----------------- */
    LU_decomp(m, n, a, lda);
    //    tsqr(LAPACK_ROW_MAJOR, m, n, a, lda, t, ldt);
    /* ----------------- END ----------------- */
    double const t2 = omp_get_wtime();
    print_matrix("res=", m, n, a, lda);
    //    print_matrix("t=", n, n, t, ldt);

    printf("Time [sec.]: %e\n", t2 - t1);
    printf("Performance [GFLOPS]: %e\n", (2.0 / 3.0) * n * n * n / (t2 - t1) * 1.0e-9);
    free(a);
    return EXIT_SUCCESS;
}
