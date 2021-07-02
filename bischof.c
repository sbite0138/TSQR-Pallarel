#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <omp.h> // for a timing routine.
#include <cblas.h>
#include <lapacke.h>
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

// QR分解はLAPACKを使って大丈夫

//https://www.hpc.nec/documents/sdk/SDK_NLC/UsersGuide/man/dsyr2k.html
void bischof(int matrix_layout, int N, double *a, int lda)
{
    int L = 10;
    int nb = 2;
    assert(MIN(N, L) >= nb && nb >= 1);
    int ldt = nb;
    assert(ldt >= nb);

    assert(N % L == 0);
    for (int k = 0; k < N / L; k++)
    {

        double *const t = calloc(ldt * MIN(N, L), sizeof(double));
        int Nk = (N / L - k) * L;
        printf("%d\n", ldt);
        LAPACKE_dgeqrt(LAPACK_ROW_MAJOR, N, L, nb, a, lda, t, ldt);
        // print_matrix("t= ", L, L, t, L + 1);
        printf("ok %p\n", t);
        free(t);
    }
    printf("ok");
}

int main(void)
{
    size_t const m = 30;

    unsigned long const random_seed = 10;

    size_t const lda = m + 4;

    double *const a = malloc(sizeof(double) * m * lda);
    srand(random_seed);
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
