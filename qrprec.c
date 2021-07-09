#include <assert.h>
#include <cblas.h>
#include <lapacke.h>
#include <math.h>
#include <omp.h> // for a timing routine.
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
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
    print_matrix("A=", m, n, A, lda);
    //    print_matrix("V=", m, n, V, n);

    // print_matrix("t= ", L, L, t, L + 1);
    // print_matrix("v=", n, n, v, n);
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
    int n = 10;
    int lda = n + 3;
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

    double *t = malloc((ldt + 1) * n * sizeof(double));
    printf("%p\n", t);
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
    free(a);
    free(b);
    free(q);
    free(v);
    free(vt);
}

int main(void)
{
    unsigned long const random_seed = 10;
    //sranddev(); //srand(time(NULL));だと最初のrand()の返り値が偏る cf:https://stackoverflow.com/questions/32489058/trouble-generating-random-numbers-in-c-on-a-mac-using-xcode
    srand(random_seed);
    printf("%lf\n", (double)(rand()) / RAND_MAX);
    test();
    return 0;
}
