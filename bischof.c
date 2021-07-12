#include <assert.h>
#include <cblas.h>
#include <lapacke.h>
#include <math.h>
#include <omp.h> // for a timing routine.
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define MIN(a, b) ((a) < (b) ? (a) : (b))

// (m,n)次元でleading dimensionがldaの行列aの要素を，msgで指定された文字列の後に出力する．（デバッグ用）
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

// (m,n)次元、leading dimensionがldaの行列aのFrobenius normを計算して返す
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

// LAPACKE_dgeqrtを行った後のA,TからQを陽に構築して、Qの先頭アドレスを返す。A,Tの内容は変更しない。
// cf. http://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_gaddcf152e87deec6123a1899f6f51101e.html
double *construct_Q(int m, int n, double *A, int lda, double *T, int ldt)
{
    printf("in construct Q\n");
    // Q = I - VTV^t と表すことができ、Aの下三角部分からVが構築できる
    double *V = malloc(sizeof(double) * m * n);
    double *VT = malloc(sizeof(double) * m * n);
    double *Q = malloc(sizeof(double) * m * m);
    // AからVを構築する。
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
    // VT = V * Tを計算
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, n, 1.0, V, n, T, ldt, 0.0, VT, n);

    // Q =I - VT * V^t を計算する
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

// (N,N)次元の実対称行列A(leading dimensionはlda)を帯行列化し、変換後の帯行列BとB=QAQ^tなるQを求める。関数終了時にAの要素はBの要素で上書きされる
void bischof(int matrix_layout, int N, double *A, int lda, double *Q)
{
    int L = 50; //帯行列化時の幅
    int nb = L;
    int ldt = L;
    assert(MIN(N, L) >= nb && nb >= 1);
    assert(ldt >= nb);
    assert(N % L == 0);
    // Qの計算に用いる領域
    double *Qnext = malloc(N * N * sizeof(double));
    double *Qtmp = malloc(N * N * sizeof(double));

    // Qを単位行列で初期化する
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            Q[j + i * N] = (i == j ? 1.0 : 0.0);
        }
    }

    for (int k = 0; k < N / L - 1; k++)
    {

        double *const t = calloc(ldt * nb, sizeof(double));
        int Nk = N - L - k * L;
        // Aの(k+1,k)ブロック以下をQR分解する
        LAPACKE_dgeqrt(LAPACK_ROW_MAJOR, Nk, L, L, &A[k * L + lda * (k + 1) * L], lda, t, ldt);
        // tの下三角部分を0クリア
        for (int i = 0; i < L; i++)
        {
            for (int j = 0; j < i; j++)
            {
                t[j + i * ldt] = 0.0;
            }
        }

        // Qを構築する
        double *tmp = construct_Q(Nk, L, &A[k * L + lda * (k + 1) * L], lda, t, ldt);
        // Qnext = [I O] とする
        //         [O Q]
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

        // QをQ*Qnextで更新する
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, N, N, N, 1.0, Qnext, N, Q, N, 0.0, Qtmp, N);
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                Q[j + i * N] = Qtmp[j + i * N];
            }
        }
        free(tmp);

        // Aを更新する
        double *V = malloc(Nk * L * sizeof(double));
        double *update_tmp = malloc(Nk * Nk * sizeof(double));
        double *update_P = malloc(Nk * L * sizeof(double));
        double *update_beta = malloc(L * L * sizeof(double));
        double *update_Q = malloc(Nk * L * sizeof(double));

        // Aの(k+1,k)ブロックの先頭アドレス
        double *a_part = &A[k * L + lda * (k + 1) * L];

        // Vにa_partの下三角部分を代入する(LAPACKE_dgeqrtを実行しているので，a_partの下三角部分にはQのcompact-WY表現の一部が入っている)
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

        // a_partの下三角部分を0クリアし，上三角部分の要素の値を，その要素の位置と対称な位置へ代入する（Aは対称行列なので）
        for (int i = L * k + L; i < N; i++)
        {
            for (int j = L * k; j < L * (k + 1); j++)
            {
                if (i > j + L)
                {
                    A[j + i * lda] = A[i + j * lda] = 0.0;
                }
                else
                {
                    A[i + j * lda] = A[j + i * lda];
                }
            }
        }

        // 山本有作先生の『キャッシュマシン向け三重対角化アルゴリズムの性能予測方式』で説明されているBischofのアルゴリズム中のPを構築する
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nk, L, Nk, 1.0, &A[(k + 1) * L + lda * (k + 1) * L], lda, V, L, 0.0, update_tmp, Nk);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nk, L, L, 1.0, update_tmp, Nk, t, ldt, 0.0, update_P, L);

        // betaを構築する
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, L, L, Nk, 0.5, V, L, update_P, L, 0.0, update_tmp, Nk);
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, L, L, L, 1.0, t, ldt, update_tmp, Nk, 0.0, update_beta, L);

        // Qを構築する
        for (int i = 0; i < Nk; i++)
        {
            for (int j = 0; j < L; j++)
            {
                update_Q[j + i * L] = update_P[j + i * L];
            }
        }
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nk, L, L, -1.0, V, L, update_beta, L, 1.0, update_Q, L);

        // Aをrank-2k更新する
        cblas_dsyr2k(CblasRowMajor, CblasUpper, CblasNoTrans, Nk, L, -1.0, V, L, update_Q, L, 1.0, &A[(k + 1) * L + lda * (k + 1) * L], lda);
        for (int i = 0; i < N; i++)
        {
            for (int j = i; j < N; j++)
            {
                A[i + j * lda] = A[j + i * lda];
            }
        }

        // AをDGEMMを使って更新する
        // cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Nk, Nk, Nk, 1.0, tmp, Nk, &a[(k + 1) * L + lda * (k + 1) * L], lda, 0.0, tmp2, Nk);
        // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nk, Nk, Nk, 1.0, tmp2, Nk, tmp, Nk, 0.0, &a[(k + 1) * L + lda * (k + 1) * L], lda);

        free(update_P);
        free(update_beta);
        free(update_Q);
        free(update_tmp);
        free(V);
        free(t);
    }
    free(Qtmp);
    free(Qnext);
}

int main(void)
{
    unsigned long const random_seed = 100;
    srand(random_seed);

    size_t const m = 3000; // 帯行列化する行列のサイズ
    size_t const lda = m + 1;

    double *const a = malloc(sizeof(double) * m * lda);   // 帯行列化する行列
    double *const b = malloc(sizeof(double) * m * lda);   // aのコピーを格納する行列
    double *Q = malloc(m * m * sizeof(double));           // Bischofのアルゴリズムで計算されるQを格納する行列
    double *const tmp = malloc(sizeof(double) * m * lda); //計算に使う一時的な行列

    // aの値を0~1の乱数で初期化する
    for (size_t i = 0; i < m; ++i)
    {
        for (size_t j = 0; j < m; ++j)
        {
            a[j + lda * i] = a[i + j * lda] = (double)(rand()) / RAND_MAX;
        }
    }

    // aのコピーをbに格納する
    for (size_t i = 0; i < m; ++i)
    {
        for (size_t j = 0; j < m; ++j)
        {
            b[j + lda * i] = a[j + i * lda];
        }
    }

    // print_matrix("A= ", m, m, a, lda);
    // aのノルムを計算しておく
    double norm_a = calc_Frobenius_norm(m, m, a, lda);

    //Bischofのアルゴリズム実行
    bischof(LAPACK_ROW_MAJOR, m, a, lda, Q);

    //Qが直交行列かQQ^Tを計算することで確かめる
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, m, m, 1.0, Q, m, Q, m, 0.0, tmp, lda);
    //直交行列ならこの値が1になっているはず
    printf("norm of QQ^T = %e\n", calc_Frobenius_norm(m, m, tmp, lda) / sqrt(m));

    //QaQ^Tを計算することで，帯行列化したaから元の行列aを復元する
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m, m, m, 1.0, Q, m, a, lda, 0.0, tmp, lda);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, m, m, 1.0, tmp, lda, Q, m, 0.0, a, lda);

    // 復元したaと元のaとの誤差を計算する
    for (size_t i = 0; i < m; ++i)
    {
        for (size_t j = 0; j < m; ++j)
        {
            b[j + lda * i] -= a[j + i * lda];
        }
    }
    printf("norm  = %e\n", calc_Frobenius_norm(m, m, b, lda) / norm_a);

    printf("END\n");

    free(a);
    free(b);
    free(tmp);
    free(Q);
    return EXIT_SUCCESS;
}
