#include <assert.h>
#include <mkl_cblas.h>
#include <mkl_lapacke.h>
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
            printf("%12.8lf ,", a[i + lda * j]);
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
            norm += a[i + lda * j] * a[i + lda * j];
        }
    }
    return sqrt(norm);
}

// LAPACKE_dgeqrtを行った後のA,TからQを陽に構築して、Qの先頭アドレスを返す。A,Tの内容は変更しない。
// cf. http://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_gaddcf152e87deec6123a1899f6f51101e.html
double *construct_Q(int m, int n, double *A, int lda, double *T, int ldt)
{
    // printf("in construct Q\n");
    // Q = I - VTV^t と表すことができ、Aの下三角部分からVが構築できる
    double *V = malloc(sizeof(double) * m * n);
    double *VT = malloc(sizeof(double) * m * n);
    double *Q = malloc(sizeof(double) * m * m);
    // print_matrix("(costruct) A = ", m, n, A, lda);

    // AからVを構築する。
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
                V[i + m * j] = 1.0;
            else if (i > j)
            {

                V[i + m * j] = A[i + lda * j];
            }
            else
            {
                V[i + m * j] = 0.0;
            }
        }
    }
    // print_matrix("(costruct) V = ", m, n, V, m);
    // print_matrix("(costruct) T = ", n, n, T, ldt);

    // VT = V * Tを計算
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, n, 1.0, V, m, T, ldt, 0.0, VT, m);
    // print_matrix("(costruct) VT = ", m, n, VT, m);

    // Q =I - VT * V^t を計算する
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, m, m, n, -1.0, VT, m, V, m, 0.0, Q, m);
    for (int i = 0; i < m; i++)
    {
        Q[i + m * i] += 1.0;
    }
    // print_matrix("(costruct) Q = ", m, m, Q, m);

    free(V);
    free(VT);
    // printf("out construct Q\n");

    return Q;
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
            assert(fabs(A[i + i * lda]) > 0.001);
            A[i + j * lda] = lu / A[i + i * lda];
        }
    }
}

// Aを半帯幅がLのN*N帯行列BとTによって構築し，その精度を調べる
void check(int N, int L, double *A, int lda, double *B, int ldb, double *T, int ldt)
{

    // print_matrix("T = ", L, N, T, ldt);
    // print_matrix("B = ", N, N, B, ldb);

    // Qを構築する
    double *tmp = malloc(N * N * sizeof(double));
    double *Q = malloc(N * N * sizeof(double));
    double *Qnext = malloc(N * N * sizeof(double));
    double *Qtmp = malloc(N * N * sizeof(double));

    // Qを単位行列で初期化する
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            Q[i + j * N] = (i == j ? 1.0 : 0.0);
        }
    }
    for (int k = 0; k < N / L - 1; k++)
    {
        int Nk = N - L - k * L;
        double *tmp = construct_Q(Nk, L, &B[(k + 1) * L + lda * k * L], lda, &T[k * ldt * L], ldt);
        // double *tmp2 = malloc(Nk * Nk * sizeof(double));
        // print_matrix("construct Q = ", Nk, Nk, tmp, Nk);
        // cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, Nk, Nk, Nk, 1.0, tmp, Nk, tmp, Nk, 0.0, tmp2, Nk);
        //// print_matrix("QQ^T = ", Nk, Nk, tmp2, Nk);
        // free(tmp2);

        // Qnext = [I O] とする
        //         [O Q]
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (j < L * (k + 1) || i < (L * (k + 1)))
                    Qnext[i + j * N] = (i == j ? 1.0 : 0.0);
                else
                    Qnext[i + j * N] = tmp[(i - L * (k + 1)) + (j - L * (k + 1)) * Nk];
            }
        }

        // print_matrix("tmp = ", N, N, Qnext, N);

        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, N, N, N, 1.0, Qnext, N, Qnext, N, 0.0, Qtmp, N);

        // print_matrix("QQ^T = ", N, N, Qtmp, N);

        for (int i = 0; i < N; i++)
        {
            Qtmp[i + i * N] -= 1.0;
        }
        //直交行列ならこの値が0になっているはず
        printf("norm of QQ^T - I = %e\n", calc_Frobenius_norm(N, N, Qtmp, N) / sqrt(N));

        // QをQ*Qnextで更新する
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, N, N, N, 1.0, Qnext, N, Q, N, 0.0, Qtmp, N);
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                Q[i + j * N] = Qtmp[i + j * N];
            }
        }
        free(tmp);
    }
    printf("N=%d\n", N);
    // print_matrix("Q = ", N, N, Q, N);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, N, N, N, 1.0, Q, N, Q, N, 0.0, Qtmp, N);
    for (int i = 0; i < N; i++)
    {
        Qtmp[i + i * N] -= 1.0;
    }
    //直交行列ならこの値が0になっているはず
    printf("norm of QQ^T - I = %e\n", calc_Frobenius_norm(N, N, Qtmp, N) / sqrt(N));
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (abs(i - j) > L)
            {
                B[j + i * ldb] = 0.0;
            }
        }
    }
    //QBQ^Tを計算することで，帯行列化したBから元の行列Aを復元する
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, N, N, N, 1.0, Q, N, B, ldb, 0.0, tmp, N);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, tmp, N, Q, N, 0.0, B, ldb);

    // 復元したaと元のaとの誤差を計算する
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j < N; ++j)
        {
            B[i + lda * j] -= A[i + j * lda];
        }
    }
    printf("norm  = %e\n", calc_Frobenius_norm(N, N, B, ldb) / calc_Frobenius_norm(N, N, A, lda));

    free(Qnext);
    free(Qtmp);
    free(tmp);
    free(Q);
}

int tsqr(int matrix_layout, int m, int n, double *a, int lda, double *t, int ldt)
{
    if (m <= n)
    //if (1 <= 1)
    {
        LAPACKE_dgeqrt(matrix_layout, m, n, n, a, lda, t, ldt);
        return 0;
    }
    int m_top = m / 2;
    int m_bottom = m - m_top;

    double *t_top = calloc(n * ldt, sizeof(double));
    double *t_bottom = calloc(n * ldt, sizeof(double));
    double *a_org = calloc(n * lda, sizeof(double));
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            a_org[i + j * lda] = a[i + j * lda];
        }
    }
    //    LAPACKE_dlacpy(matrix_layout, 'A', m, n, a, lda, a_org, lda);

    double *a_top = a;
    double *a_bottom = &a[m_top];
    //// print_matrix(m_top, n, a_top, lda);
    //// print_matrix(m_bottom, n, a_bottom, lda);

    LAPACKE_dgeqrt(matrix_layout, m_bottom, n, n, a_bottom, lda, t_bottom, ldt);
    LAPACKE_dgeqrt(matrix_layout, m_top, n, n, a_top, lda, t_top, ldt);
    // print_matrix("a=", m, n, a, lda);

    // TODO: forの範囲を変えて高速化
    for (int i = 0; i < m_top; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i > j)
            {
                a_top[i + j * lda] = 0.0;
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
                a_bottom[i + j * lda] = 0.0;
            }
        }
    }
    // (R1;R2) -> Q'R
    LAPACKE_dgeqrt(matrix_layout, m, n, n, a, lda, t, ldt);
    // a_org=a_org-R
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i <= j)
            {
                a_org[i + j * lda] -= a[i + j * lda];
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
            t[i + j * ldt] = R[i + j * lda];
        }
    }

    // a_org=LU
    print_matrix("a= ", m, n, a_org, lda);
    LU_decomp(m, n, a_org, lda);
    // int *ipiv = calloc(n, sizeof(int));
    // assert(LAPACKE_dgetrf(matrix_layout, m, n, a_org, lda, ipiv) == 0);
    // for (int i = 0; i < n; i++)
    // {
    //     // TODO:行列の置換が起きているときも動作するようにする
    //     assert(ipiv[i] == i + 1);
    // }
    print_matrix("lu= ", m, n, a_org, lda);
    //// print_matrix(n, n, t, ldt);
    // t = t^-1 = R_1^-1

    LAPACKE_dtrtri(matrix_layout, 'U', 'N', n, t, ldt);

    // print_matrix("Line 107", n, n, t, ldt);
    double *U = calloc(n * ldt, sizeof(double));

    // TODO: LAPACKEのルーチンを使ってできないか考える
    for (int i = 0; i < n; i++)
    {
        for (int j = i; j < n; j++)
        {
            U[i + j * ldt] = a_org[i + j * lda];
        }
    }
    // LAPACKE_dlacpy(matrix_layout, 'U', n, n, a_org, lda, U, ldt);

    double *Yinv = calloc(n * ldt, sizeof(double));
    // Yinv = L
    // print_matrix("Yinv 0", n, n, Yinv, ldt);

    // Yinv = a_orgの下三角部分
    // TODO: LAPACKEのルーチンを使ってできないか考える
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            Yinv[i + j * ldt] = a_org[i + j * lda];
        }
    }
    // print_matrix("Yinv 1", n, n, Yinv, ldt);

    for (int i = 0; i < n; i++)
    {
        Yinv[i + i * ldt] = 1.0;
    }
    // print_matrix("Yinv 2", n, n, Yinv, ldt);

    // Yinv = Yinv^-1 = L^-1
    LAPACKE_dtrtri(matrix_layout, 'L', 'U', n, Yinv, ldt);
    // print_matrix("Yinv 3", n, n, Yinv, ldt);

    // t2 =  R_1^-1 * Yinv^T
    double *t2 = calloc(n * ldt, sizeof(double));

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, n, 1.0, t, ldt, Yinv, ldt, 0.0, t2, ldt);
    // print_matrix("Line 123", n, n, t2, ldt);

    // t = -Ut2
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, -1.0, U, ldt, t2, ldt, 0.0, t, ldt);
    // print_matrix("Line 127", n, n, t, ldt);

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i > j)
            {
                a[i + j * lda] = R[i + j * lda];
            }
        }
    }
    free(U);
    free(t2);
    free(Yinv);
    //free(ipiv);
    free(t_top);
    free(t_bottom);
    free(a_org);

    return 0;
}

// (N,N)次元の実対称行列A(leading dimensionはlda)を帯行列化し、変換後の帯行列BとB=QAQ^tなるQを求める。関数終了時にAの要素はBの要素で上書きされる
void bischof(int N, int L, double *A, int lda, double *T, int ldt)
{
    int nb = L;
    int ldt_iter = nb;
    assert(MIN(N, L) >= nb && nb >= 1);
    assert(ldt_iter >= nb);
    assert(N % L == 0);
    double *const T_iter = calloc(ldt_iter * L, sizeof(double));

    for (int k = 0; k < N / L - 1; k++)
    {
        printf("iteration %d/%d\n", k + 1, N / L - 1);
        int Nk = N - L - k * L;
        // Aの(k+1,k)ブロック以下をQR分解する
        //LAPACKE_dgeqrt(LAPACK_COL_MAJOR, Nk, L, L, &A[(k + 1) * L + lda * k * L], lda, T_iter, ldt_iter);
        tsqr(LAPACK_COL_MAJOR, Nk, L, &A[(k + 1) * L + lda * k * L], lda, T_iter, ldt_iter);
        //// print_matrix("T iter = ", L, L, T_iter, ldt_iter);
        // TにT_iterを代入
        // print_matrix("T iter = ", L, L, T_iter, ldt_iter);

        for (int i = 0; i < L; i++)
        {
            for (int j = 0; j < L; j++)
            {
                if (j >= i)
                    T[i + (j + L * k) * ldt] = T_iter[i + j * ldt_iter];
                else
                    T[i + (j + L * k) * ldt] = T_iter[i + j * ldt_iter] = 0.0;
            }
        }
        // print_matrix("T = ", L, N, T, ldt);

        // Aを更新する
        double *V = malloc(Nk * L * sizeof(double));
        double *update_tmp = malloc(Nk * Nk * sizeof(double));
        double *update_P = malloc(Nk * L * sizeof(double));
        double *update_beta = malloc(L * L * sizeof(double));
        double *update_Q = malloc(Nk * L * sizeof(double));

        // Aの(k+1,k)ブロックの先頭アドレス
        double *a_part = &A[(k + 1) * L + lda * k * L];
        // print_matrix("A part = ", Nk, L, a_part, lda);

        // Vにa_partの下三角部分を代入する(LAPACKE_dgeqrtを実行しているので，a_partの下三角部分にはQのcompact-WY表現の一部が入っている)
        for (int i = 0; i < Nk; i++)
        {
            for (int j = 0; j < L; j++)
            {
                if (i == j)
                {
                    V[i + j * Nk] = 1.0;
                }
                else if (i > j)
                {
                    V[i + j * Nk] = a_part[i + j * lda];
                }
                else
                {
                    V[i + j * Nk] = 0.0;
                }
            }
        }

        // a_partの下三角部分を0クリアし，上三角部分の要素の値を，その要素の位置と対称な位置へ代入する（Aは対称行列なので）
        for (int i = L * k + L; i < N; i++)
        {
            for (int j = L * k; j < L * (k + 1); j++)
            {
                A[j + i * lda] = A[i + j * lda];
            }
        }
        // 山本有作先生の『キャッシュマシン向け三重対角化アルゴリズムの性能予測方式』で説明されているBischofのアルゴリズム中のPを構築する
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Nk, L, Nk, 1.0, &A[(k + 1) * L + lda * (k + 1) * L], lda, V, Nk, 0.0, update_tmp, Nk);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Nk, L, L, 1.0, update_tmp, Nk, T_iter, ldt_iter, 0.0, update_P, Nk);

        // betaを構築する
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, L, L, Nk, 0.5, V, Nk, update_P, Nk, 0.0, update_tmp, Nk);
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, L, L, L, 1.0, T_iter, ldt_iter, update_tmp, Nk, 0.0, update_beta, L);

        // Qを構築する
        for (int i = 0; i < Nk; i++)
        {
            for (int j = 0; j < L; j++)
            {
                update_Q[i + j * Nk] = update_P[i + j * Nk];
            }
        }
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Nk, L, L, -1.0, V, Nk, update_beta, L, 1.0, update_Q, Nk);

        // Aをrank-2k更新する
        cblas_dsyr2k(CblasColMajor, CblasUpper, CblasNoTrans, Nk, L, -1.0, V, Nk, update_Q, Nk, 1.0, &A[(k + 1) * L + lda * (k + 1) * L], lda);

        for (int i = 0; i < N; i++)
        {
            for (int j = i; j < N; j++)
            {
                A[j + i * lda] = A[i + j * lda];
            }
        }

        free(update_P);
        free(update_beta);
        free(update_Q);
        free(update_tmp);
        free(V);
    }
    free(T_iter);
}

int main(int argc, char **argv)
{
    unsigned long const random_seed = 100;
    srand(random_seed);

    // size_t const m = 10000; // 帯行列化する行列のサイズ
    size_t const m = atoi(argv[1]); // 帯行列化する行列のサイズ
    size_t const lda = m + 1;
    int L = 3; //帯行列化時の幅

    double *const a = malloc(sizeof(double) * m * lda); // 帯行列化する行列
    double *const b = malloc(sizeof(double) * m * lda); // aのコピーを格納する行列
    // double *Q = malloc(m * m * sizeof(double));           // Bischofのアルゴリズムで計算されるQを格納する行列
    double *T = calloc(L * m, sizeof(double));            // Bischofのアルゴリズムで計算されるQを格納する行列
    double *const tmp = malloc(sizeof(double) * m * lda); //計算に使う一時的な行列

    // aの値を0~1の乱数で初期化する
    for (size_t i = 0; i < m; ++i)
    {
        for (size_t j = 0; j < m; ++j)
        {
            a[i + j * lda] = a[j + i * lda] = (double)(rand()) / RAND_MAX;
        }
    }

    // aのコピーをbに格納する
    for (size_t i = 0; i < m; ++i)
    {
        for (size_t j = 0; j < m; ++j)
        {
            b[i + lda * j] = a[i + j * lda];
        }
    }

    //// print_matrix("A= ", m, m, a, lda);
    // aのノルムを計算しておく
    //    double norm_a = calc_Frobenius_norm(m, m, a, lda);
    //Bischofのアルゴリズム実行
    double t1 = omp_get_wtime();
    bischof(m, L, a, lda, T, L);
    printf("time : %lf [sec]\n", omp_get_wtime() - t1);
    // print_matrix("A = ", m, m, b, lda);

    // print_matrix("B = ", m, m, a, lda);

    // bはaのコピーなので，この順番が正しい

    //check(m, L, b, lda, a, lda, T, L);

    //Qが直交行列かQQ^Tを計算することで確かめる
    // cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, m, m, m, 1.0, Q, m, Q, m, 0.0, tmp, lda);

    // //直交行列ならこの値が1になっているはず
    // printf("norm of QQ^T = %e\n", calc_Frobenius_norm(m, m, tmp, lda) / sqrt(m));

    // //QaQ^Tを計算することで，帯行列化したaから元の行列aを復元する
    // cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, m, m, m, 1.0, Q, m, a, lda, 0.0, tmp, lda);
    // cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, m, m, 1.0, tmp, lda, Q, m, 0.0, a, lda);

    // // 復元したaと元のaとの誤差を計算する
    // for (size_t i = 0; i < m; ++i)
    // {
    //     for (size_t j = 0; j < m; ++j)
    //     {
    //         b[i + lda * j] -= a[i + j * lda];
    //     }
    // }
    // printf("norm  = %e\n", calc_Frobenius_norm(m, m, b, lda) / norm_a);

    printf("END\n");

    free(a);
    free(b);
    free(tmp);
    free(T);

    return EXIT_SUCCESS;
}
