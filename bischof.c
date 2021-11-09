#include <assert.h>
#define ON_LOCAL
#ifdef ON_LOCAL
#include <cblas.h>
#include <lapacke.h>
#else
#include <mkl_cblas.h>
#include <mkl_lapacke.h>
#endif
#include <math.h>
#include <omp.h> // for a timing routine.
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define measure_time(x)                 \
    do                                  \
    {                                   \
        double start = omp_get_wtime(); \
        {                               \
            x;                          \
        }                               \
        double end = omp_get_wtime();   \
    } while (0)
// (m,n)次元でleading dimensionがldaの行列aの要素を，msgで指定された文字列の後に出力する．（デバッグ用）
void print_matrix(char *msg, int m, int n, double *a, int lda)
{
    ;
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
    print_matrix("(costruct) A = ", m, n, A, lda);

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
    print_matrix("(costruct) V = ", m, n, V, m);
    print_matrix("(costruct) T = ", n, n, T, ldt);

    // VT = V * Tを計算
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, n, 1.0, V, m, T, ldt, 0.0, VT, m);
    print_matrix("(costruct) VT = ", m, n, VT, m);

    // Q =I - VT * V^t を計算する
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, m, m, n, -1.0, VT, m, V, m, 0.0, Q, m);
    for (int i = 0; i < m; i++)
    {
        Q[i + m * i] += 1.0;
    }
    print_matrix("(costruct) Q = ", m, m, Q, m);

    free(V);
    free(VT);
    // printf("out construct Q\n");

    return Q;
}

// 半帯幅LのN*N帯行列Bから，帯行列化前の行列Aを復元する
// 復元には帯行列化で用いた直交行列のcompact-WY表現を格納した行列T，Yを用いる
void restore(int N, int L, double *B, int ldb, double *T, int ldt, double *Y, int ldy)
{
    int nb = L;
    int ldt_iter = nb;
    assert(MIN(N, L) >= nb && nb >= 1);
    assert(ldt_iter >= nb);
    assert(N % L == 0);
    double *T_iter; //= calloc(ldt_iter * L, sizeof(double));
    for (int k = N / L - 2; k >= 0; k--)
    {

        // print_matrix("A = ", N, N, B, ldb);

        printf("iteration %d/%d\n", k + 1, N / L - 1);
        int Nk = N - L - k * L;
        double *T_iter = &T[0 + L * k * ldt];
        double *V; // = malloc(Nk * L * sizeof(double));
        V = &Y[0 + (L * k) * ldy];
        // print_matrix("T iter = ", L, L, T_iter, ldt_iter);
        // print_matrix("V = ", Nk, L, V, ldy);

        double *update_tmp = malloc(Nk * Nk * sizeof(double));
        double *update_P = malloc(Nk * L * sizeof(double));
        double *update_beta = malloc(L * L * sizeof(double));
        double *update_Q = malloc(Nk * L * sizeof(double));
        //	);
        // // Aの(k+1,k)ブロックの先頭アドレス

        double *a_part = &B[(k + 1) * L + ldb * k * L];
        // print_matrix("A part = ", Nk, L, a_part, ldb);
        double *update_QR = malloc(Nk * L * sizeof(double));

        // QR=(I-VTV^T)a_partを計算する
        // update_QR=V^Ta_part
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, Nk, L, Nk, 1.0, V, ldy, a_part, ldb, 0.0, update_QR, Nk);
        // update_tmp=Tupdate_QR
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, L, L, L, 1.0, T_iter, ldt, update_QR, Nk, 0.0, update_tmp, Nk);
        // update_QR=Vupdate_tmp
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Nk, L, L, 1.0, V, ldy, update_tmp, Nk, 0.0, update_QR, Nk);
        // a_part=a_part-update_QR
        for (int i = 0; i < Nk; i++)
        {
            for (int j = 0; j < L; j++)
            {
                a_part[i + j * ldb] -= update_QR[i + j * Nk];
            }
        }
        for (int i = 0; i < N; i++)
        {
            for (int j = i; j < N; j++)
            {
                //  B[j + i * ldb] = B[i + j * ldb];
                B[i + j * ldb] = B[j + i * ldb];
            }
        }
        // print_matrix("A part = ", Nk, L, a_part, ldb);
        double *b_part = &B[(k + 1) * L + ldb * (k + 1) * L];

        for (int i = 0; i < Nk; i++)
            for (int j = 0; j < Nk; j++)
            {
                assert(b_part[i + j * ldb] == b_part[j + i * ldb]);
            }
        // 山本有作先生の『キャッシュマシン向け三重対角化アルゴリズムの性能予測方式』で説明されているBischofのアルゴリズム中のPを構築する
        cblas_dsymm(CblasColMajor, CblasLeft, CblasLower, Nk, L, 1.0, b_part, ldb, V, ldy, 0.0, update_tmp, Nk);
        // print_matrix("!b_part", Nk, Nk, b_part, ldb);
        // print_matrix("!V", Nk, L, V, ldy);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, Nk, L, L, 1.0, update_tmp, Nk, T_iter, ldt_iter, 0.0, update_P, Nk);

        // betaを構築する
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, L, L, Nk, 0.5, V, ldy, update_P, Nk, 0.0, update_tmp, Nk);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, L, L, L, 1.0, T_iter, ldt_iter, update_tmp, Nk, 0.0, update_beta, L);

        // Qを構築する
        for (int i = 0; i < Nk; i++)
        {
            for (int j = 0; j < L; j++)
            {
                update_Q[i + j * Nk] = update_P[i + j * Nk];
            }
        }
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Nk, L, L, -1.0, V, ldy, update_beta, L, 1.0, update_Q, Nk);
        //  print_matrix("update_Q", Nk, L, update_Q, Nk);
        // Aをrank-2k更新する
        cblas_dsyr2k(CblasColMajor, CblasLower, CblasNoTrans, Nk, L, -1.0, V, ldy, update_Q, Nk, 1.0, b_part, ldb);

        // measure_time(
        free(update_P);
        free(update_beta);
        free(update_Q);
        free(update_QR);
        free(update_tmp);
        // free(V););
    }
    for (int j = 0; j < N; j++)
    {
        for (int i = j; i < N; i++)

        {
            B[j + i * ldb] = B[i + j * ldb];
        }
    }
}

// (N,N)次元の実対称行列A(leading dimensionはlda)を帯行列化し、変換後の帯行列BとB=QAQ^tなるQを求める。関数終了時にAの要素はBの要素で上書きされる
void bischof(int N, int L, double *A, int lda, double *T, int ldt, double *Y, int ldy)
{
    int nb = L;
    int ldt_iter = nb;
    assert(MIN(N, L) >= nb && nb >= 1);
    assert(ldt_iter >= nb);
    assert(N % L == 0);
    double *const T_iter = calloc(ldt_iter * L, sizeof(double));

    for (int k = 0; k < N / L - 1; k++)
    {
        // print_matrix("A = ", N, N, A, lda);
        //  printf("iteration %d/%d\n", k + 1, N / L - 1);
        int Nk = N - L - k * L;
        // Aの(k+1,k)ブロック以下をQR分解する
        measure_time(LAPACKE_dgeqrt(LAPACK_COL_MAJOR, Nk, L, L, &A[(k + 1) * L + lda * k * L], lda, T_iter, ldt_iter));
        //        print_matrix("A=", N, N, A, lda);
        // print_matrix("T iter = ", L, L, T_iter, ldt_iter);
        // TにT_iterを代入
        // print_matrix("T iter = ", L, L, T_iter, ldt_iter);
        measure_time(
            for (int i = 0; i < L; i++) {
                for (int j = 0; j < L; j++)
                {
                    if (j >= i)
                        T[i + (j + L * k) * ldt] = T_iter[i + j * ldt_iter];
                    else
                        T[i + (j + L * k) * ldt] = T_iter[i + j * ldt_iter] = 0.0;
                }
            });
        //        print_matrix("T = ", L, N, T, ldt);

        // Aを更新する

        //        measure_time(
        double *V = malloc(Nk * L * sizeof(double));
        double *update_tmp = malloc(Nk * Nk * sizeof(double));
        double *update_P = malloc(Nk * L * sizeof(double));
        double *update_beta = malloc(L * L * sizeof(double));
        double *update_Q = malloc(Nk * L * sizeof(double));
        //	);
        // Aの(k+1,k)ブロックの先頭アドレス
        double *a_part = &A[(k + 1) * L + lda * k * L];
        // print_matrix("A part = ", Nk, L, a_part, lda);

        // Vにa_partの下三角部分を代入する(LAPACKE_dgeqrtを実行しているので，a_partの下三角部分にはQのcompact-WY表現の一部が入っている)
        measure_time(
            for (int i = 0; i < Nk; i++) {
                for (int j = 0; j < L; j++)
                {
                    if (i == j)
                    {
                        V[i + j * Nk] = 1.0;
                    }
                    else if (i > j)
                    {
                        V[i + j * Nk] = a_part[i + j * lda];
                        a_part[i + j * lda] = 0.0;
                    }
                    else
                    {
                        V[i + j * Nk] = 0.0;
                    }
                }
            });
        // print_matrix("V = ", Nk, L, V, Nk);

        // YにVを代入
        for (int j = 0; j < L; j++)
        {
            for (int i = 0; i < Nk; i++)
            {
                Y[i + (j + L * k) * ldy] = V[i + j * Nk];
            }
        }
        // print_matrix("Y=", Nk, L, Y, ldy);

        // a_partの下三角部分を0クリアし，上三角部分の要素の値を，その要素の位置と対称な位置へ代入する（Aは対称行列なので）
        measure_time(
            for (int i = L * k + L; i < N; i++) {
                for (int j = L * k; j < L * (k + 1); j++)
                {
                    A[j + i * lda] = A[i + j * lda];
                }
            });
        // 山本有作先生の『キャッシュマシン向け三重対角化アルゴリズムの性能予測方式』で説明されているBischofのアルゴリズム中のPを構築する
        measure_time(
            cblas_dsymm(CblasColMajor, CblasLeft, CblasLower, Nk, L, 1.0, &A[(k + 1) * L + lda * (k + 1) * L], lda, V, Nk, 0.0, update_tmp, Nk);

            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Nk, L, L, 1.0, update_tmp, Nk, T_iter, ldt_iter, 0.0, update_P, Nk););
        // print_matrix("update_P=", Nk, L, update_P, L);

        // betaを構築する
        measure_time(
            cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, L, L, Nk, 0.5, V, Nk, update_P, Nk, 0.0, update_tmp, Nk);
            cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, L, L, L, 1.0, T_iter, ldt_iter, update_tmp, Nk, 0.0, update_beta, L););
        //  print_matrix("update_beta=", L, L, update_beta, L);

        // Qを構築する
        measure_time(
            for (int i = 0; i < Nk; i++) {
                for (int j = 0; j < L; j++)
                {
                    update_Q[i + j * Nk] = update_P[i + j * Nk];
                }
            });
        // print_matrix("update_Q=", Nk, L, update_Q, Nk);

        measure_time(
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Nk, L, L, -1.0, V, Nk, update_beta, L, 1.0, update_Q, Nk););
        //  print_matrix("update_Q", Nk, L, update_Q, Nk);
        //  print_matrix("V=", Nk, L, V, Nk);
        //  print_matrix("update_Q=", Nk, L, update_Q, Nk);

        // Aをrank-2k更新する
        measure_time(
            cblas_dsyr2k(CblasColMajor, CblasLower, CblasNoTrans, Nk, L, -1.0, V, Nk, update_Q, Nk, 1.0, &A[(k + 1) * L + lda * (k + 1) * L], lda));

        measure_time(
            free(update_P);
            free(update_beta);
            free(update_Q);
            free(update_tmp);
            free(V););
    }
    free(T_iter);
    measure_time(
        for (int i = 0; i < N; i++) {
            for (int j = i; j < N; j++)
            {
                A[i + j * lda] = A[j + i * lda];
            }
        });
}

int main(int argc, char **argv)
{
    unsigned long const random_seed = 100;
    srand(random_seed);

    // size_t const m = 10000; // 帯行列化する行列のサイズ
    if (argc != 2)
    {
        printf("Usage %s matrix_size\n", argv[0]);
        return 0;
    }
    size_t const m = atoi(argv[1]); // 帯行列化する行列のサイズ
    size_t const lda = m + 1;
    int L = 2; //帯行列化時の幅

    double *const a = malloc(sizeof(double) * m * lda); // 帯行列化する行列
    double *const b = malloc(sizeof(double) * m * lda); // aのコピーを格納する行列
    // double *Q = malloc(m * m * sizeof(double));           // Bischofのアルゴリズムで計算されるQを格納する行列
    double *T = calloc(L * m, sizeof(double)); // Bischofのアルゴリズムで計算されるQのcompact-WY表現を格納する行列
    double *Y = calloc(m * m, sizeof(double)); // Bischofのアルゴリズムで計算されるQのcompact-WY表現を格納する行列

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

    //  print_matrix("A= ", m, m, a, lda);
    // aのノルムを計算しておく
    //    double norm_a = calc_Frobenius_norm(m, m, a, lda);
    // Bischofのアルゴリズム実行
    // double t1=omp_get_wtime();
    measure_time(bischof(m, L, a, lda, T, L, Y, m));
    // print_matrix("B= ", m, m, a, lda);
    //  printf("time : %lf [sec]\n",omp_get_wtime()-t1);
    //  print_matrix("L = ", L, m, T, L);
    //  print_matrix("Y = ", m, m, Y, m);
    /// restore(m, L, a, lda, T, L, Y, m);
    // print_matrix("rescored A= ", m, m, a, lda);

    // print_matrix("B = ", m, m, a, lda);

    // bはaのコピーなので，この順番が正しい

    // check(m, L, b, lda, a, lda, T, L);

    // Qが直交行列かQQ^Tを計算することで確かめる
    //  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, m, m, m, 1.0, Q, m, Q, m, 0.0, tmp, lda);

    // //直交行列ならこの値が1になっているはず
    // printf("norm of QQ^T = %e\n", calc_Frobenius_norm(m, m, tmp, lda) / sqrt(m));

    // //QaQ^Tを計算することで，帯行列化したaから元の行列aを復元する
    // cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, m, m, m, 1.0, Q, m, a, lda, 0.0, tmp, lda);
    // cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, m, m, 1.0, tmp, lda, Q, m, 0.0, a, lda);

    // double a_norm = calc_Frobenius_norm(m, m, b, lda);
    // // 復元したaと元のaとの誤差を計算する
    // for (size_t i = 0; i < m; ++i)
    // {
    //     for (size_t j = 0; j < m; ++j)
    //     {
    //         b[i + lda * j] -= a[i + j * lda];
    //     }
    // }

    // printf("norm  = %e\n", calc_Frobenius_norm(m, m, b, lda) / a_norm);

    // printf("END\n");

    free(a);
    free(b);
    free(tmp);
    free(T);
    free(Y);
    return 0;

    return EXIT_SUCCESS;
}
