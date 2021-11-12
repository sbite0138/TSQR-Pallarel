@ {"rank":0, "line":474, "cmd":"for (size_t i = 0; i < A->global_row; ++i) { for (size_t j = 0; j < A->global_col; ++j) { double r = (double)(rand()) / RAND_MAX; set(A, i, j, r); set(A, j, i, r); } }", "time":0.001063}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000020}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000046}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000080}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000126}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000108}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000010}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000098}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000083}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000594}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000033}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000019}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000082}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000395}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000120}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000003}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000363}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000116}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.001201}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000019}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000065}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000124}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000140}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000003}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000079}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000078}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000592}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000018}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000064}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000129}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000135}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000003}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000078}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000075}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000628}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000052}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000122}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000154}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000003}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000081}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000085}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000635}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000018}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000020}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000065}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000144}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000130}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000003}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000079}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000077}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000648}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000054}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000133}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000129}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000003}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000079}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000079}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000643}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.003301}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000016}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.003337}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000126}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000144}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000003}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000079}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000076}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.003913}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000017}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000058}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000107}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000187}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000003}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000079}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000077}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000615}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000051}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000110}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000089}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000082}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000081}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000532}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000016}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000054}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000084}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000078}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000080}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000088}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000503}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000013}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000051}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000083}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000078}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000078}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000076}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000486}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000014}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000017}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000053}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000162}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000081}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000079}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000078}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000567}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000014}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000053}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.007297}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000082}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000141}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000079}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.007940}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000014}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000016}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000086}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000191}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000083}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000208}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000083}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000936}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000014}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000083}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000213}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000087}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000079}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000185}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000918}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000046}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000023}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000125}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000089}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000081}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000190}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000080}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000838}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.004811}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000016}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.004885}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000197}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000086}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000211}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000080}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.005731}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000014}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000082}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000083}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000189}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000080}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000079}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000781}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000014}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000081}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000216}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000081}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.003286}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000132}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.004073}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000061}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000025}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000148}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000090}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000105}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000081}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000100}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.001143}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000014}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000276}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000084}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000080}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000081}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000080}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.001266}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000238}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000086}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.008147}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000148}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000098}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.009425}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000019}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000256}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000106}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000117}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000080}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000079}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.001363}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000013}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000016}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000256}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000089}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000103}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000003}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000086}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000081}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.005026}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000080}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000085}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000081}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000076}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000076}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000505}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000014}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000041}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000081}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000078}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000077}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000075}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000447}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000013}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000041}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000085}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000077}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000078}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000075}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000448}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000014}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000042}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000084}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000079}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000076}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000075}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000430}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000014}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000042}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000082}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000079}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000090}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000076}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000467}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000014}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000042}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000088}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000078}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000077}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000076}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000435}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000017}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000016}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000046}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000083}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000077}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000078}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000077}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000443}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000016}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000043}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000082}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000081}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000078}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000075}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000453}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000014}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000041}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000084}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000079}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000077}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000075}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000438}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000014}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000041}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000082}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000078}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000077}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000077}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000428}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000014}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000042}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000082}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000078}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000078}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000076}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000436}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000014}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000042}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000086}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000078}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000078}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000077}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000434}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000014}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000014}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000040}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000082}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000078}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000079}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000077}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000437}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000015}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000014}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000041}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000084}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000080}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000077}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000078}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000441}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000011}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000018}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000043}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000083}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000077}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000077}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000074}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000001}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000427}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000011}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000007}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000059}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000099}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000078}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000084}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000076}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000471}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000011}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000006}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000035}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000080}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000077}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000003}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000078}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000076}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000421}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000011}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000006}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000036}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000081}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000078}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000078}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000075}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000424}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000010}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000006}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000035}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000081}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000078}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000077}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000076}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000003}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000421}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000011}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000006}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000035}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000081}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000078}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000078}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000076}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000422}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000012}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000006}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000037}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000080}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000077}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000078}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000076}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000424}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000011}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000006}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000036}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000085}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000077}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000080}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000076}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000429}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000011}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000006}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000036}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000080}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000078}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000001}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000077}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000077}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000422}
@ {"rank":0, "line":208, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000011}
@ {"rank":0, "line":213, "cmd":"pdgeqrf_(&m, &n, matrix->data, &row, &col, matrix->desc, tau, work, &lwork, &info)", "time":0.000006}
@ {"rank":0, "line":229, "cmd":"pdgeqrf_wrap(m, n, matrix, row, col, tau)", "time":0.000033}
@ {"rank":0, "line":249, "cmd":"pdgemr2d_wrap(m - 1, 1, matrix, row + 1, col, Y, 1, 0)", "time":0.000081}
@ {"rank":0, "line":256, "cmd":"pdlaset_wrap('A', j, 1, 0.0, 0.0, y, 0, 0)", "time":0.000000}
@ {"rank":0, "line":257, "cmd":"set(y, j, 0, 1.0)", "time":0.000000}
@ {"rank":0, "line":258, "cmd":"pdgemr2d_wrap(m - j - 1, 1, matrix, row + j + 1, col + j, y, j + 1, 0)", "time":0.000000}
@ {"rank":0, "line":260, "cmd":"pdgemv_wrap('T', m, n, 1.0, Y, 0, 0, y, 0, 0, 1, 0.0, tmp, 0, 0, 1)", "time":0.000002}
@ {"rank":0, "line":261, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":263, "cmd":"pdgemv_wrap('N', j, j, -val, T, 0, 0, tmp, 0, 0, 1, 0.0, z, 0, 0, 1)", "time":0.000000}
@ {"rank":0, "line":268, "cmd":"pdgemr2d_wrap(m, 1, y, 0, 0, Y, 0, j)", "time":0.000077}
@ {"rank":0, "line":269, "cmd":"pdgemr2d_wrap(j, 1, z, 0, 0, T, 0, j)", "time":0.000076}
@ {"rank":0, "line":270, "cmd":"pdelget_(ADDR(char, 'A'), ADDR(char, ' '), &val, tau, ADDR(int, 1), ADDR(int, 1 + col + j), desc)", "time":0.000002}
@ {"rank":0, "line":272, "cmd":"set(T, j, j, val)", "time":0.000000}
@ {"rank":0, "line":296, "cmd":"pdgeqrt_wrap(rank, nproc_row, nproc_col, Nk, L, A, (k + 1) * L, k * L, T_iter)", "time":0.000338}
@ {"rank":0, "line":491, "cmd":"bischof(rank, nproc_row, nproc_col, m, L, A, T, Y)", "time":0.220885}
0.315598
