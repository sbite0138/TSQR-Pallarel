CC = gcc -std=c11 -Wall -Wextra -O3 -march=native -funroll-loops -fopenmp
OPENBLAS_DIR = ./OpenBLAS-0.3.13

all: TSQR_test 


TSQR_test: code.c
	$(CC) -isystem $(OPENBLAS_DIR)/lapack-netlib/LAPACKE/include $< -L$(OPENBLAS_DIR) -lopenblas -lm -o $@



