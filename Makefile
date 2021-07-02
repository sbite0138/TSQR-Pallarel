CC = gcc -std=c11 -Wall -Wextra -g -march=native -funroll-loops -fopenmp
OPENBLAS_DIR = ./OpenBLAS-0.3.13

all: TSQR_test bischof


TSQR_test: code.c
	$(CC) -isystem $(OPENBLAS_DIR)/lapack-netlib/LAPACKE/include $< -L$(OPENBLAS_DIR)/lib -lopenblas -lm -o $@

bischof: bischof.c
	$(CC) -isystem $(OPENBLAS_DIR)/lapack-netlib/LAPACKE/include $< -L$(OPENBLAS_DIR)/lib -lopenblas -lm -o $@



