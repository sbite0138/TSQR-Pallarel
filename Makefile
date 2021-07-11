UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
CC = /opt/homebrew/bin/gcc-11  -std=c11 -Wall -Wextra -O3 -march=native -funroll-loops -fopenmp 
else
ifeq ($(UNAME), Linux)
CC = gcc -std=c11 -Wall -Wextra -g -march=native -funroll-loops -fopenmp -pg
endif
endif
OPENBLAS_DIR = ./OpenBLAS-0.3.13

all: TSQR_test bischof qrprec


TSQR_test: code.c
	$(CC) -isystem $(OPENBLAS_DIR)/lapack-netlib/LAPACKE/include $< -L$(OPENBLAS_DIR)/lib -lopenblas -lm -o $@

bischof: bischof.c
	$(CC) -isystem $(OPENBLAS_DIR)/lapack-netlib/LAPACKE/include $< -L$(OPENBLAS_DIR)/lib -lopenblas -lm -o $@

qrprec: qrprec.c
	$(CC) -isystem $(OPENBLAS_DIR)/lapack-netlib/LAPACKE/include $< -L$(OPENBLAS_DIR)/lib -lopenblas -lm -o $@

clean:
	rm -f qrprec
	rm -f bischof
	

