UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
CC = /opt/homebrew/bin/gcc-11  -std=c11 -Wall -Wextra -O3 -march=native -funroll-loops -fopenmp 
else
ifeq ($(UNAME), Linux)
CC = gcc -std=c11 -Wall -Wextra -g -march=native -funroll-loops -fopenmp -pg -fsanitize=address
endif
endif
OPENBLAS_DIR = ./OpenBLAS-0.3.13

all: tsqr_test bischof qrprec parallel


tsqr_test: tsqr_test.c
	mpicc  $< -lm -L'/home/sbite/github/scalapack-2.1.0/lib' -lscalapack -llapacke -llapack -lgfortran -g -fopenmp  -std=c11  -o $@

bischof: bischof.c
	$(CC) -isystem $(OPENBLAS_DIR)/lapack-netlib/LAPACKE/include $< -L$(OPENBLAS_DIR)/lib -lopenblas -lm -o $@

qrprec: qrprec.c
	$(CC) -isystem $(OPENBLAS_DIR)/lapack-netlib/LAPACKE/include $< -L$(OPENBLAS_DIR)/lib -lopenblas -lm -o $@

parallel: bischof_parallel.c
	 mpicc bischof_parallel.c -L'/home/sbite/github/scalapack-2.1.0/lib' -lscalapack -llapack -lgfortran  -fopenmp  -g -std=c11  -o $@


clean:
	rm -f qrprec
	rm -f bischof
	

