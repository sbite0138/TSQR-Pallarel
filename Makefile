UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
CC = /opt/homebrew/bin/gcc-11  -std=c11 -Wall -Wextra -O3 -march=native -funroll-loops -fopenmp 
else
ifeq ($(UNAME), Linux)
CC = gcc -std=c11 -Wall -Wextra -g -march=native -funroll-loops -fopenmp -pg -fsanitize=address
endif
endif
OPENBLAS_DIR = ./OpenBLAS-0.3.13

all: band_parallel

band_parallel: code.c
	mpicc  $< -lm -L'/home/sbite/github/scalapack-2.1.0/lib' -lscalapack -lopenblas -llapacke -llapack -lgfortran -g -DDEBUG -fopenmp  -std=c11  -o $@

clean:
	rm -f band_parallel
	

