#!/bin/bash
mkdir -p run
cd run
mpiicc -std=c99 -o ../band_parallel  -O3 -xCORE-AVX512 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -mkl  ../code.c