#!/bin/sh
#------ pjsub option --------#
#PJM -L rscgrp=ea
#PJM -L node=4
#PJM --mpi proc=8
#PJM -L elapse=1:00:00
#PJM -g n22247
#PJM -j
#------- Program execution -------#
#PJM -g s23365

# Note: #PJM -L rscgrp=ea
# Note: #PJM -L node=4
# Note: #PJM --mpi proc=8

# Note: #PJM -L rscgrp=ea
# Note: #PJM -L node=4
# Note: #PJM --mpi proc=4

# Note: #PJM -L rscgrp=n22247a
# Note: #PJM -L node=1
# Note: #PJM --mpi proc=2

# 20C/40T * 2 CPUs / Node
#export OMP_NUM_THREADS=20
export OMP_NUM_THREADS=40
#export OMP_NUM_THREADS=80

ulimit -s unlimited
date
cd $PJM_O_WORKDIR
set -x
for i in 1 2 4 8
do
    for N in 32768 16384 8192 4096 
    do
        mpiexec.hydra -n $i ../parallel $N  > out_bischof_`printf "%06d" "${N}"`_`printf "%04d" "${i}"`.txt
    done
done

 
