#!/bin/bash
#============ PBS Options ============
#QSUB -q gr10037b
#QSUB -ug gr10037
#QSUB -W 9:00
#QSUB -A p=8:t=36:c=18:m=50G

#============ Shell Script ============
ulimit -s unlimited
date
cd $QSUB_WORKDIR
set -x
for i in 1 2 3 4 5 6 7 8
do
    for N in 32768 16384 8192 4096 
    do
        mpiexec.hydra -n $i ../parallel $N  > out_bischof_`printf "%06d" "${N}"`_`printf "%04d" "${i}"`.txt
    done
done

 
