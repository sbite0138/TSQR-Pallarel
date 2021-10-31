#!/bin/bash
#============ PBS Options ============
#QSUB -q gr10037b
#QSUB -ug gr10037
#QSUB -W 9:00
#QSUB -A p=10:t=4:c=2

#============ Shell Script ============
cd $QSUB_WORKDIR
set -x

mpiexec.hydra ../parallel 4000 4000 64 64 > out.txt
 
