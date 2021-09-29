#!/bin/bash
#============ PBS Options ============
#QSUB -q gr10037b
#QSUB -ug gr10037
#QSUB -W 9:00
#QSUB -A p=1:t=72:c=36:m=64130M

#============ Shell Script ============
cd $QSUB_WORKDIR
set -x

../bischof 1000 > 1.txt
../bischof 2000 > 2.txt 
../bischof 4000 > 4.txt
../bischof 8000 > 8.txt
../bischof 16000 > 16.txt

