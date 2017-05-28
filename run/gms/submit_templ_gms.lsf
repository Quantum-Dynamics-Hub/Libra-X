#!/bin/bash
#BSUB -R "span[hosts=1]"
#BSUB -n 1
#BSUB -J "L-X_H2O"

echo "LSB_JOBID="$LSB_JOBID
echo "LSB_QUEUE="$LSB_QUEUE
echo "LS_SUBCWD="$LS_SUBCWD
#echo "LSB_PROCS="$LSB_PROCS

/usr/bin/time python run_gms.py >> namd.out 2>&1
