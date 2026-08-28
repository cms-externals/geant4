#!/usr/bin/env bash
#

node_name=`uname -n`

echo " node name = ${node_name} "

echo " SLURM_SUBMIT_DIR = ${SLURM_SUBMIT_DIR}"

srun -l $1 $2


