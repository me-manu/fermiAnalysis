#!/usr/bin/env bash

source $HOME/.bashrc

conda activate {0[conda_env]:s}

# run
echo "Running $SLURM_JOB_NAME with job id $SLURM_JOB_ID"
date

echo "calling python: {0[script]:s} -c {0[config]:s} {0[add_opt]:s} 1> {0[logdir]:s}/out.$SLURM_ARRAY_TASK_ID 2> {0[logdir]:s}/err.$SLURM_ARRAY_TASK_ID" 1>&2
python {0[script]:s} -c {0[config]:s} {0[add_opt]:s} 1> {0[logdir]:s}/out.$SLURM_ARRAY_TASK_ID 2> {0[logdir]:s}/err.$SLURM_ARRAY_TASK_ID 
if [ $? -eq 0 ]; then
    echo "python OK" 1>&2

