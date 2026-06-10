#!/bin/bash

if [ $# -ne 2 ]; then
    echo "Usage: ${0} <input_file <log_file>"
    exit
fi

#export OMP_NUM_THREADS=40
#export OMP_NUM_THREADS=8

tasks_per_step=1
cpus_per_task=40

PROJECT_DIR=${PROJECT_DIR:-/work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations}
taskscript="${PROJECT_DIR}/cluster/run_script_su3/compute_cl2qcd.sh"

input_file="${1}"
log_file="${2}"

date
echo -e "input:\t${input_file}"
echo -e "log:\t${log_file}"

srun -n${tasks_per_step} --cpus-per-task=${cpus_per_task} "${taskscript}" "${input_file}" "${log_file}"

date
