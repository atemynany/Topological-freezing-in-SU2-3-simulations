#!/bin/bash
#SBATCH --account=agmisc
#SBATCH --job-name=su2_continue
#SBATCH --partition=fuchs
#SBATCH --qos=fuchs
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=10
#SBATCH --mem=0
#SBATCH --time=96:00:00
#SBATCH --output=logs/continue_%j.out
#SBATCH --error=logs/continue_%j.err

cd /work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations
mkdir -p logs

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=close
export OMP_PLACES=cores

srun -n1 --exclusive --cpus-per-task=$SLURM_CPUS_PER_TASK \
    ./build/bin/mc_heatbath -i data/results/T42_L30_b2.7_open_seed2082463229_su2/input.txt &
srun -n1 --exclusive --cpus-per-task=$SLURM_CPUS_PER_TASK \
    ./build/bin/mc_heatbath -i data/results/T42_L30_b2.7_open_seed49907711_su2/input.txt &
wait
