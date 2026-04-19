#!/bin/bash
#SBATCH --account=mesonqcd
#SBATCH --job-name=su2_continue
#SBATCH --partition=general2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=11-00:00:00
#SBATCH --output=/work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations/logs/continue_%j.out
#SBATCH --error=/work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations/logs/continue_%j.err

# Purge the environment to prevent dynamic linking conflicts
module purge

cd /work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations
mkdir -p logs

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#./build/bin/mc_heatbath -i data/results/T65_L65_b2.95_periodic_seed195149869_su2/input.txt &
./build/bin/mc_heatbath -i data/results/T65_L65_b2.95_periodic_seed2114821679_su2/input.txt &
wait
