#!/bin/bash
#SBATCH --account=agmisc
#SBATCH --job-name=su2_continue
#SBATCH --partition=fuchs
#SBATCH --qos=fuchs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=0
#SBATCH --time=48:00:00
#SBATCH --output=logs/continue_%j.out
#SBATCH --error=logs/continue_%j.err

cd /work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations
mkdir -p logs

export OMP_NUM_THREADS=20

#./build/bin/mc_heatbath -i data/results/T65_L65_b2.95_periodic_seed195149869_su2/input.txt &
./build/bin/mc_heatbath -i data/results/T65_L65_b2.95_periodic_seed2114821679_su2/input.txt &
wait
