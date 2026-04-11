#!/bin/bash
#SBATCH --account=agmisc
#SBATCH --job-name=su2_continue
#SBATCH --partition=fuchs
#SBATCH --qos=fuchs
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --time=48:00:00
#SBATCH --output=logs/continue_%j.out
#SBATCH --error=logs/continue_%j.err

cd /work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations
mkdir -p logs

export OMP_NUM_THREADS=1

./build/bin/mc_heatbath -i data/results/T30_L30_b2.7_periodic_seed1924684401_su2/input.txt &
./build/bin/mc_heatbath -i data/results/T30_L30_b2.7_periodic_seed2039470449_su2/input.txt &
./build/bin/mc_heatbath -i data/results/T38_L30_b2.7_open_seed6973043_su2/input.txt &
./build/bin/mc_heatbath -i data/results/T38_L30_b2.7_open_seed122058099_su2/input.txt &
wait
