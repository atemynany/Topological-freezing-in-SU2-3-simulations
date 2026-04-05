#!/bin/bash
#SBATCH --account=agmisc
#SBATCH --job-name=lqcd_heatbath
#SBATCH --partition=fuchs
#SBATCH --qos=fuchs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=logs/heatbath_%j.out
#SBATCH --error=logs/heatbath_%j.err

cd /work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations
mkdir -p logs

./build/bin/mc_heatbath -i input/run_input.txt
