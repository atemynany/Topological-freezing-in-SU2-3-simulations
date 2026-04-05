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
#SBATCH --output=heatbath_%j.out  
#SBATCH --error=heatbath_%j.err    

# Navigate to your working directory
cd /work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations

# Execute the program
srun --exclusive -n 1 -c 1 ./build/bin/mc_heatbath -i input/run_input.txt
