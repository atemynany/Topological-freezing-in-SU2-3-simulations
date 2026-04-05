#!/bin/bash
#SBATCH --account=agmisc
#SBATCH --job-name=lqcd_scan
#SBATCH --partition=fuchs
#SBATCH --qos=fuchs
#SBATCH --nodes=1
#SBATCH --ntasks=20           
#SBATCH --cpus-per-task=1
#SBATCH --mem=0               
#SBATCH --time=48:00:00
#SBATCH --output=scan_%j.out
#SBATCH --error=scan_%j.err

cd /work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations

bash run_heatbath_scan.sh --su2 --skip-build --jobs 20