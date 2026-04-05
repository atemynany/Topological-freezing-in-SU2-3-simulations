#!/bin/bash
#SBATCH --account=agmisc
#SBATCH --job-name=lqcd_topcharge
#SBATCH --partition=fuchs
#SBATCH --qos=fuchs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=topcharge_%j.out
#SBATCH --error=topcharge_%j.err

cd /work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations

srun --exclusive -n 1 -c 1 bash run_topcharge_scan.sh --su2 --skip-build