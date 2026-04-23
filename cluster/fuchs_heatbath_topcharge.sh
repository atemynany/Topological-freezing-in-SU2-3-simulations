#!/bin/bash
#SBATCH --account=agmisc
#SBATCH --job-name=lqcd_hb_tc
#SBATCH --partition=fuchs
#SBATCH --qos=fuchs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=0
#SBATCH --time=240:00:00
#SBATCH --output=/home/mesonqcd/barros/SU2_Calc/Topological-freezing-in-SU2-3-simulations/logs/hb_tc_%j.out
#SBATCH --error=/home/mesonqcd/barros/SU2_Calc/Topological-freezing-in-SU2-3-simulations/logs/hb_tc_%j.err

# Fused SU(2) heatbath + in-memory topological charge scan.
# Working directory and results live in $HOME (only plaquette.dat / topcharge.dat /
# topcharge_timeslice.dat per ensemble — no configs on disk). Matches the threading
# of fuchs_heatbath.sh (one multi-threaded job) because smearing + topcharge are
# also OMP-parallel.

export HEATBATH_RESULTS_DIR=/home/mesonqcd/barros/SU2_Calc/Topological-freezing-in-SU2-3-simulations/data/results
mkdir -p /home/mesonqcd/barros/SU2_Calc/Topological-freezing-in-SU2-3-simulations/logs

cd /home/mesonqcd/barros/SU2_Calc/Topological-freezing-in-SU2-3-simulations

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

bash run_heatbath_topcharge_scan.sh --skip-build --jobs 1
# Discovers all input/heatbath_topcharge_input/setup_*_su2.txt files automatically.
# Use --exclude "T65\|T81" to skip specific lattices, same pattern as fuchs_topcharge.sh.
