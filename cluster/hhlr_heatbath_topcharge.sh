#!/bin/bash
#SBATCH --account=mesonqcd
#SBATCH --job-name=lqcd_hb_tc
#SBATCH --partition=general1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=14-00:00:00
#SBATCH --output=/home/mesonqcd/barros/SU2_Calc/Topological-freezing-in-SU2-3-simulations/logs/hb_tc_%j.out
#SBATCH --error=/home/mesonqcd/barros/SU2_Calc/Topological-freezing-in-SU2-3-simulations/logs/hb_tc_%j.err
#SBATCH --open-mode=append

# Fused SU(2) heatbath + in-memory topological charge scan (HHLR).
# Working directory and results live in $HOME — no gauge configs are written;
# per ensemble only plaquette.dat, action_density.dat, topcharge.dat,
# topcharge_timeslice.dat.
# Uses the same thread layout as hhlr_heatbath.sh because smearing + topcharge
# are also OMP-parallel.

module purge

export HEATBATH_RESULTS_DIR=/home/mesonqcd/barros/SU2_Calc/Topological-freezing-in-SU2-3-simulations/data/results
mkdir -p /home/mesonqcd/barros/SU2_Calc/Topological-freezing-in-SU2-3-simulations/logs

cd /home/mesonqcd/barros/SU2_Calc/Topological-freezing-in-SU2-3-simulations

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

bash run_heatbath_topcharge_scan.sh --skip-build --jobs 1 --exclude "_later|_finished"
# --exclude uses bash ERE (plain | for alternation). "_later" skips the larger lattices
# (41^4 and up); "_finished" is a defensive guard for setup files flagged as done.
