#!/bin/bash
#SBATCH --account=mesonqcd
#SBATCH --job-name=lqcd_hot_q1
#SBATCH --partition=general2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=14-00:00:00
#SBATCH --output=/work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations/logs/hot_q1_%j.out
#SBATCH --error=/work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations/logs/hot_q1_%j.err

# SU(2) hot-start loop for a Q ~= 1 accepted initial configuration, followed by
# config-saving heatbath production from that accepted field. Writes configs,
# plaquette.dat, and action_density.dat under HEATBATH_RESULTS_DIR.

module purge

# Output/results go to work; binary and script run from home.
export HEATBATH_RESULTS_DIR=/work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations/data/results
mkdir -p /work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations/logs

cd /home/mesonqcd/barros/SU2_Calc/Topological-freezing-in-SU2-3-simulations

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

bash scripts/loop_hot_ensemble.sh \
    --setup input/heatbath_topcharge_input/setup_80x80_periodic_su2.txt \
    --beta 3.15 \
    --target-q 1 \
    --q-tol 0.2 \
    --candidate-sweeps 50 \
    --production-sweeps 10000 \
    --save-interval 100 \
    --max-tries 100 \
    --skip-build
