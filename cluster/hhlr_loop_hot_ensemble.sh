#!/bin/bash
#SBATCH --account=mesonqcd
#SBATCH --job-name=lqcd_hot_q
#SBATCH --partition=general2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=14-00:00:00
#SBATCH --output=/work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations/logs/hot_q_%j.out
#SBATCH --error=/work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations/logs/hot_q_%j.err

# SU(2) hot-start loop for a Q within tolerance around a target integer value,
# followed by config-saving heatbath production from that accepted field.
# Writes configs, plaquette.dat, and action_density.dat under HEATBATH_RESULTS_DIR.

module purge

# Output/results, binary, and script run from the work checkout.
export HEATBATH_RESULTS_DIR=/work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations/data/results
mkdir -p /work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations/logs

cd /work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

TARGET_Q_SPEC=${1:-}
Q_TOL=${2:-0.2}

ARGS=(
    --setup input/heatbath_topcharge_input/setup_80x80_periodic_su2.txt
    --beta 3.15
)
if [ -n "$TARGET_Q_SPEC" ]; then
    if [[ "$TARGET_Q_SPEC" == *,* ]]; then
        ARGS+=(--target-q-list "$TARGET_Q_SPEC")
    else
        ARGS+=(--target-q "$TARGET_Q_SPEC")
    fi
fi
ARGS+=(
    --q-tol "$Q_TOL"
    --candidate-sweeps 50
    --production-sweeps 10000
    --save-interval 100
    --max-tries 100
    --skip-build
)

bash scripts/loop_hot_ensemble.sh "${ARGS[@]}"
