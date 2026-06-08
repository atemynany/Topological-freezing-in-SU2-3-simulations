#!/bin/bash
#SBATCH --account=mesonqcd
#SBATCH --job-name=su3_E_obc_b680
#SBATCH --partition=general2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=14-00:00:00
#SBATCH --array=0-3
#SBATCH --output=/work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations/logs/su3_E_obc_b680_%A_%a.out
#SBATCH --error=/work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations/logs/su3_E_obc_b680_%A_%a.err

set -euo pipefail

# Pure-gauge SU(3) CL2QCD heatbath setup for ensemble E(OBC):
#   beta=6.80, a=0.0307 fm, lattice=40^3 x 140
#   N_sim=4, N_total=80000, N_or=15, N_therm=40000,
#   N_sep=100, N_meas=400, exclude=30 time slices per side.
#
# Important: this wrapper deliberately does not use CL2QCD's
# nThermalizationSteps loop. That loop measures gauge observables every
# thermalization step and is very slow on the CPU OpenCL backend.
#
# Instead, each stream runs two CL2QCD generation phases:
#   1. warmup:     normal heatbath/overrelax sweeps, no measurements,
#                  one checkpoint at the warmup boundary;
#   2. production: continue from that checkpoint, save/measure every N_sep.
#
# Submit with:
#   sbatch cluster/hhlr_su3_cl2qcd_heatbath_obc.sh
#
# Override parameters with:
#   sbatch --export=ALL,NAME=value cluster/hhlr_su3_cl2qcd_heatbath_obc.sh

module purge

PROJECT_DIR=/work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations
CL2QCD_DIR=/home/mesonqcd/public/barros/cl2qcd_obc_barros
CL2QCD_BIN=${CL2QCD_BIN:-$CL2QCD_DIR/build/su3heatbath}
RESULTS_ROOT=${RESULTS_ROOT:-$PROJECT_DIR/data/results_su3_ensembles}
LOG_DIR=$PROJECT_DIR/logs

ENSEMBLE_LABEL=${ENSEMBLE_LABEL:-E_OBC_b6p80_L40_T140}
N_SIM=${N_SIM:-4}
STREAM_ID=${SLURM_ARRAY_TASK_ID:-${STREAM_ID:-0}}

NSPACE=${NSPACE:-40}
NTIME=${NTIME:-140}
BETA=${BETA:-6.80}
LATTICE_SPACING_FM=${LATTICE_SPACING_FM:-0.0307}
START_CONDITION=${START_CONDITION:-hot}
USE_OPEN_TEMPORAL_GAUGE_BOUNDARY=${USE_OPEN_TEMPORAL_GAUGE_BOUNDARY:-true}
EXCLUDE_BOUNDARY_SLICES=${EXCLUDE_BOUNDARY_SLICES:-30}

N_TOTAL=${N_TOTAL:-80000}
N_THERM_TOTAL=${N_THERM_TOTAL:-40000}
N_SEP=${N_SEP:-100}
N_MEAS=${N_MEAS:-400}
N_OVERRELAXATION_STEPS=${N_OVERRELAXATION_STEPS:-15}

PRECISION=${PRECISION:-64}
USE_GPU=${USE_GPU:-0}
USE_CPU=${USE_CPU:-1}
DEVICE_ID=${DEVICE_ID:-}
LOG_LEVEL=${LOG_LEVEL:-INFO}
N_DIGITS_IN_CONF_CHECKPOINT=${N_DIGITS_IN_CONF_CHECKPOINT:-5}
DRY_RUN=${DRY_RUN:-false}
ALLOW_ODD_LATTICE=${ALLOW_ODD_LATTICE:-false}

INITIAL_CONF=${INITIAL_CONF:-}
INITIAL_PRNG=${INITIAL_PRNG:-}
HOST_SEED_BASE=${HOST_SEED_BASE:-137}
HOST_SEED_STRIDE=${HOST_SEED_STRIDE:-97}

if [ "$ALLOW_ODD_LATTICE" != "true" ] && [ "$ALLOW_ODD_LATTICE" != "1" ]; then
    if [ $((NSPACE % 2)) -ne 0 ] || [ $((NTIME % 2)) -ne 0 ]; then
        echo "Refusing to run CL2QCD with odd local lattice extents: L=$NSPACE, T=$NTIME" >&2
        echo "This CL2QCD build warns that odd local lattice sizes are known to cause problems." >&2
        echo "The previous L=39 production attempt segfaulted before writing configurations." >&2
        echo "Use an even lattice extent or a code path that supports odd lattice sizes." >&2
        exit 1
    fi
fi

if ! [[ "$STREAM_ID" =~ ^[0-9]+$ ]]; then
    echo "STREAM_ID must be a non-negative integer, got: $STREAM_ID" >&2
    exit 1
fi
if [ "$STREAM_ID" -ge "$N_SIM" ]; then
    echo "STREAM_ID=$STREAM_ID is outside N_SIM=$N_SIM" >&2
    exit 1
fi
if [ $((N_TOTAL % N_SIM)) -ne 0 ]; then
    echo "N_TOTAL=$N_TOTAL must be divisible by N_SIM=$N_SIM" >&2
    exit 1
fi
if [ $((N_THERM_TOTAL % N_SIM)) -ne 0 ]; then
    echo "N_THERM_TOTAL=$N_THERM_TOTAL must be divisible by N_SIM=$N_SIM" >&2
    exit 1
fi
if ! [[ "$EXCLUDE_BOUNDARY_SLICES" =~ ^[0-9]+$ ]]; then
    echo "EXCLUDE_BOUNDARY_SLICES must be a non-negative integer, got: $EXCLUDE_BOUNDARY_SLICES" >&2
    exit 1
fi
if [ "$EXCLUDE_BOUNDARY_SLICES" -lt 1 ]; then
    echo "EXCLUDE_BOUNDARY_SLICES must be positive for open boundaries" >&2
    exit 1
fi
if [ $((2 * EXCLUDE_BOUNDARY_SLICES)) -ge "$NTIME" ]; then
    echo "EXCLUDE_BOUNDARY_SLICES=$EXCLUDE_BOUNDARY_SLICES removes the full temporal extent T=$NTIME" >&2
    exit 1
fi

N_TOTAL_PER_STREAM=$((N_TOTAL / N_SIM))
N_WARMUP_STEPS=${N_WARMUP_STEPS:-$((N_THERM_TOTAL / N_SIM))}
N_PRODUCTION_STEPS=${N_HEATBATH_STEPS:-$((N_TOTAL_PER_STREAM - N_WARMUP_STEPS))}

ONLINE_MEASURE_EVERY=${ONLINE_MEASURE_EVERY:-$N_SEP}
CREATE_CHECKPOINT_EVERY=${CREATE_CHECKPOINT_EVERY:-$N_SEP}
OVERWRITE_TEMPORARY_CHECKPOINT_EVERY=${OVERWRITE_TEMPORARY_CHECKPOINT_EVERY:-0}

if [ "$N_WARMUP_STEPS" -lt 0 ]; then
    echo "N_WARMUP_STEPS=$N_WARMUP_STEPS must be non-negative" >&2
    exit 1
fi
if [ "$N_PRODUCTION_STEPS" -le 0 ]; then
    echo "N_PRODUCTION_STEPS=$N_PRODUCTION_STEPS must be positive" >&2
    exit 1
fi
if [ $((N_WARMUP_STEPS + N_PRODUCTION_STEPS)) -ne "$N_TOTAL_PER_STREAM" ]; then
    echo "Warmup + production mismatch: $N_WARMUP_STEPS + $N_PRODUCTION_STEPS != $N_TOTAL_PER_STREAM" >&2
    exit 1
fi
if [ $((N_PRODUCTION_STEPS % N_SEP)) -ne 0 ]; then
    echo "N_PRODUCTION_STEPS=$N_PRODUCTION_STEPS must be divisible by N_SEP=$N_SEP" >&2
    exit 1
fi
if [ "$N_WARMUP_STEPS" -gt 0 ] && [ $((N_WARMUP_STEPS % N_SEP)) -ne 0 ]; then
    echo "N_WARMUP_STEPS=$N_WARMUP_STEPS must be divisible by N_SEP=$N_SEP for clean production numbering" >&2
    exit 1
fi

N_MEAS_PER_STREAM=$((N_PRODUCTION_STEPS / N_SEP))
N_MEAS_TOTAL_COMPUTED=$((N_MEAS_PER_STREAM * N_SIM))
if [ "$N_MEAS_TOTAL_COMPUTED" -ne "$N_MEAS" ]; then
    echo "Measurement mismatch: computed $N_MEAS_TOTAL_COMPUTED but N_MEAS=$N_MEAS" >&2
    exit 1
fi

VOLUME=$((NSPACE * NSPACE * NSPACE * NTIME))
HOST_SEED_MAX=$((10000000000 / VOLUME))
if [ "$HOST_SEED_MAX" -lt 1 ]; then
    echo "Lattice volume $VOLUME is too large for CL2QCD ranluxcl seed initialization" >&2
    exit 1
fi

HOST_SEED=${HOST_SEED:-$((HOST_SEED_BASE + STREAM_ID * HOST_SEED_STRIDE))}
if [ "$HOST_SEED" -lt 1 ] || [ "$HOST_SEED" -gt "$HOST_SEED_MAX" ]; then
    echo "HOST_SEED=$HOST_SEED is outside CL2QCD allowed range 1..$HOST_SEED_MAX for volume $VOLUME" >&2
    exit 1
fi

BOUNDARY_LABEL=open
if [ "$USE_OPEN_TEMPORAL_GAUGE_BOUNDARY" != "true" ] && [ "$USE_OPEN_TEMPORAL_GAUGE_BOUNDARY" != "1" ]; then
    BOUNDARY_LABEL=periodic
fi
if [ "$BOUNDARY_LABEL" != "open" ]; then
    echo "This OBC setup requires USE_OPEN_TEMPORAL_GAUGE_BOUNDARY=true" >&2
    exit 1
fi

RUN_LABEL=${RUN_LABEL:-T${NTIME}_L${NSPACE}_b${BETA}_${BOUNDARY_LABEL}_seed${HOST_SEED}_stream${STREAM_ID}_su3_cl2qcd}
RUN_DIR=$RESULTS_ROOT/$RUN_LABEL
WARMUP_INPUT_FILE=$RUN_DIR/input_su3heatbath_warmup.txt
PRODUCTION_INPUT_FILE=$RUN_DIR/input_su3heatbath_production.txt
INPUT_FILE=$PRODUCTION_INPUT_FILE

mkdir -p "$LOG_DIR" "$RUN_DIR/configs" "$RUN_DIR/prng"

if [ ! -x "$CL2QCD_BIN" ]; then
    echo "CL2QCD su3heatbath binary is missing or not executable: $CL2QCD_BIN" >&2
    exit 1
fi

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-40}
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

format_checkpoint_number() {
    printf "%0${N_DIGITS_IN_CONF_CHECKPOINT}d" "$1"
}

write_input_file() {
    local input_file=$1
    local phase=$2
    local start_condition=$3
    local n_heatbath_steps=$4
    local online_measure_every=$5
    local create_checkpoint_every=$6
    local initial_conf=${7:-}
    local initial_prng=${8:-}

    cat > "$input_file" <<EOF_INPUT
# Auto-generated CL2QCD SU(3) heatbath input
# Generated: $(date)
# Ensemble: $ENSEMBLE_LABEL
# Stream: $STREAM_ID / $N_SIM
# Phase: $phase
# Topcharge boundary exclusion per side: $EXCLUDE_BOUNDARY_SLICES
# CL2QCD thermalization is disabled intentionally; warmup uses normal heatbath generation.

nSpace=$NSPACE
nTime=$NTIME
useGPU=$USE_GPU
useCPU=$USE_CPU
precision=$PRECISION
hostSeed=$HOST_SEED
startCondition=$start_condition
beta=$BETA
useOpenTemporalGaugeBoundary=$USE_OPEN_TEMPORAL_GAUGE_BOUNDARY

nThermalizationSteps=0
nHeatbathSteps=$n_heatbath_steps
nOverrelaxationSteps=$N_OVERRELAXATION_STEPS
onlineMeasureEvery=$online_measure_every
createCheckpointEvery=$create_checkpoint_every
overwriteTemporaryCheckpointEvery=$OVERWRITE_TEMPORARY_CHECKPOINT_EVERY
nDigitsInConfCheckpoint=$N_DIGITS_IN_CONF_CHECKPOINT

confPrefix=configs/conf.
confPostfix=
PRNGPrefix=prng/prng.
PRNGPostfix=
gaugeObsInSingleFile=true
gaugeObsPrefix=gaugeObs
gaugeObsPostfix=.dat
logLevel=$LOG_LEVEL
EOF_INPUT

    if [ -n "$DEVICE_ID" ]; then
        echo "deviceId=$DEVICE_ID" >> "$input_file"
    fi
    if [ -n "$initial_conf" ]; then
        echo "initialConf=$initial_conf" >> "$input_file"
    fi
    if [ -n "$initial_prng" ]; then
        echo "initialPRNG=$initial_prng" >> "$input_file"
    fi
}

WARMUP_CONF=
WARMUP_PRNG=
if [ "$N_WARMUP_STEPS" -gt 0 ]; then
    WARMUP_CHECKPOINT_NUMBER=$(format_checkpoint_number "$N_WARMUP_STEPS")
    WARMUP_CONF=configs/conf.$WARMUP_CHECKPOINT_NUMBER
    WARMUP_PRNG=prng/prng.$WARMUP_CHECKPOINT_NUMBER
    WARMUP_ONLINE_MEASURE_EVERY=$((N_WARMUP_STEPS + 1))

    write_input_file \
        "$WARMUP_INPUT_FILE" \
        "warmup_without_cl2qcd_thermalization" \
        "$START_CONDITION" \
        "$N_WARMUP_STEPS" \
        "$WARMUP_ONLINE_MEASURE_EVERY" \
        "$N_WARMUP_STEPS" \
        "$INITIAL_CONF" \
        "$INITIAL_PRNG"

    write_input_file \
        "$PRODUCTION_INPUT_FILE" \
        "production_continue_after_warmup" \
        "continue" \
        "$N_PRODUCTION_STEPS" \
        "$ONLINE_MEASURE_EVERY" \
        "$CREATE_CHECKPOINT_EVERY" \
        "$WARMUP_CONF" \
        "$WARMUP_PRNG"
else
    write_input_file \
        "$PRODUCTION_INPUT_FILE" \
        "production_no_warmup" \
        "$START_CONDITION" \
        "$N_PRODUCTION_STEPS" \
        "$ONLINE_MEASURE_EVERY" \
        "$CREATE_CHECKPOINT_EVERY" \
        "$INITIAL_CONF" \
        "$INITIAL_PRNG"
fi

cat > "$RUN_DIR/run_info.txt" <<EOF_INFO
generated=$(date)
ensemble_label=$ENSEMBLE_LABEL
cl2qcd_dir=$CL2QCD_DIR
cl2qcd_bin=$CL2QCD_BIN
run_dir=$RUN_DIR
input_file=$INPUT_FILE
warmup_input_file=$WARMUP_INPUT_FILE
production_input_file=$PRODUCTION_INPUT_FILE
stream_id=$STREAM_ID
nSim=$N_SIM
nTotal=$N_TOTAL
nTotalPerStream=$N_TOTAL_PER_STREAM
nThermTotal=$N_THERM_TOTAL
nSep=$N_SEP
nMeas=$N_MEAS
nMeasPerStream=$N_MEAS_PER_STREAM
latticeSpacingFm=$LATTICE_SPACING_FM
nSpace=$NSPACE
nTime=$NTIME
beta=$BETA
boundary=$BOUNDARY_LABEL
excludeBoundarySlices=$EXCLUDE_BOUNDARY_SLICES
hostSeed=$HOST_SEED
startCondition=$START_CONDITION
thermalizationMode=normal_heatbath_warmup_no_cl2qcd_thermalize
nThermalizationStepsCL2QCD=0
nWarmupHeatbathSteps=$N_WARMUP_STEPS
nProductionHeatbathSteps=$N_PRODUCTION_STEPS
nOverrelaxationSteps=$N_OVERRELAXATION_STEPS
onlineMeasureEveryProduction=$ONLINE_MEASURE_EVERY
createCheckpointEveryProduction=$CREATE_CHECKPOINT_EVERY
overwriteTemporaryCheckpointEvery=$OVERWRITE_TEMPORARY_CHECKPOINT_EVERY
warmupCheckpointConf=$WARMUP_CONF
warmupCheckpointPrng=$WARMUP_PRNG
useGPU=$USE_GPU
useCPU=$USE_CPU
deviceId=$DEVICE_ID
EOF_INFO

echo "Starting CL2QCD SU(3) heatbath"
echo "Run directory: $RUN_DIR"
echo "Warmup input file: $WARMUP_INPUT_FILE"
echo "Production input file: $PRODUCTION_INPUT_FILE"
echo "Binary: $CL2QCD_BIN"
echo "Boundary: $BOUNDARY_LABEL"
echo "Exclude boundary slices per side: $EXCLUDE_BOUNDARY_SLICES"
echo "Ensemble: $ENSEMBLE_LABEL"
echo "Stream: $STREAM_ID / $N_SIM"
echo "Per-stream total sweeps: $N_TOTAL_PER_STREAM"
echo "Warmup heatbath sweeps: $N_WARMUP_STEPS"
echo "CL2QCD thermalization sweeps: 0"
echo "Production heatbath sweeps: $N_PRODUCTION_STEPS"
echo "Expected saved measurements: $N_MEAS_PER_STREAM"

if [ "$DRY_RUN" = "true" ] || [ "$DRY_RUN" = "1" ]; then
    echo "DRY_RUN=$DRY_RUN, so CL2QCD was not started."
    exit 0
fi

cd "$RUN_DIR"
if [ "$N_WARMUP_STEPS" -gt 0 ]; then
    echo "Starting warmup phase without CL2QCD thermalization"
    "$CL2QCD_BIN" "$WARMUP_INPUT_FILE" 2>&1 | tee heatbath_warmup.log
    if [ ! -f "$WARMUP_CONF" ] || [ ! -f "$WARMUP_PRNG" ]; then
        echo "Warmup checkpoint missing: $WARMUP_CONF / $WARMUP_PRNG" >&2
        exit 1
    fi
fi

echo "Starting production phase"
"$CL2QCD_BIN" "$PRODUCTION_INPUT_FILE" 2>&1 | tee heatbath_production.log
