#!/bin/bash
#SBATCH --account=mesonqcd
#SBATCH --job-name=su3_obc_hb
#SBATCH --partition=general2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=14-00:00:00
#SBATCH --output=/work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations/logs/su3_obc_hb_%j.out
#SBATCH --error=/work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations/logs/su3_obc_hb_%j.err

set -euo pipefail

# Generic pure-gauge SU(3) CL2QCD heatbath run with open temporal gauge
# boundaries. Override parameters with sbatch --export=ALL,NAME=value.

module purge

PROJECT_DIR=/work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations
CL2QCD_DIR=/home/mesonqcd/public/barros/cl2qcd_obc_barros
CL2QCD_BIN=${CL2QCD_BIN:-$CL2QCD_DIR/build/su3heatbath}
RESULTS_ROOT=${RESULTS_ROOT:-$PROJECT_DIR/data/results_su3_ensembles}
LOG_DIR=$PROJECT_DIR/logs

NSPACE=${NSPACE:-16}
NTIME=${NTIME:-16}
BETA=${BETA:-6.0}
START_CONDITION=${START_CONDITION:-cold}
USE_OPEN_TEMPORAL_GAUGE_BOUNDARY=${USE_OPEN_TEMPORAL_GAUGE_BOUNDARY:-true}

N_THERMALIZATION_STEPS=${N_THERMALIZATION_STEPS:-0}
N_HEATBATH_STEPS=${N_HEATBATH_STEPS:-1000}
N_OVERRELAXATION_STEPS=${N_OVERRELAXATION_STEPS:-1}
ONLINE_MEASURE_EVERY=${ONLINE_MEASURE_EVERY:-1}
CREATE_CHECKPOINT_EVERY=${CREATE_CHECKPOINT_EVERY:-100}
OVERWRITE_TEMPORARY_CHECKPOINT_EVERY=${OVERWRITE_TEMPORARY_CHECKPOINT_EVERY:-100}

PRECISION=${PRECISION:-64}
USE_GPU=${USE_GPU:-0}
USE_CPU=${USE_CPU:-1}
DEVICE_ID=${DEVICE_ID:-}
LOG_LEVEL=${LOG_LEVEL:-INFO}
N_DIGITS_IN_CONF_CHECKPOINT=${N_DIGITS_IN_CONF_CHECKPOINT:-5}

INITIAL_CONF=${INITIAL_CONF:-}
INITIAL_PRNG=${INITIAL_PRNG:-}

generate_host_seed() {
    python3 - "$NSPACE" "$NTIME" <<'PY'
import random
import sys

ns = int(sys.argv[1])
nt = int(sys.argv[2])
volume = ns * ns * ns * nt
max_seed = int(10_000_000_000 // volume)
if max_seed < 1:
    raise SystemExit("Lattice volume is too large for CL2QCD ranluxcl seed initialization")
print(random.randint(1, max_seed))
PY
}

HOST_SEED=${HOST_SEED:-$(generate_host_seed)}

BOUNDARY_LABEL=obc
if [ "$USE_OPEN_TEMPORAL_GAUGE_BOUNDARY" != "true" ] && [ "$USE_OPEN_TEMPORAL_GAUGE_BOUNDARY" != "1" ]; then
    BOUNDARY_LABEL=periodic
fi

RUN_LABEL=${RUN_LABEL:-T${NTIME}_L${NSPACE}_b${BETA}_${BOUNDARY_LABEL}_seed${HOST_SEED}_su3_cl2qcd}
RUN_DIR=$RESULTS_ROOT/$RUN_LABEL
INPUT_FILE=$RUN_DIR/input_su3heatbath.txt

mkdir -p "$LOG_DIR" "$RUN_DIR/configs" "$RUN_DIR/prng"

if [ ! -x "$CL2QCD_BIN" ]; then
    echo "CL2QCD su3heatbath binary is missing or not executable: $CL2QCD_BIN" >&2
    exit 1
fi

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-40}
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

cat > "$INPUT_FILE" <<EOF_INPUT
# Auto-generated CL2QCD SU(3) heatbath input
# Generated: $(date)

nSpace=$NSPACE
nTime=$NTIME
useGPU=$USE_GPU
useCPU=$USE_CPU
precision=$PRECISION
hostSeed=$HOST_SEED
startCondition=$START_CONDITION
beta=$BETA
useOpenTemporalGaugeBoundary=$USE_OPEN_TEMPORAL_GAUGE_BOUNDARY

nThermalizationSteps=$N_THERMALIZATION_STEPS
nHeatbathSteps=$N_HEATBATH_STEPS
nOverrelaxationSteps=$N_OVERRELAXATION_STEPS
onlineMeasureEvery=$ONLINE_MEASURE_EVERY
createCheckpointEvery=$CREATE_CHECKPOINT_EVERY
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
    echo "deviceId=$DEVICE_ID" >> "$INPUT_FILE"
fi

if [ -n "$INITIAL_CONF" ]; then
    echo "initialConf=$INITIAL_CONF" >> "$INPUT_FILE"
fi

if [ -n "$INITIAL_PRNG" ]; then
    echo "initialPRNG=$INITIAL_PRNG" >> "$INPUT_FILE"
fi

cat > "$RUN_DIR/run_info.txt" <<EOF_INFO
generated=$(date)
cl2qcd_dir=$CL2QCD_DIR
cl2qcd_bin=$CL2QCD_BIN
run_dir=$RUN_DIR
input_file=$INPUT_FILE
nSpace=$NSPACE
nTime=$NTIME
beta=$BETA
boundary=$BOUNDARY_LABEL
hostSeed=$HOST_SEED
startCondition=$START_CONDITION
nThermalizationSteps=$N_THERMALIZATION_STEPS
nHeatbathSteps=$N_HEATBATH_STEPS
nOverrelaxationSteps=$N_OVERRELAXATION_STEPS
onlineMeasureEvery=$ONLINE_MEASURE_EVERY
createCheckpointEvery=$CREATE_CHECKPOINT_EVERY
overwriteTemporaryCheckpointEvery=$OVERWRITE_TEMPORARY_CHECKPOINT_EVERY
useGPU=$USE_GPU
useCPU=$USE_CPU
deviceId=$DEVICE_ID
EOF_INFO

echo "Starting CL2QCD SU(3) heatbath"
echo "Run directory: $RUN_DIR"
echo "Input file: $INPUT_FILE"
echo "Binary: $CL2QCD_BIN"

cd "$RUN_DIR"
"$CL2QCD_BIN" "$INPUT_FILE" 2>&1 | tee heatbath.log
