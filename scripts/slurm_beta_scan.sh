#!/bin/bash
#SBATCH --job-name=su2_beta_scan
#SBATCH --output=slurm_logs/job_%A_%a.out
#SBATCH --error=slurm_logs/job_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=8G

# adjust --cpus-per-task, --time, --mem, --partition for your cluster

set -euo pipefail

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

PARAMS_FILE="${PROJECT_DIR}/input/base_params_su2.txt"
BETA_FILE="${PROJECT_DIR}/input/beta_scan_su2.txt"
SEEDS_FILE="${PROJECT_DIR}/input/seeds_su2.txt"

if [[ ! -f "$PARAMS_FILE" ]]; then echo "ERROR: Missing $PARAMS_FILE"; exit 1; fi
if [[ ! -f "$BETA_FILE" ]]; then echo "ERROR: Missing $BETA_FILE"; exit 1; fi
if [[ ! -f "$SEEDS_FILE" ]]; then echo "ERROR: Missing $SEEDS_FILE"; exit 1; fi

# read betas and seeds into arrays
mapfile -t BETAS < <(grep -v '^#' "$BETA_FILE" | grep -v '^$')
mapfile -t SEEDS < <(grep -v '^#' "$SEEDS_FILE" | grep -v '^$')

N_BETAS=${#BETAS[@]}
N_SEEDS=${#SEEDS[@]}

# map SLURM_ARRAY_TASK_ID to (beta_idx, seed_idx)
TASK_ID=${SLURM_ARRAY_TASK_ID:-0}
BETA_IDX=$((TASK_ID / N_SEEDS))
SEED_IDX=$((TASK_ID % N_SEEDS))

if [[ $BETA_IDX -ge $N_BETAS ]]; then
    echo "ERROR: TASK_ID=$TASK_ID out of range (${N_BETAS} betas x ${N_SEEDS} seeds)"
    exit 1
fi

BETA="${BETAS[$BETA_IDX]}"
SEED="${SEEDS[$SEED_IDX]}"

# read base params
T_SIZE=$(grep '^T ' "$PARAMS_FILE" | awk '{print $2}')
L_SIZE=$(grep '^L ' "$PARAMS_FILE" | awk '{print $2}')
START_TYPE=$(grep '^start_type' "$PARAMS_FILE" | awk '{print $2}')
BOUNDARY=$(grep '^boundary' "$PARAMS_FILE" | awk '{print $2}')
NUM_SWEEPS=$(grep '^num_sweeps' "$PARAMS_FILE" | awk '{print $2}')
SAVE_INTERVAL=$(grep '^save_interval' "$PARAMS_FILE" | awk '{print $2}')
SMEAR_STEPS=$(grep '^smear_steps' "$PARAMS_FILE" | awk '{print $2}')
SMEAR_INTERVAL=$(grep '^smear_interval' "$PARAMS_FILE" | awk '{print $2}')
SMEAR_ALPHA=$(grep '^smear_alpha' "$PARAMS_FILE" | awk '{print $2}')
EXCLUDE_SLICES=$(grep '^exclude_boundary_slices' "$PARAMS_FILE" | awk '{print $2}')
BOUNDARY=${BOUNDARY:-periodic}
EXCLUDE_SLICES=${EXCLUDE_SLICES:-0}

# verify critical params are set
for VAR_NAME in T_SIZE L_SIZE NUM_SWEEPS SAVE_INTERVAL SMEAR_STEPS SMEAR_ALPHA; do
    if [[ -z "${!VAR_NAME}" ]]; then
        echo "ERROR: Missing parameter $VAR_NAME in $PARAMS_FILE"
        exit 1
    fi
done

RUN_DIR="${PROJECT_DIR}/data/results/T${T_SIZE}_L${L_SIZE}_b${BETA}_${BOUNDARY}_seed${SEED}"
CONFIG_DIR="${RUN_DIR}/configs"
OUTPUT_DIR="${RUN_DIR}/output"

mkdir -p "$CONFIG_DIR" "$OUTPUT_DIR" "${PROJECT_DIR}/slurm_logs"

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $1" | tee -a "${RUN_DIR}/job_errors.log"
}

echo "=== Job $TASK_ID: beta=$BETA seed=$SEED ==="
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "Run dir: $RUN_DIR"

# generate heatbath input file
HB_INPUT="${RUN_DIR}/heatbath_input.txt"
cat > "$HB_INPUT" << EOF
T               ${T_SIZE}
L               ${L_SIZE}
beta            ${BETA}
seed            ${SEED}
start_type      ${START_TYPE:-cold}
boundary        ${BOUNDARY}
num_sweeps      ${NUM_SWEEPS}
save_interval   ${SAVE_INTERVAL}
output_dir      ${RUN_DIR}/
config_dir      ${CONFIG_DIR}/
EOF

# step 1: heatbath
echo "--- Running heatbath ---"
if ! "${PROJECT_DIR}/build/bin/mc_heatbath" -i "$HB_INPUT"; then
    log_error "mc_heatbath failed for beta=$BETA seed=$SEED"
    exit 1
fi

# step 2: detect thermalization
echo "--- Detecting thermalization ---"
PLAQ_FILE="${RUN_DIR}/plaquette.dat"
if [[ ! -f "$PLAQ_FILE" ]]; then
    log_error "plaquette.dat not found"
    exit 1
fi

START_CONF=$(python3 "${PROJECT_DIR}/scripts/detect_thermalization.py" \
    --plaquette-file "$PLAQ_FILE" \
    --save-interval "$SAVE_INTERVAL" \
    2>&1) || {
    log_error "detect_thermalization.py failed: $START_CONF"
    exit 1
}
echo "Thermalization detected at conf $START_CONF"

# step 3: measure topological charge
echo "--- Measuring topological charge ---"
MEAS_INPUT="${RUN_DIR}/meas_input.txt"
cat > "$MEAS_INPUT" << EOF
config_dir      ${CONFIG_DIR}/
output_dir      ${OUTPUT_DIR}/
T               ${T_SIZE}
L               ${L_SIZE}
beta            ${BETA}
start_conf      ${START_CONF}
end_conf        ${NUM_SWEEPS}
conf_step       ${SAVE_INTERVAL}
smear_steps     ${SMEAR_STEPS}
smear_interval  ${SMEAR_INTERVAL:-${SMEAR_STEPS}}
smear_alpha     ${SMEAR_ALPHA}
seed            ${SEED}
boundary        ${BOUNDARY}
exclude_boundary_slices ${EXCLUDE_SLICES}
EOF

if ! "${PROJECT_DIR}/build/bin/meas_topcharge" -i "$MEAS_INPUT"; then
    log_error "meas_topcharge failed for beta=$BETA seed=$SEED"
    exit 1
fi

echo "=== Job $TASK_ID complete ==="
