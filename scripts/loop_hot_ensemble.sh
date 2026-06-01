#!/bin/bash
# ==============================================================================
# loop_hot_ensemble.sh
#
# Build one SU(2) config-writing hot ensemble after finding a hot-start
# candidate whose measured topological charge has Q ~= target_q.
#
# The accepted hot trial is copied into the production config directory as
# conf.0000, then mc_heatbath is run with start_type=continue.
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="$PROJECT_DIR/build"
RESULTS_DIR="${HEATBATH_RESULTS_DIR:-$PROJECT_DIR/data/results}"
TOPCHARGE_PARAMS="$PROJECT_DIR/input/heatbath_input/topcharge_params_su2.txt"
if [ ! -f "$TOPCHARGE_PARAMS" ] && [ -f "${TOPCHARGE_PARAMS%.txt}_finished.txt" ]; then
    TOPCHARGE_PARAMS="${TOPCHARGE_PARAMS%.txt}_finished.txt"
fi

SETUP_FILE=""
BETA=""
TARGET_Q=""
TARGET_Q_LIST=""
Q_TOL="0.1"
MAX_TRIES="100"
CANDIDATE_SWEEPS="100"
PRODUCTION_SWEEPS=""
SAVE_INTERVAL=""
SMEAR_STEPS=""
SMEAR_INTERVAL=""
SMEAR_ALPHA=""
EXCLUDE_BOUNDARY_SLICES=""
SKIP_BUILD=false
DRY_RUN=false

usage() {
    cat << EOF
Usage:
  $0 --setup FILE --beta BETA [options]

Options:
  --target-q Q              Accepted hot-sector target (default: any integer)
  --target-q-list LIST      Accept if Q is within tolerance of any value in LIST
                            Example: --target-q-list "-2,-1,0,1,2"
  --q-tol EPS               Accept if |Q - nearest int| <= EPS (default: 0.1)
  --max-tries N             Maximum hot candidates to try (default: 100)
  --candidate-sweeps N      Sweeps before checking each hot candidate (default: 100)
  --production-sweeps N     Production sweeps (default: setup num_sweeps)
  --save-interval N         Save interval (default: setup save_interval)
  --smear-steps N           Acceptance smearing steps (default: topcharge params)
  --smear-interval N        Acceptance smearing interval (default: smear_steps)
  --smear-alpha A           Acceptance APE alpha (default: topcharge params)
  --exclude-boundary N      Boundary slices excluded for Q (default: setup/topcharge)
  --results-dir DIR         Output root (default: HEATBATH_RESULTS_DIR or data/results)
  --skip-build              Do not build mc_heatbath/meas_topcharge
  --dry-run                 Print planned directories and exit
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --setup)              SETUP_FILE="$2"; shift 2 ;;
        --beta)               BETA="$2"; shift 2 ;;
        --target-q)           TARGET_Q="$2"; shift 2 ;;
        --target-q-list)      TARGET_Q_LIST="$2"; shift 2 ;;
        --q-tol)              Q_TOL="$2"; shift 2 ;;
        --max-tries)          MAX_TRIES="$2"; shift 2 ;;
        --candidate-sweeps)   CANDIDATE_SWEEPS="$2"; shift 2 ;;
        --production-sweeps)  PRODUCTION_SWEEPS="$2"; shift 2 ;;
        --save-interval)      SAVE_INTERVAL="$2"; shift 2 ;;
        --smear-steps)        SMEAR_STEPS="$2"; shift 2 ;;
        --smear-interval)     SMEAR_INTERVAL="$2"; shift 2 ;;
        --smear-alpha)        SMEAR_ALPHA="$2"; shift 2 ;;
        --exclude-boundary)   EXCLUDE_BOUNDARY_SLICES="$2"; shift 2 ;;
        --results-dir)        RESULTS_DIR="$2"; shift 2 ;;
        --skip-build)         SKIP_BUILD=true; shift ;;
        --dry-run)            DRY_RUN=true; shift ;;
        -h|--help)            usage; exit 0 ;;
        *) echo "Unknown option: $1"; usage; exit 1 ;;
    esac
done

if [ -z "$SETUP_FILE" ] || [ -z "$BETA" ]; then
    usage
    exit 1
fi

if [ ! -f "$SETUP_FILE" ]; then
    echo "Error: setup file not found: $SETUP_FILE"
    exit 1
fi

read_param() {
    grep -E "^${1}[[:space:]]+" "$2" | awk '{print $2}' | head -1
}

generate_seed() {
    python3 -c "import time; print(int(time.time()*1e9) % 2147483646 + 1)" 2>/dev/null \
        || echo $(($(date +%s) * $$ % 2147483647))
}

abs_diff_leq() {
    awk -v q="$1" -v target="$2" -v tol="$3" 'BEGIN {
        d = q - target
        if (d < 0) d = -d
        exit(d <= tol ? 0 : 1)
    }'
}

abs_diff_nearest_int_leq() {
    awk -v q="$1" -v tol="$2" 'BEGIN {
        n = q < 0 ? int(q - 0.5) : int(q + 0.5)
        d = q - n
        if (d < 0) d = -d
        exit(d <= tol ? 0 : 1)
    }'
}

abs_diff_target_list_leq() {
    awk -v q="$1" -v targets="$2" -v tol="$3" 'BEGIN {
        gsub(/[[:space:]]+/, "", targets)
        n_targets = split(targets, target, ",")
        for (i = 1; i <= n_targets; i++) {
            if (target[i] == "") continue
            d = q - target[i]
            if (d < 0) d = -d
            if (d <= tol) exit 0
        }
        exit 1
    }'
}

T_SIZE=$(read_param "T" "$SETUP_FILE"); T_SIZE=${T_SIZE:-16}
L_SIZE=$(read_param "L" "$SETUP_FILE"); L_SIZE=${L_SIZE:-16}
BOUNDARY=$(read_param "boundary" "$SETUP_FILE"); BOUNDARY=${BOUNDARY:-periodic}
PRODUCTION_SWEEPS=${PRODUCTION_SWEEPS:-$(read_param "num_sweeps" "$SETUP_FILE")}
PRODUCTION_SWEEPS=${PRODUCTION_SWEEPS:-1000}
SAVE_INTERVAL=${SAVE_INTERVAL:-$(read_param "save_interval" "$SETUP_FILE")}
SAVE_INTERVAL=${SAVE_INTERVAL:-10}

if [ -f "$TOPCHARGE_PARAMS" ]; then
    SMEAR_STEPS=${SMEAR_STEPS:-$(read_param "smear_steps" "$TOPCHARGE_PARAMS")}
    SMEAR_ALPHA=${SMEAR_ALPHA:-$(read_param "smear_alpha" "$TOPCHARGE_PARAMS")}
    if [ "$BOUNDARY" = "open" ]; then
        EXCLUDE_BOUNDARY_SLICES=${EXCLUDE_BOUNDARY_SLICES:-$(read_param "exclude_boundary_slices_open" "$TOPCHARGE_PARAMS")}
    fi
fi

SMEAR_STEPS=${SMEAR_STEPS:-40}
SMEAR_INTERVAL=${SMEAR_INTERVAL:-$SMEAR_STEPS}
SMEAR_ALPHA=${SMEAR_ALPHA:-0.5}
EXCLUDE_BOUNDARY_SLICES=${EXCLUDE_BOUNDARY_SLICES:-$(read_param "exclude_boundary_slices" "$SETUP_FILE")}
EXCLUDE_BOUNDARY_SLICES=${EXCLUDE_BOUNDARY_SLICES:-0}

RUN_STAMP="$(date +%Y%m%d_%H%M%S)"
BASE_NAME="qtarget_T${T_SIZE}_L${L_SIZE}_b${BETA}_${BOUNDARY}_${RUN_STAMP}_su2"
BASE_DIR="$RESULTS_DIR/$BASE_NAME"
CANDIDATE_ROOT="$BASE_DIR/hot_candidates"
if [ -n "$TARGET_Q_LIST" ]; then
    TARGET_LABEL="list_${TARGET_Q_LIST//,/to}"
elif [ -n "$TARGET_Q" ]; then
    TARGET_LABEL="$TARGET_Q"
else
    TARGET_LABEL="int"
fi
HOT_DIR="$BASE_DIR/hot_${TARGET_LABEL}"

write_heatbath_input() {
    local input_file="$1"
    local seed="$2"
    local start_type="$3"
    local sweeps="$4"
    local output_dir="$5"
    local config_dir="$6"
    local save_interval="$7"

    cat > "$input_file" << EOF
# Auto-generated by scripts/loop_hot_ensemble.sh
T                   $T_SIZE
L                   $L_SIZE
beta                $BETA
seed                $seed
start_type          $start_type
boundary            $BOUNDARY
num_sweeps          $sweeps
save_interval       $save_interval

output_dir          ${output_dir}/
config_dir          ${config_dir}/
EOF
}

write_topcharge_input() {
    local input_file="$1"
    local config_dir="$2"
    local output_dir="$3"
    local conf="$4"
    local seed="$5"

    cat > "$input_file" << EOF
# Auto-generated by scripts/loop_hot_ensemble.sh
config_dir              ${config_dir}/
output_dir              ${output_dir}/
beta                    $BETA
T                       $T_SIZE
L                       $L_SIZE
seed                    $seed
boundary                $BOUNDARY
start_conf              $conf
end_conf                $conf
conf_step               1
smear_steps             $SMEAR_STEPS
smear_interval          $SMEAR_INTERVAL
smear_alpha             $SMEAR_ALPHA
exclude_boundary_slices $EXCLUDE_BOUNDARY_SLICES
EOF
}

echo "SU(2) accepted-hot Q-sector ensemble driver"
echo "  setup:              $SETUP_FILE"
echo "  beta:               $BETA"
echo "  lattice:            ${T_SIZE}x${L_SIZE}^3 $BOUNDARY"
if [ "$DRY_RUN" = true ]; then
    if [ -n "$TARGET_Q" ]; then
        echo "  target hot Q:       $TARGET_Q +/- $Q_TOL"
    elif [ -n "$TARGET_Q_LIST" ]; then
        echo "  target hot Q list:  $TARGET_Q_LIST +/- $Q_TOL"
    else
        echo "  target hot Q:       any integer +/- $Q_TOL"
    fi
    exit 0
fi

mkdir -p "$BUILD_DIR" "$BASE_DIR" "$CANDIDATE_ROOT" "$HOT_DIR/output" "$HOT_DIR/configs"

if [ "$SKIP_BUILD" = false ]; then
    cd "$BUILD_DIR"
    cmake .. -DCMAKE_BUILD_TYPE=Release >/dev/null
    make -j$(nproc 2>/dev/null || echo 4) mc_heatbath meas_topcharge >/dev/null
    cd "$PROJECT_DIR"
fi

ACCEPTED_Q=""
ACCEPTED_SEED=""
ACCEPTED_CONF=""
ACCEPTED_DIR=""

for TRY in $(seq 1 "$MAX_TRIES"); do
    SEED=$(generate_seed)
    TRY_DIR="$CANDIDATE_ROOT/try_$(printf '%04d' "$TRY")_seed${SEED}"
    mkdir -p "$TRY_DIR/output" "$TRY_DIR/configs" "$TRY_DIR/topcharge"

    write_heatbath_input "$TRY_DIR/input_heatbath.txt" "$SEED" "hot" "$CANDIDATE_SWEEPS" "$TRY_DIR/output" "$TRY_DIR/configs" "$CANDIDATE_SWEEPS"
    echo "[hot try $TRY/$MAX_TRIES] generating candidate with seed=$SEED..."
    "$BUILD_DIR/bin/mc_heatbath" -i "$TRY_DIR/input_heatbath.txt" > "$TRY_DIR/heatbath.log" 2>&1

    CANDIDATE_CONF=$(printf "%04d" "$CANDIDATE_SWEEPS")
    write_topcharge_input "$TRY_DIR/input_topcharge.txt" "$TRY_DIR/configs" "$TRY_DIR/topcharge" "$CANDIDATE_SWEEPS" "$SEED"
    "$BUILD_DIR/bin/meas_topcharge" -i "$TRY_DIR/input_topcharge.txt" > "$TRY_DIR/topcharge.log" 2>&1

    Q_VALUE=$(awk 'NF >= 3 && $1 !~ /^#/ {q=$3} END {print q}' "$TRY_DIR/topcharge/topcharge.dat")
    echo "[hot try $TRY/$MAX_TRIES] Q=$Q_VALUE"

    if [ -n "$Q_VALUE" ]; then
        if [ -n "$TARGET_Q_LIST" ]; then
            accept=false
            abs_diff_target_list_leq "$Q_VALUE" "$TARGET_Q_LIST" "$Q_TOL" && accept=true
        elif [ -n "$TARGET_Q" ]; then
            accept=false
            abs_diff_leq "$Q_VALUE" "$TARGET_Q" "$Q_TOL" && accept=true
        else
            accept=false
            abs_diff_nearest_int_leq "$Q_VALUE" "$Q_TOL" && accept=true
        fi
        if [ "$accept" = true ]; then
            ACCEPTED_Q="$Q_VALUE"
            ACCEPTED_SEED="$SEED"
            ACCEPTED_CONF="$TRY_DIR/configs/conf.${CANDIDATE_CONF}"
            ACCEPTED_DIR="$TRY_DIR"
            break
        fi
    fi
done

if [ -z "$ACCEPTED_Q" ]; then
    echo "Error: no hot candidate accepted after $MAX_TRIES tries."
    exit 1
fi

cp "$ACCEPTED_CONF" "$HOT_DIR/configs/conf.0000"
HOT_SEED=$(generate_seed)
write_heatbath_input "$HOT_DIR/input.txt" "$HOT_SEED" "continue" "$PRODUCTION_SWEEPS" "$HOT_DIR/output" "$HOT_DIR/configs" "$SAVE_INTERVAL"
printf "# Sweep  Plaquette  ActionDensity\n" > "$HOT_DIR/output/plaquette.dat"
printf "# Sweep  ActionDensity\n" > "$HOT_DIR/output/action_density.dat"

echo "[hot accepted] Q=$ACCEPTED_Q from $ACCEPTED_DIR"
echo "[hot accepted] running production from accepted config..."
"$BUILD_DIR/bin/mc_heatbath" -i "$HOT_DIR/input.txt" > "$HOT_DIR/heatbath.log" 2>&1

cat > "$BASE_DIR/run_info.txt" << EOF
# Q-targeted hot ensemble run
generated              $(date)
setup_file             $SETUP_FILE
beta                   $BETA
T                      $T_SIZE
L                      $L_SIZE
boundary               $BOUNDARY
target_q               ${TARGET_Q:-any integer}
target_q_list          ${TARGET_Q_LIST:-none}
q_tolerance            $Q_TOL
accepted_q             $ACCEPTED_Q
accepted_seed          $ACCEPTED_SEED
accepted_candidate     $ACCEPTED_DIR
hot_run                $HOT_DIR

# Action density files
hot_action_density_file $HOT_DIR/output/action_density.dat
hot_plaquette_file     $HOT_DIR/output/plaquette.dat
EOF

echo "Done."
echo "  hot action density: $HOT_DIR/output/action_density.dat"
echo "  hot plaquette:      $HOT_DIR/output/plaquette.dat"
echo "  accepted Q:         $ACCEPTED_Q"
