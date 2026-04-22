#!/bin/bash
# ==============================================================================
# test_parity_heatbath_topcharge.sh
# ==============================================================================
# Runs the legacy pipeline (mc_heatbath → meas_topcharge) and the fused
# pipeline (heatbath_topcharge_su2) with identical seed/lattice/beta/smearing
# on a small lattice, then diffs their outputs column-by-column.
#
# Pinned to OMP_NUM_THREADS=1 by design: with multiple threads the heatbath
# trajectory is nondeterministic (thread scheduling affects the order of RNG
# draws in the checkerboard loop), so two independent runs diverge after a
# few sweeps even with the same seed. Single-threaded parity is the real
# correctness check for the fusion.
#
# Usage:
#   scripts/test_parity_heatbath_topcharge.sh [--threads N] [--keep]
#                                             [--T N] [--L N] [--beta B]
#                                             [--sweeps N] [--save N]
#                                             [--smear N] [--alpha A]
#                                             [--boundary periodic|open]
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BIN="$PROJECT_DIR/build/bin"

# Defaults — small enough for <30s on a laptop
THREADS=1     # forced below; --threads accepted but warns
KEEP=false
T=8
L=8
BETA=2.4
NUM_SWEEPS=200
SAVE_INTERVAL=20
SMEAR_STEPS=15
SMEAR_ALPHA=0.3
BOUNDARY=periodic
SEED=20260421

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads)
            if [ "$2" != "1" ]; then
                echo "WARNING: --threads $2 requested, but parity forces OMP_NUM_THREADS=1."
                echo "         Multi-threaded heatbath trajectories are nondeterministic."
            fi
            shift 2 ;;
        --keep)     KEEP=true;          shift ;;
        --T)        T="$2";             shift 2 ;;
        --L)        L="$2";             shift 2 ;;
        --beta)     BETA="$2";          shift 2 ;;
        --sweeps)   NUM_SWEEPS="$2";    shift 2 ;;
        --save)     SAVE_INTERVAL="$2"; shift 2 ;;
        --smear)    SMEAR_STEPS="$2";   shift 2 ;;
        --alpha)    SMEAR_ALPHA="$2";   shift 2 ;;
        --boundary) BOUNDARY="$2";      shift 2 ;;
        --seed)     SEED="$2";          shift 2 ;;
        -h|--help)
            sed -n '2,17p' "$0"
            exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

EXCLUDE_BOUNDARY=0
if [ "$BOUNDARY" = "open" ]; then EXCLUDE_BOUNDARY=2; fi

GREEN='\033[0;32m'; RED='\033[0;31m'; YELLOW='\033[1;33m'; CYAN='\033[0;36m'; NC='\033[0m'

# ------------------------------------------------------------------------------
# Preconditions
# ------------------------------------------------------------------------------
for bin in mc_heatbath meas_topcharge heatbath_topcharge_su2; do
    if [ ! -x "$BIN/$bin" ]; then
        echo -e "${RED}[ERROR]${NC} Missing binary: $BIN/$bin"
        echo "Build first: cmake --build build --parallel"
        exit 1
    fi
done

TMPDIR=$(mktemp -d -t hb_tc_parity.XXXXXX)
trap '[ "$KEEP" = true ] || rm -rf "$TMPDIR"' EXIT

echo -e "${CYAN}============================================================${NC}"
echo -e "${CYAN}  heatbath_topcharge_su2 parity test${NC}"
echo -e "${CYAN}============================================================${NC}"
printf "  lattice        : %d x %d^3 (%s)\n" "$T" "$L" "$BOUNDARY"
printf "  beta / seed    : %s / %s\n" "$BETA" "$SEED"
printf "  sweeps / save  : %d / %d\n" "$NUM_SWEEPS" "$SAVE_INTERVAL"
printf "  smear / alpha  : %d / %s\n" "$SMEAR_STEPS" "$SMEAR_ALPHA"
printf "  OMP threads    : %d\n" "$THREADS"
printf "  workdir        : %s%s\n" "$TMPDIR" "$( [ "$KEEP" = true ] && echo '  (kept)')"
echo ""

export OMP_NUM_THREADS=1

# ------------------------------------------------------------------------------
# Legacy: mc_heatbath → meas_topcharge
# ------------------------------------------------------------------------------
echo -e "${YELLOW}[1/3]${NC} legacy mc_heatbath..."
mkdir -p "$TMPDIR/legacy/output" "$TMPDIR/legacy/configs"
cat > "$TMPDIR/legacy/heatbath.txt" <<EOF
T              $T
L              $L
beta           $BETA
seed           $SEED
start_type     cold
boundary       $BOUNDARY
num_sweeps     $NUM_SWEEPS
save_interval  $SAVE_INTERVAL
output_dir     $TMPDIR/legacy/output/
config_dir     $TMPDIR/legacy/configs/
EOF
"$BIN/mc_heatbath" -i "$TMPDIR/legacy/heatbath.txt" > "$TMPDIR/legacy/heatbath.log" 2>&1

echo -e "${YELLOW}[2/3]${NC} legacy meas_topcharge..."
cat > "$TMPDIR/legacy/topcharge.txt" <<EOF
config_dir               $TMPDIR/legacy/configs/
output_dir               $TMPDIR/legacy/output/
beta                     $BETA
T                        $T
L                        $L
start_conf               $SAVE_INTERVAL
end_conf                 $NUM_SWEEPS
conf_step                $SAVE_INTERVAL
smear_steps              $SMEAR_STEPS
smear_interval           $SMEAR_STEPS
smear_alpha              $SMEAR_ALPHA
seed                     $SEED
boundary                 $BOUNDARY
exclude_boundary_slices  $EXCLUDE_BOUNDARY
EOF
"$BIN/meas_topcharge" -i "$TMPDIR/legacy/topcharge.txt" > "$TMPDIR/legacy/topcharge.log" 2>&1

# ------------------------------------------------------------------------------
# Fused: heatbath_topcharge_su2
# ------------------------------------------------------------------------------
echo -e "${YELLOW}[3/3]${NC} fused heatbath_topcharge_su2..."
mkdir -p "$TMPDIR/fused/output"
cat > "$TMPDIR/fused/input.txt" <<EOF
T                        $T
L                        $L
beta                     $BETA
seed                     $SEED
start_type               cold
boundary                 $BOUNDARY
num_sweeps               $NUM_SWEEPS
save_interval            $SAVE_INTERVAL
smear_steps              $SMEAR_STEPS
smear_interval           $SMEAR_STEPS
smear_alpha              $SMEAR_ALPHA
exclude_boundary_slices  $EXCLUDE_BOUNDARY
output_dir               $TMPDIR/fused/output/
EOF
"$BIN/heatbath_topcharge_su2" -i "$TMPDIR/fused/input.txt" > "$TMPDIR/fused/run.log" 2>&1

# ------------------------------------------------------------------------------
# Compare
# ------------------------------------------------------------------------------
echo ""
echo -e "${CYAN}--- comparison ---${NC}"

FAIL=0

# 1. Plaquette trace
LEG_PLAQ="$TMPDIR/legacy/output/plaquette.dat"
FUS_PLAQ="$TMPDIR/fused/output/plaquette.dat"

if [ "$THREADS" = 1 ]; then
    # Bit-identical expected
    if diff -q "$LEG_PLAQ" "$FUS_PLAQ" > /dev/null; then
        echo -e "  plaquette.dat         ${GREEN}bit-identical${NC}"
    else
        echo -e "  plaquette.dat         ${RED}DIFFERS${NC}"
        diff "$LEG_PLAQ" "$FUS_PLAQ" | head -6
        FAIL=1
    fi
else
    # Numeric tolerance
    MAX_DIFF=$(paste <(grep -v '^#' "$LEG_PLAQ") <(grep -v '^#' "$FUS_PLAQ") \
        | awk '{d=$2-$4; if(d<0)d=-d; if(d>m)m=d} END{printf "%.3e", m+0}')
    echo "  plaquette.dat         max |Δ<P>| = $MAX_DIFF"
fi

# 2. Q column — legacy row: "smear_steps conf Q plaq"; fused row: "smear_steps sweep Q plaq"
LEG_Q=$(mktemp); FUS_Q=$(mktemp)
grep -v '^#' "$TMPDIR/legacy/output/topcharge.dat" | awk '{print $2,$3}' | sort -n > "$LEG_Q"
grep -v '^#' "$TMPDIR/fused/output/topcharge.dat"  | awk '{print $2,$3}' | sort -n > "$FUS_Q"

N_LEG=$(wc -l < "$LEG_Q"); N_FUS=$(wc -l < "$FUS_Q")
if [ "$N_LEG" -ne "$N_FUS" ]; then
    echo -e "  topcharge row counts  ${RED}legacy=$N_LEG  fused=$N_FUS${NC}"
    FAIL=1
fi

MAX_Q_DIFF=$(paste "$LEG_Q" "$FUS_Q" | awk '{
    if ($1 != $3) {print "sweep mismatch:", $1, "vs", $3 > "/dev/stderr"; exit 2}
    d=$2-$4; if(d<0)d=-d; if(d>m)m=d
} END{printf "%.3e", m+0}')
if [ "$THREADS" = 1 ]; then
    THRESHOLD=1e-10
else
    THRESHOLD=1e-7
fi
STATUS=$(awk -v d="$MAX_Q_DIFF" -v t="$THRESHOLD" 'BEGIN{print (d+0 < t+0) ? "OK" : "FAIL"}')
if [ "$STATUS" = "OK" ]; then
    echo -e "  topcharge Q           ${GREEN}max |ΔQ| = $MAX_Q_DIFF (≤ $THRESHOLD)${NC}"
else
    echo -e "  topcharge Q           ${RED}max |ΔQ| = $MAX_Q_DIFF (> $THRESHOLD)${NC}"
    FAIL=1
fi
rm -f "$LEG_Q" "$FUS_Q"

# 3. Timeslice q(t) — rows: "smear_steps sweep t q_t", keep (sweep,t,q_t)
LEG_T=$(mktemp); FUS_T=$(mktemp)
grep -v '^#' "$TMPDIR/legacy/output/topcharge_timeslice.dat" | awk '{print $2,$3,$4}' | sort -n -k1,1 -k2,2 > "$LEG_T"
grep -v '^#' "$TMPDIR/fused/output/topcharge_timeslice.dat"  | awk '{print $2,$3,$4}' | sort -n -k1,1 -k2,2 > "$FUS_T"

MAX_QT_DIFF=$(paste "$LEG_T" "$FUS_T" | awk '{
    if ($1 != $4 || $2 != $5) {print "row mismatch" > "/dev/stderr"; exit 2}
    d=$3-$6; if(d<0)d=-d; if(d>m)m=d
} END{printf "%.3e", m+0}')
STATUS=$(awk -v d="$MAX_QT_DIFF" -v t="$THRESHOLD" 'BEGIN{print (d+0 < t+0) ? "OK" : "FAIL"}')
if [ "$STATUS" = "OK" ]; then
    echo -e "  timeslice q(t)        ${GREEN}max |Δq(t)| = $MAX_QT_DIFF (≤ $THRESHOLD)${NC}"
else
    echo -e "  timeslice q(t)        ${RED}max |Δq(t)| = $MAX_QT_DIFF (> $THRESHOLD)${NC}"
    FAIL=1
fi
rm -f "$LEG_T" "$FUS_T"

echo ""
if [ "$FAIL" -eq 0 ]; then
    echo -e "${GREEN}============================================================${NC}"
    echo -e "${GREEN}  PARITY OK${NC}"
    echo -e "${GREEN}============================================================${NC}"
else
    echo -e "${RED}============================================================${NC}"
    echo -e "${RED}  PARITY FAILED — inspect $TMPDIR (use --keep)${NC}"
    echo -e "${RED}============================================================${NC}"
fi
exit $FAIL
