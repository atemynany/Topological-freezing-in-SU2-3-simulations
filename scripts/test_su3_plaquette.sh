#!/bin/bash
# Test SU(3) plaquette against literature values

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$(dirname "$SCRIPT_DIR")"

echo "=== SU(3) Plaquette Validation ==="
echo "Reference: Montvay & Münster, 'Quantum Fields on a Lattice'"
echo ""

for beta in 5.7 6.0 6.2; do
  cat > /tmp/test_beta.txt << EOF
T 8
L 8
beta $beta
seed 54321
start_type cold
num_sweeps 200
save_interval 200
overrelax_steps 0
EOF
  
  result=$(./build/bin/mc_heatbath_su3 -i /tmp/test_beta.txt 2>&1 | grep "Final <P>" | awk '{print $4}')
  echo "beta=$beta: <P> = $result"
done

echo ""
echo "Expected (literature, infinite volume):"
echo "  beta=5.7: ~0.547"
echo "  beta=6.0: ~0.593"  
echo "  beta=6.2: ~0.614"
