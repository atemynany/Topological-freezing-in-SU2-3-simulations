##### THIS README WAS AUTOMATICALY GENRATED WITH CLAUDE OPUS 4.6 AND MAY HAVE ERRORS SINCE ITS NOT REVIEWED OR CORRECTLY MAINTAINED !!!!!!!! I WILL REWRITE IT AS SOON AS THE PROJECT IS FINISHED!!!!!!! ####



# SU(2) Lattice Gauge Theory: Topological Charge Simulation 

Numerical simulation of topological charge in SU(2) lattice gauge theory using Monte Carlo methods with heatbath updates and APE smearing.

## Overview

This project implements:
- Monte Carlo configuration generation using the heatbath algorithm for SU(2) gauge theory
- Topological charge measurement using the clover (field strength tensor) definition
- APE smearing for UV fluctuation reduction
- Statistical analysis and visualization tools

## Requirements

### Build Requirements
- CMake >= 3.16
- C++17 compatible compiler (GCC >= 7, Clang >= 5)
- C11 compatible compiler
- Math library (libm)

### Optional
- Python >= 3.7 (for analysis and plotting)
- NumPy, Matplotlib (for visualization)
- SLURM (for cluster job submission)

## Directory Structure

```
.
├── CMakeLists.txt          # Main CMake configuration
├── _Utility/               # Core SU(2) library
│   ├── include/            # Header files
│   │   ├── fields.hh       # Gauge field allocation/manipulation
│   │   ├── geometry.hh     # Lattice indexing
│   │   ├── io.hh           # Configuration I/O
│   │   ├── linear_algebra.hh # SU(2) matrix operations
│   │   ├── ranlux.hh       # Random number generator
│   │   └── smearing_techniques.hh # APE smearing
│   └── src/                # Implementation files
├── include/                # Project headers
│   ├── Wilson_loops.hh     # Plaquette calculations
│   └── topcharge_su2.hh    # Topological charge measurement
├── src/                    # Main source files
│   ├── MC_heatbath.cc      # Configuration generator
│   └── meas_topcharge_su2.cc # Measurement program
├── tests/                  # Unit tests (Catch2)
│   ├── test_heatbath.cc
│   ├── test_linear_algebra.cc
│   ├── test_topcharge.cc
│   └── test_smearing.cc
├── scripts/                # Shell scripts
│   ├── build.sh            # Build script
│   ├── run_heatbath.sh     # Configuration generation
│   ├── run_topcharge.sh    # Measurement job
│   ├── run_tests.sh        # Test runner
│   └── run_full_analysis.sh # Complete workflow
├── analysis/               # Python analysis tools
│   ├── plot_topcharge.py   # Visualization
│   └── analyze_topcharge.py # Statistical analysis
├── params/                 # Parameter files
│   └── example.conf        # Example configuration
└── input/                  # Input file templates
    └── topcharge_example.txt
```

## Build Instructions

### Quick Build

```bash
./scripts/build.sh release
```

### Manual Build

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --parallel
```

### Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `CMAKE_BUILD_TYPE` | Release | Build type (Debug/Release) |
| `BUILD_TESTING` | ON | Build unit tests |

## Usage

### 1. Generate Gauge Configurations

```bash
./build/bin/mc_heatbath <output_path> <beta> <T> <L> <seed> <cold/hot> <num_sweeps> <output_interval> <periodic/open>
```

Example:
```bash
./build/bin/mc_heatbath data/configs/ 2.4 16 16 12345 cold 10000 500 periodic
```

### 2. Measure Topological Charge

Create an input file (see `input/topcharge_example.txt`):
```
config_dir          data/configs/T16_L16_b2.4
output_dir          data/results/T16_L16_b2.4
beta                2.4
T                   16
L                   16
start_conf          500
end_conf            10000
conf_step           500
smear_steps         200
smear_output_interval 20
smear_alpha         0.5
```

Run measurement:
```bash
./build/bin/meas_topcharge -i input/topcharge_example.txt
```

### 3. Analyze and Plot Results

```bash
python3 analysis/plot_topcharge.py --input data/results/topcharge_*.dat --output plots/
python3 analysis/analyze_topcharge.py --input data/results/topcharge_*.dat --output analysis_results.txt
```

### 4. Complete Workflow

```bash
./scripts/run_full_analysis.sh params/example.conf
```

## Physical Background

### Topological Charge

The topological charge Q is defined as:

$$Q = \frac{1}{32\pi^2} \int d^4x \, \epsilon_{\mu\nu\rho\sigma} \text{Tr}(F_{\mu\nu} F_{\rho\sigma})$$

On the lattice, using the clover definition:

$$Q = \frac{1}{4\pi^2} \sum_x \sum_{\mu<\nu,\rho<\sigma} \epsilon_{\mu\nu\rho\sigma} \text{Tr}(\text{Im}(C_{\mu\nu}) \text{Im}(C_{\rho\sigma}))$$

where $C_{\mu\nu}$ is the clover (4-leaf plaquette) average.

### APE Smearing

APE smearing reduces UV fluctuations while preserving topological information:

$$U_\mu^{\text{new}}(x) = \text{Proj}_{SU(2)}\left[ U_\mu(x) + \alpha \sum_{\nu \neq \mu} S_{\mu\nu}(x) \right]$$

where $S_{\mu\nu}$ is the staple sum.

## Testing

Run all tests:
```bash
./scripts/run_tests.sh
```

Run specific test:
```bash
./scripts/run_tests.sh test_topcharge
```

Tests verify:
- SU(2) matrix operations
- Geometry and indexing
- Plaquette and clover calculations
- Smearing SU(2) preservation
- Topological charge properties

## Cluster Usage (SLURM)

Submit heatbath job:
```bash
sbatch scripts/run_heatbath.sh params/production.conf
```

Submit measurement job:
```bash
sbatch scripts/run_topcharge.sh input/measurement.txt
```

## Output Format

### Configuration Files
Binary format with header containing simulation parameters.

### Topological Charge Data
```
# smear_steps  config_number  Q  plaquette
    0    500    0.123456    0.567890
   20    500   -0.054321    0.789012
  ...
```

## References

1. M. Luscher, "Computational Strategies in Lattice QCD"
2. C. Morningstar, M. Peardon, "Analytic smearing of SU(3) link variables in lattice QCD"
3. APE Collaboration, "Glueball Masses and String Tension in Lattice QCD"

## License

For academic use. See accompanying license file.

## Author

Alexander de Barros Noll, 2026
