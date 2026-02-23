# CL2QCD Workflow for SU(3) Configuration Generation

This workflow uses CL2QCD (an OpenCL-based Lattice QCD code) to generate SU(3) gauge configurations, then converts them to our raw binary format for topological charge measurement.

## Overview

```
                        ┌─────────────────┐
                        │ cl2qcd_heatbath │
                        │   (our input)   │
                        └────────┬────────┘
                                 │
                                 ▼
┌─────────────────────────────────────────────────────────┐
│              convert_input_to_cl2qcd.py                 │
│           Converts to CL2QCD native format              │
└─────────────────────────┬───────────────────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────┐
│                   CL2QCD su3heatbath                    │
│       Runs Monte Carlo simulation on GPU/CPU            │
│           Outputs configs in LIME/ILDG format           │
└─────────────────────────┬───────────────────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────┐
│                     lime_to_raw                         │
│     Converts LIME format to our raw binary format       │
│        Handles index convention conversion              │
└─────────────────────────┬───────────────────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────┐
│                  meas_topcharge_su3                     │
│         Measures Q using APE smearing                   │
└─────────────────────────────────────────────────────────┘
```

## Quick Start

### 1. Build Everything

```bash
# Build CL2QCD and lime_to_raw converter
./scripts/build_cl2qcd.sh
```

This builds:
- CL2QCD `su3heatbath` executable
- `lime_to_raw` format converter
- Installs c-lime library if needed

### 2. Configure Simulation

Edit `input/cl2qcd_heatbath.txt`:

```
# Lattice dimensions
nSpace 16
nTime 16

# Action parameters
beta 6.0

# Monte Carlo parameters
nThermalizationSteps 1000
nHeatbathSteps 10000
saveFrequency 100

# Output
outputDir output/configs_cl2qcd

# Topological charge measurement
smearSteps 40
smearRho 0.12
```

### 3. Run Workflow

```bash
# Generate configs only
./scripts/run_cl2qcd_workflow.sh input/cl2qcd_heatbath.txt

# Generate configs AND measure topological charge
./scripts/run_cl2qcd_workflow.sh input/cl2qcd_heatbath.txt --measure-topcharge
```

### 4. Analyze Results

```bash
python analysis/run_analysis_su3.py output/configs_cl2qcd/topcharge_results.dat
```

## File Descriptions

| File | Description |
|------|-------------|
| `input/cl2qcd_heatbath.txt` | Main input file (unified format) |
| `scripts/build_cl2qcd.sh` | Builds CL2QCD and dependencies |
| `scripts/run_cl2qcd_workflow.sh` | Main workflow orchestration |
| `scripts/convert_input_to_cl2qcd.py` | Converts our input to CL2QCD format |
| `include/lime_reader.hh` | LIME format reader header |
| `src/lime_reader.cc` | LIME reader + format converter |

## Index Convention

CL2QCD and our code use different conventions:

| Aspect | CL2QCD | Our Code |
|--------|--------|----------|
| Site index | `x + y*L + z*L² + t*L³` | `((t*L + x)*L + y)*L + z` |
| Link index | `mu + 4*site` | `(4*site + mu) * 18` |
| Byte order | Big-endian (LIME) | Native (little-endian on x86) |

The `lime_to_raw` converter handles all these transformations automatically.

## Dependencies

- **CL2QCD**: https://github.com/ChristianReisinger/cl2qcd
- **c-lime**: https://github.com/usqcd-software/c-lime
- **OpenCL** (optional): For GPU acceleration

## Notes

- CL2QCD currently only supports **periodic boundary conditions**
- Open boundary support may be available in future CL2QCD versions
- For CPU-only builds, CL2QCD will use the CPU via OpenCL CPU runtime

## Troubleshooting

### CL2QCD build fails
```bash
# Check OpenCL installation
apt-get install ocl-icd-opencl-dev  # Ubuntu/Debian
# or
brew install opencl-headers         # macOS
```

### lime_to_raw can't find lime library
```bash
# Set library path
export LD_LIBRARY_PATH=$HOME/lime_install/lib:$LD_LIBRARY_PATH
```
