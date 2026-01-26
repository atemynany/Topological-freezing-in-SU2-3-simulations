# SU2 Includes Library

Static library containing utility functions for SU(2) lattice gauge theory simulations.

## Build Instructions

```bash
cd Includes
mkdir -p build
cd build
cmake ..
make
```

## Build Output

| File | Location |
|------|----------|
| `libsu2_includes.a` | `build/` |
| Object files (`.o`) | `build/CMakeFiles/su2_includes.dir/` |

## Contents

### Source Files

| File | Description |
|------|-------------|
| `fields.cc` | Field operations and manipulations |
| `io.cc` | Input/output utilities for reading/writing configurations |
| `ranlux.cc` | C++ wrapper for RANLUX random number generator |
| `ranlxd.c` | RANLUX double precision RNG |
| `ranlxs.c` | RANLUX single precision RNG |
| `smearing_techniques_all.cc` | all available Smearing algorithms (APE, HYP, etc.) |
| `smearing_techniques.cc` | used Smearning algorithm |

### Header Files

| File | Description |
|------|-------------|
| `fields.hh` | Field class definitions |
| `geometry.hh` | Lattice geometry and indexing |
| `io.hh` | I/O function declarations |
| `linear_algebra.hh` | SU(2) matrix operations|
| `progressbar.hh` | Progress bar utility |
| `ranlux.hh` | Random number generator interface |
| `smearing_techniques.hh` | Smearing function declarations |

## Requirements

- CMake ≥ 3.16
- C++17 compatible compiler
- Math library (`-lm`)

## Notes

- macOS users: Uncomment `set(CMAKE_OSX_ARCHITECTURES "arm64")` in `CMakeLists.txt` for Apple Silicon