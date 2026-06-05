# SU(2) Utility Library

`_Utility/` contains the shared SU(2)-oriented support code used by the root
CMake project. It is built as the static library `su2_utility` and linked into
the SU(2) executables, the SU(3) executables where they reuse common
infrastructure, and the Catch2 tests.

Build it through the repository root:

```bash
./scripts/build.sh release
```

The archive is written to `build/lib/libsu2_utility.a`; executables are written
to `build/bin/`.

## Contents

| File | Purpose |
|------|---------|
| `src/fields.cc` / `include/fields.hh` | Gauge-field allocation, copying, and basic field helpers |
| `src/io.cc` / `include/io.hh` | Parameter parsing and binary gauge-configuration I/O |
| `src/ranlux.cc`, `src/ranlxd.c`, `src/ranlxs.c` / `include/ranlux.hh` | RANLUX random-number generators |
| `src/smearing_techniques.cc` / `include/smearing_techniques.hh` | APE smearing used in production measurements |
| `src/smearing_techniques_all.cc` | Additional smearing variants retained for experiments |
| `include/geometry.hh` | Lattice indexing and neighbor-table helpers |
| `include/linear_algebra.hh` | SU(2) matrix operations |
| `include/progressbar.hh` | Lightweight terminal progress bar |

## Notes

- The library assumes the global lattice-state pattern used by the simulation
  binaries (`T`, `L`, gauge-field pointer, boundary flags, neighbor tables).
- OpenMP is enabled when CMake finds it. On macOS, install it with
  `brew install libomp`; otherwise the code still builds and runs sequentially.
- SU(3)-specific linear algebra and topcharge code lives in the root
  `include/` directory (`su3_linear_algebra.hh`, `topcharge_su3.hh`,
  `su3_heatbath.hh`).
