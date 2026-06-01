# Minimal cl2qcd Open Temporal Boundary Patch

This directory contains only the code patch needed to add open temporal gauge
boundary support to an existing `cl2qcd` checkout for the pure-gauge
`su3heatbath` workflow.

Patch:

```text
obc_patch/obc_runtime_minimal.patch
```

Apply it inside a clean/base `cl2qcd` source tree:

```bash
cd /path/to/cl2qcd
git apply /path/to/SU23_freezing/su3_static_potential_check/obc_patch/obc_runtime_minimal.patch
```

Then rebuild `su3heatbath`.

Enable OBC in the input file:

```text
useOpenTemporalGaugeBoundary=true
```

or on the command line:

```bash
./su3heatbath input_file --useOpenTemporalGaugeBoundary=true
```

The patch contains only:

- parameter parsing for `useOpenTemporalGaugeBoundary`
- forwarding that option into OpenCL kernel parameters
- the `_OPEN_TEMPORAL_GAUGE_BC_` kernel define
- temporal-boundary plaquette/staple guards
- skipping final-time temporal-link updates
- open-boundary plaquette normalization
- a sample input line documenting the option
