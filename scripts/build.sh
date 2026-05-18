#!/bin/bash
# ==============================================================================
# build.sh
# ==============================================================================
# Build script for SU(2)/SU(3) Topological Charge simulation.
# Compiles the project using CMake.
#
# Usage:
#   ./scripts/build.sh [debug|release]
#
# Author: Alexander de Barros Noll
# Date: January 2026
# ==============================================================================

set -e  # Exit on error

# Get script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="${PROJECT_ROOT}/build"

# Parse arguments
BUILD_TYPE="${1:-Release}"
MARCH_FLAG=""
for arg in "$@"; do
    case "$arg" in
        --march=*) MARCH_FLAG="${arg#--march=}" ;;
    esac
done
case "$BUILD_TYPE" in
    debug|Debug)
        BUILD_TYPE="Debug"
        ;;
    release|Release)
        BUILD_TYPE="Release"
        ;;
    *)
        echo "Unknown build type: $BUILD_TYPE"
        echo "Usage: $0 [debug|release] [--march=<arch>]"
        exit 1
        ;;
esac

echo "=============================================="
echo "Building SU(2)/SU(3) Topological Charge Simulation"
echo "=============================================="
echo "Build type: ${BUILD_TYPE}"
echo "Project root: ${PROJECT_ROOT}"
echo "Build directory: ${BUILD_DIR}"
echo ""

# Create build directory
mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

# Configure with CMake
echo "Configuring with CMake..."
CMAKE_EXTRA=""
[ -n "$MARCH_FLAG" ] && CMAKE_EXTRA="-DMARCH=${MARCH_FLAG}"
cmake -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
      -DBUILD_TESTING=ON \
      $CMAKE_EXTRA \
      "${PROJECT_ROOT}"

# Build
echo ""
echo "Building..."
cmake --build . --parallel ${BUILD_JOBS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)}

echo ""
echo "=============================================="
echo "Build complete!"
echo "=============================================="
echo "Executables:"
echo "  - ${BUILD_DIR}/bin/mc_heatbath"
echo "  - ${BUILD_DIR}/bin/mc_heatbath_su3"
echo "  - ${BUILD_DIR}/bin/meas_topcharge"
echo "  - ${BUILD_DIR}/bin/meas_topcharge_su3"
echo "  - ${BUILD_DIR}/bin/heatbath_topcharge_su2"
echo "  - ${BUILD_DIR}/bin/compute_instanton_Q"
echo ""
echo "Test executables:"
echo "  - ${BUILD_DIR}/bin/test_heatbath"
echo "  - ${BUILD_DIR}/bin/test_linear_algebra"
echo "  - ${BUILD_DIR}/bin/test_topcharge"
echo "  - ${BUILD_DIR}/bin/test_topcharge_su3"
echo "  - ${BUILD_DIR}/bin/test_smearing"
echo "  - ${BUILD_DIR}/bin/test_plaquette"
echo "  - ${BUILD_DIR}/bin/test_gauge_invariance"
echo ""
echo "To run tests: cd ${BUILD_DIR} && ctest --output-on-failure"
echo "=============================================="
