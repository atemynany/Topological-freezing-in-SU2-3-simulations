#!/bin/bash
# ==============================================================================
# run_tests.sh
# ==============================================================================
# Run all unit tests for the SU(2) topological charge simulation.
#
# Usage:
#   ./scripts/run_tests.sh [test_name]
#
# Author: Alexander de Barros Noll
# Date: January 2026
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="${PROJECT_ROOT}/build"

# Build if necessary
if [ ! -d "$BUILD_DIR" ]; then
    echo "Build directory not found. Running build script..."
    "${SCRIPT_DIR}/build.sh"
fi

cd "${BUILD_DIR}"

# Check if tests were built
if [ ! -f "bin/test_heatbath" ]; then
    echo "Tests not built. Rebuilding with tests enabled..."
    cmake -DBUILD_TESTING=ON "${PROJECT_ROOT}"
    cmake --build . --parallel $(nproc 2>/dev/null || echo 4)
fi

echo "=============================================="
echo "Running SU(2) Unit Tests"
echo "=============================================="
echo ""

if [ -n "$1" ]; then
    # Run specific test
    echo "Running test: $1"
    ctest --output-on-failure -R "$1"
else
    # Run all tests
    echo "Running all tests..."
    ctest --output-on-failure
fi

echo ""
echo "=============================================="
echo "Tests complete!"
echo "=============================================="
