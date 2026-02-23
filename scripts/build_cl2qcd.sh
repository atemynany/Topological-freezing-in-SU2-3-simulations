#!/bin/bash
# build_cl2qcd.sh - Build CL2QCD and dependencies
#
# This script builds:
# 1. c-lime library (for reading LIME format files)
# 2. CL2QCD su3heatbath executable
# 3. lime_to_raw converter
#
# Usage: ./scripts/build_cl2qcd.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="$PROJECT_DIR/build"
DEPS_DIR="$PROJECT_DIR/su3_imports"
CL2QCD_DIR="$DEPS_DIR/cl2qcd_reisinger"

echo "=============================================="
echo "Building CL2QCD and Dependencies"
echo "=============================================="
echo "Project dir: $PROJECT_DIR"
echo "Build dir: $BUILD_DIR"
echo

# -----------------------------------------------------------------------------
# 1. Check/Install c-lime library
# -----------------------------------------------------------------------------
echo "[1/3] Checking c-lime library..."

LIME_INSTALLED=false

# Check if lime is available
if pkg-config --exists lime 2>/dev/null; then
    echo "  c-lime found via pkg-config"
    LIME_INSTALLED=true
elif [ -f "/usr/local/lib/liblime.a" ] || [ -f "/usr/lib/liblime.a" ]; then
    echo "  c-lime found in system"
    LIME_INSTALLED=true
fi

if [ "$LIME_INSTALLED" = false ]; then
    echo "  c-lime not found, building from source..."
    
    LIME_SRC="$DEPS_DIR/c-lime"
    LIME_INSTALL="$DEPS_DIR/lime_install"
    
    if [ ! -d "$LIME_SRC" ]; then
        echo "  Cloning c-lime..."
        cd "$DEPS_DIR"
        git clone https://github.com/usqcd-software/c-lime.git
    fi
    
    # Use CMake to build c-lime (easier than autotools)
    LIME_BUILD="$LIME_SRC/build"
    mkdir -p "$LIME_BUILD"
    cd "$LIME_BUILD"
    
    echo "  Configuring c-lime with CMake..."
    cmake .. -DCMAKE_INSTALL_PREFIX="$LIME_INSTALL" -DCMAKE_BUILD_TYPE=Release
    
    echo "  Building c-lime..."
    make -j$(nproc)
    make install
    
    export LIME_DIR="$LIME_INSTALL"
    export PKG_CONFIG_PATH="$LIME_DIR/lib/pkgconfig:$PKG_CONFIG_PATH"
    export LD_LIBRARY_PATH="$LIME_DIR/lib:$LD_LIBRARY_PATH"
    
    echo "  c-lime installed to $LIME_DIR"
fi

echo

# -----------------------------------------------------------------------------
# 2. Build CL2QCD
# -----------------------------------------------------------------------------
echo "[2/3] Building CL2QCD su3heatbath..."

if [ ! -d "$CL2QCD_DIR" ]; then
    echo "  Error: CL2QCD not found at $CL2QCD_DIR"
    echo "  Please clone it first:"
    echo "    cd $DEPS_DIR"
    echo "    git clone https://github.com/ChristianReisinger/cl2qcd.git cl2qcd_reisinger"
    exit 1
fi

# Check for required system dependencies
echo "  Checking system dependencies..."
MISSING_DEPS=""
pkg-config --exists libxml-2.0 || MISSING_DEPS="$MISSING_DEPS libxml2-dev"
[ -f "/usr/include/gmp.h" ] || MISSING_DEPS="$MISSING_DEPS libgmp-dev"
[ -f "/usr/include/mpfr.h" ] || MISSING_DEPS="$MISSING_DEPS libmpfr-dev"
[ -f "/usr/include/nettle/nettle-types.h" ] || MISSING_DEPS="$MISSING_DEPS nettle-dev"
[ -f "/usr/include/CL/cl.h" ] || MISSING_DEPS="$MISSING_DEPS ocl-icd-opencl-dev"

if [ -n "$MISSING_DEPS" ]; then
    echo "  Missing dependencies:$MISSING_DEPS"
    echo "  Install with: sudo apt-get install -y$MISSING_DEPS"
    exit 1
fi

CL2QCD_BUILD="$CL2QCD_DIR/build"
mkdir -p "$CL2QCD_BUILD"
cd "$CL2QCD_BUILD"

# Get LIME paths
LIME_INSTALL="$DEPS_DIR/lime_install"

# Configure with CMake
echo "  Configuring CL2QCD..."
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DLIME_INCLUDE_DIR="$LIME_INSTALL/include" \
    -DLIME_LIBRARIES="$LIME_INSTALL/lib/liblime.a" \
    2>&1 | tail -10

# Build su3heatbath
echo "  Building su3heatbath executable..."
make -j$(nproc) su3heatbath 2>&1 | tail -10

if [ -f "$CL2QCD_BUILD/su3heatbath" ]; then
    echo "  Success: $CL2QCD_BUILD/su3heatbath"
    
    # Create symlink in project bin
    mkdir -p "$BUILD_DIR/bin"
    ln -sf "$CL2QCD_BUILD/su3heatbath" "$BUILD_DIR/bin/cl2qcd_heatbath"
    echo "  Symlinked to: $BUILD_DIR/bin/cl2qcd_heatbath"
else
    echo "  Error: su3heatbath build failed"
    echo "  Check the CL2QCD build output for errors"
    exit 1
fi

echo

# -----------------------------------------------------------------------------
# 3. Build lime_to_raw converter
# -----------------------------------------------------------------------------
echo "[3/3] Building lime_to_raw converter..."

cd "$PROJECT_DIR"

# Get lime flags
if [ -n "$LIME_DIR" ]; then
    LIME_CFLAGS="-I$LIME_DIR/include"
    LIME_LIBS="-L$LIME_DIR/lib -llime"
else
    LIME_CFLAGS=$(pkg-config --cflags lime 2>/dev/null || echo "")
    LIME_LIBS=$(pkg-config --libs lime 2>/dev/null || echo "-llime")
fi

echo "  Compiling lime_reader.cc..."
g++ -O3 -std=c++11 -DLIME_READER_MAIN \
    -I"$PROJECT_DIR/include" \
    $LIME_CFLAGS \
    "$PROJECT_DIR/src/lime_reader.cc" \
    -o "$BUILD_DIR/bin/lime_to_raw" \
    $LIME_LIBS

if [ -f "$BUILD_DIR/bin/lime_to_raw" ]; then
    echo "  Success: $BUILD_DIR/bin/lime_to_raw"
else
    echo "  Error: lime_to_raw build failed"
    exit 1
fi

echo
echo "=============================================="
echo "Build Complete!"
echo "=============================================="
echo
echo "Executables:"
echo "  $BUILD_DIR/bin/cl2qcd_heatbath  - CL2QCD su3heatbath"
echo "  $BUILD_DIR/bin/lime_to_raw      - LIME to raw converter"
echo
echo "Next steps:"
echo "  1. Run the workflow:"
echo "     ./scripts/run_cl2qcd_workflow.sh input/cl2qcd_heatbath.txt"
echo
