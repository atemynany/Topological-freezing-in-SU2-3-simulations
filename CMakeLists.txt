# ==============================================================================
# SU(2) Lattice Gauge Theory - Topological Charge Simulation
# ==============================================================================
# CMake build configuration for SU(2) lattice QCD simulation with
# topological charge measurement using APE smearing.
#
# Author: Alexander de Barros Noll
# Date: January 2026
# ==============================================================================

cmake_minimum_required(VERSION 3.16)

project(SU2_TopologicalCharge
    VERSION 1.0.0
    DESCRIPTION "SU(2) Lattice Gauge Theory Topological Charge Simulation"
    LANGUAGES C CXX
)

# ==============================================================================
# Build Options
# ==============================================================================

option(BUILD_TESTING "Build the testing tree" ON)
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(ENABLE_WARNINGS "Enable compiler warnings" ON)

# macOS platform support
if(APPLE)
    option(BUILD_APPLE_SILICON "Build for Apple Silicon (arm64)" OFF)
    if(BUILD_APPLE_SILICON)
        set(CMAKE_OSX_ARCHITECTURES "arm64")
    endif()
endif()

# ==============================================================================
# Compiler Configuration
# ==============================================================================

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build type defaults to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Compiler flags
if(ENABLE_WARNINGS)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        add_compile_options(-Wall -Wextra -Wpedantic)
    elseif(MSVC)
        add_compile_options(/W4)
    endif()
endif()

# OpenMP support for parallelization
# On macOS, Apple Clang does not ship with OpenMP.
# Install via: brew install libomp
if(APPLE)
    find_program(BREW_FOUND brew)
    if(BREW_FOUND)
        execute_process(
            COMMAND brew --prefix libomp
            OUTPUT_VARIABLE LIBOMP_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
            RESULT_VARIABLE BREW_LIBOMP_RESULT
        )
        if(BREW_LIBOMP_RESULT EQUAL 0 AND EXISTS "${LIBOMP_PREFIX}")
            set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include")
            set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include")
            set(OpenMP_C_LIB_NAMES "omp")
            set(OpenMP_CXX_LIB_NAMES "omp")
            set(OpenMP_omp_LIBRARY "${LIBOMP_PREFIX}/lib/libomp.dylib")
        endif()
    endif()
endif()

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found - enabling parallel computation")
else()
    message(WARNING "OpenMP not found - parallel loops will run sequentially")
    if(APPLE)
        message(STATUS "  Tip: install OpenMP on macOS with: brew install libomp")
    endif()
endif()

# Optimization flags for Release builds
# -O3: Maximum optimization
# -march=native: Use CPU-specific instructions (AVX, AVX2, etc.)
# -ffast-math: Faster floating-point operations (may reduce precision slightly)
# -funroll-loops: Unroll loops for better pipelining
# Note: On Apple Silicon, -march=native is accepted but is effectively a no-op.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(COMMON_RELEASE_FLAGS "-O3 -ffast-math -funroll-loops -DNDEBUG")
    if(NOT APPLE)
        set(COMMON_RELEASE_FLAGS "${COMMON_RELEASE_FLAGS} -march=native")
    elseif(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
        set(COMMON_RELEASE_FLAGS "${COMMON_RELEASE_FLAGS} -march=native")
    endif()
    set(CMAKE_C_FLAGS_RELEASE "${COMMON_RELEASE_FLAGS}")
    set(CMAKE_CXX_FLAGS_RELEASE "${COMMON_RELEASE_FLAGS}")
elseif(MSVC)
    set(CMAKE_C_FLAGS_RELEASE "/O2 /arch:AVX2 /fp:fast /DNDEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /arch:AVX2 /fp:fast /DNDEBUG")
else()
    set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
endif()

# Debug flags
set(CMAKE_C_FLAGS_DEBUG "-g -O0 -DDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG")

# ==============================================================================
# Output Directories
# ==============================================================================

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ==============================================================================
# Dependencies
# ==============================================================================

# On macOS, libm is part of libSystem and doesn't need explicit linking
if(NOT APPLE)
    find_library(MATH_LIBRARY m)
else()
    set(MATH_LIBRARY "")
endif()

# ==============================================================================
# Utility Library (SU2 Core Functions)
# ==============================================================================

add_library(su2_utility STATIC
    _Utility/src/fields.cc
    _Utility/src/io.cc
    _Utility/src/ranlux.cc
    _Utility/src/ranlxd.c
    _Utility/src/ranlxs.c
    _Utility/src/smearing_techniques.cc
    _Utility/src/smearing_techniques_all.cc
)

target_include_directories(su2_utility
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/_Utility/include
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(su2_utility
    PUBLIC
        ${MATH_LIBRARY}
        $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
)

# Add alias for cleaner cmake usage
add_library(SU2::Utility ALIAS su2_utility)

# ==============================================================================
# Main Executables
# ==============================================================================

# --- Monte Carlo Heatbath Configuration Generator ---
add_executable(mc_heatbath
    src/MC_heatbath.cc
)

target_include_directories(mc_heatbath
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/_Utility/include
)

target_link_libraries(mc_heatbath
    PRIVATE
        su2_utility
        ${MATH_LIBRARY}
        $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
)

# --- Topological Charge Measurement ---
add_executable(meas_topcharge
    src/meas_topcharge_su2.cc
)

target_include_directories(meas_topcharge
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/_Utility/include
)

target_link_libraries(meas_topcharge
    PRIVATE
        su2_utility
        ${MATH_LIBRARY}
        $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
)

# --- Instanton Topological Charge Computation ---
add_executable(compute_instanton_Q
    src/compute_instanton_Q.cc
    src/instanton.cc
)

target_include_directories(compute_instanton_Q
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/_Utility/include
)

target_link_libraries(compute_instanton_Q
    PRIVATE
        su2_utility
        ${MATH_LIBRARY}
        $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
)

# --- SU(3) Monte Carlo Heatbath ---
add_executable(mc_heatbath_su3
    src/MC_heatbath_su3.cc
)

target_include_directories(mc_heatbath_su3
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/_Utility/include
)

target_link_libraries(mc_heatbath_su3
    PRIVATE
        su2_utility
        ${MATH_LIBRARY}
        $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
)

# --- SU(3) Topological Charge Measurement ---
add_executable(meas_topcharge_su3
    src/meas_topcharge_su3.cc
)

target_include_directories(meas_topcharge_su3
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/_Utility/include
)

target_link_libraries(meas_topcharge_su3
    PRIVATE
        su2_utility
        ${MATH_LIBRARY}
        $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
)

# ==============================================================================
# Testing Configuration (Catch2)
# ==============================================================================

if(BUILD_TESTING)
    # Fetch Catch2 from GitHub
    include(FetchContent)
    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.5.2
    )
    FetchContent_MakeAvailable(Catch2)
    
    # Include Catch2's CMake integration
    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
    include(CTest)
    include(Catch)
    
    # --- Test: Heatbath Algorithm ---
    add_executable(test_heatbath
        tests/test_heatbath.cc
    )
    
    target_include_directories(test_heatbath
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/include
            ${CMAKE_CURRENT_SOURCE_DIR}/_Utility/include
    )
    
    target_link_libraries(test_heatbath
        PRIVATE
            su2_utility
            Catch2::Catch2WithMain
            ${MATH_LIBRARY}
    )
    
    # --- Test: SU(2) Linear Algebra ---
    add_executable(test_linear_algebra
        tests/test_linear_algebra.cc
    )
    
    target_include_directories(test_linear_algebra
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/include
            ${CMAKE_CURRENT_SOURCE_DIR}/_Utility/include
    )
    
    target_link_libraries(test_linear_algebra
        PRIVATE
            su2_utility
            Catch2::Catch2WithMain
            ${MATH_LIBRARY}
    )
    
    # --- Test: Topological Charge ---
    add_executable(test_topcharge
        tests/test_topcharge.cc
    )
    
    target_include_directories(test_topcharge
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/include
            ${CMAKE_CURRENT_SOURCE_DIR}/_Utility/include
    )
    
    target_link_libraries(test_topcharge
        PRIVATE
            su2_utility
            Catch2::Catch2WithMain
            ${MATH_LIBRARY}
    )
    
    # --- Test: Smearing ---
    add_executable(test_smearing
        tests/test_smearing.cc
    )
    
    target_include_directories(test_smearing
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/include
            ${CMAKE_CURRENT_SOURCE_DIR}/_Utility/include
    )
    
    target_link_libraries(test_smearing
        PRIVATE
            su2_utility
            Catch2::Catch2WithMain
            ${MATH_LIBRARY}
    )
    
    # --- Test: Plaquette ---
    add_executable(test_plaquette
        tests/test_plaquette.cc
    )
    
    target_include_directories(test_plaquette
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/include
            ${CMAKE_CURRENT_SOURCE_DIR}/_Utility/include
    )
    
    target_link_libraries(test_plaquette
        PRIVATE
            su2_utility
            Catch2::Catch2WithMain
            ${MATH_LIBRARY}
    )
    
    # Register tests with CTest
    # Using simple test registration to avoid issues with special characters in test names
    add_test(NAME test_heatbath COMMAND test_heatbath)
    add_test(NAME test_linear_algebra COMMAND test_linear_algebra)
    add_test(NAME test_topcharge COMMAND test_topcharge)
    add_test(NAME test_smearing COMMAND test_smearing)
    add_test(NAME test_plaquette COMMAND test_plaquette)
endif()

# ==============================================================================
# Installation
# ==============================================================================

install(TARGETS mc_heatbath meas_topcharge compute_instanton_Q
    RUNTIME DESTINATION bin
)

install(TARGETS su2_utility
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)

install(DIRECTORY _Utility/include/
    DESTINATION include/su2
    FILES_MATCHING PATTERN "*.hh" PATTERN "*.h"
)

install(DIRECTORY include/
    DESTINATION include/su2
    FILES_MATCHING PATTERN "*.hh" PATTERN "*.h"
)

# ==============================================================================
# Summary
# ==============================================================================

message(STATUS "")
message(STATUS "=== SU(2) Topological Charge Build Configuration ===")
message(STATUS "Build type:        ${CMAKE_BUILD_TYPE}")
message(STATUS "C compiler:        ${CMAKE_C_COMPILER}")
message(STATUS "C++ compiler:      ${CMAKE_CXX_COMPILER}")
message(STATUS "Build testing:     ${BUILD_TESTING}")
message(STATUS "Install prefix:    ${CMAKE_INSTALL_PREFIX}")
message(STATUS "System:            ${CMAKE_SYSTEM_NAME} (${CMAKE_SYSTEM_PROCESSOR})")
if(APPLE)
    message(STATUS "macOS support:     enabled")
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP (macOS):    found via Homebrew")
    else()
        message(STATUS "OpenMP (macOS):    not found (install with: brew install libomp)")
    endif()
endif()
message(STATUS "")
