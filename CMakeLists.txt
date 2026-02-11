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

# macOS Apple Silicon support (uncomment if needed)
# set(CMAKE_OSX_ARCHITECTURES "arm64")

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

# Optimization flags for Release builds
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

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

find_library(MATH_LIBRARY m)

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
message(STATUS "")
