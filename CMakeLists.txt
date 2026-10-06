cmake_minimum_required(VERSION 3.20)
project(htrip VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Target-based modern options
if(MSVC)
    set(_icra_opts /W4 /O2 /Oi /Ot)
else()
    set(_icra_opts -Wall -Wextra -Wpedantic -O3 -march=native)
endif()

enable_testing()

# Core static library - ZERO dependencies on other hbrick modules
find_package(Threads REQUIRED)

add_library(htrip STATIC
    src/grid_map.cpp
    src/frontier.cpp
    src/leaf_tile.cpp
    src/internal_node.cpp
    src/minplus_oracle.cpp
    src/bfs_oracle.cpp
    src/astar_oracle.cpp
    src/fuel_oracle.cpp
    src/dijkstra_oracle.cpp
    src/tour_planner.cpp
    src/robot_simulator.cpp
    src/mission_logger.cpp
    src/arena_recorder.cpp
    src/trial_manifest.cpp
    src/campaign_support.cpp
)

target_include_directories(htrip PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_options(htrip PRIVATE ${_icra_opts})
target_link_libraries(htrip PUBLIC Threads::Threads)

# Automated verification executable
add_executable(test_exactness apps/test_exactness.cpp)
target_link_libraries(test_exactness PRIVATE htrip)
target_compile_options(test_exactness PRIVATE ${_icra_opts})

add_test(NAME test_exactness COMMAND test_exactness)

# FUEL baseline comparison executable
add_executable(test_fuel_baseline apps/test_fuel_baseline.cpp)
target_link_libraries(test_fuel_baseline PRIVATE htrip)
target_compile_options(test_fuel_baseline PRIVATE ${_icra_opts})

add_test(NAME test_fuel_baseline COMMAND test_fuel_baseline)

# Benchmark runner executable
add_executable(bench_frontier_oracle apps/bench_frontier_oracle.cpp)
target_link_libraries(bench_frontier_oracle PRIVATE htrip)
target_compile_options(bench_frontier_oracle PRIVATE ${_icra_opts})

# Property-based differential test (randomized map mutations)
add_executable(test_property_based apps/test_property_based.cpp)
target_link_libraries(test_property_based PRIVATE htrip)
target_compile_options(test_property_based PRIVATE ${_icra_opts})

add_test(NAME test_property_based COMMAND test_property_based --seeds 6 --size 96 --steps 20)

# Google Benchmark micro/macro benchmarks (fetched on demand)
find_package(benchmark CONFIG QUIET)
if(NOT benchmark_FOUND)
    include(FetchContent)
    set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
    set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
    set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(
        benchmark
        URL https://github.com/google/benchmark/archive/refs/tags/v1.8.5.zip
    )
    FetchContent_MakeAvailable(benchmark)
endif()

add_executable(bench_micro apps/bench_micro.cpp)
target_link_libraries(bench_micro PRIVATE htrip benchmark::benchmark)
target_compile_options(bench_micro PRIVATE ${_icra_opts})

# Weighted graph unit test
add_executable(test_weighted apps/test_weighted.cpp)
target_link_libraries(test_weighted PRIVATE htrip)
target_compile_options(test_weighted PRIVATE ${_icra_opts})

add_test(NAME test_weighted COMMAND test_weighted)

# Weighted graph benchmark (construction, single cell update, 10-cell tile batch update)
add_executable(bench_weighted apps/bench_weighted.cpp)
target_link_libraries(bench_weighted PRIVATE htrip)
target_compile_options(bench_weighted PRIVATE ${_icra_opts})

# Systematic BFS comparison on large-scale mazes (256, 512, 1024)
add_executable(bench_systematic_bfs apps/bench_systematic_bfs.cpp)
target_link_libraries(bench_systematic_bfs PRIVATE htrip)
target_compile_options(bench_systematic_bfs PRIVATE ${_icra_opts})

# Reproducer suite for code review bugs
add_executable(test_reproducers apps/test_reproducers.cpp)
target_link_libraries(test_reproducers PRIVATE htrip)
target_compile_options(test_reproducers PRIVATE ${_icra_opts})

add_test(NAME test_reproducers COMMAND test_reproducers)

# Multi-robot continuous SE(2) exploration unit test
add_executable(test_robot_simulator apps/test_robot_simulator.cpp)
target_link_libraries(test_robot_simulator PRIVATE htrip)
target_compile_options(test_robot_simulator PRIVATE ${_icra_opts})

add_test(NAME test_robot_simulator COMMAND test_robot_simulator)

# Multi-robot exploration CLI benchmark runner
add_executable(run_multi_robot_sim apps/run_multi_robot_sim.cpp)
target_link_libraries(run_multi_robot_sim PRIVATE htrip)
target_compile_options(run_multi_robot_sim PRIVATE ${_icra_opts})

# ICRA 2027 results campaign (E0–E6). Cold-path CSV drivers; not hot query kernels.
add_executable(campaign_e0 apps/campaign_e0.cpp)
target_link_libraries(campaign_e0 PRIVATE htrip)
target_compile_options(campaign_e0 PRIVATE ${_icra_opts})

add_executable(campaign_exactness apps/campaign_exactness.cpp)
target_link_libraries(campaign_exactness PRIVATE htrip)
target_compile_options(campaign_exactness PRIVATE ${_icra_opts})

add_executable(bench_query_vs_k apps/bench_query_vs_k.cpp)
target_link_libraries(bench_query_vs_k PRIVATE htrip)
target_compile_options(bench_query_vs_k PRIVATE ${_icra_opts})

add_executable(bench_update_latency apps/bench_update_latency.cpp)
target_link_libraries(bench_update_latency PRIVATE htrip)
target_compile_options(bench_update_latency PRIVATE ${_icra_opts})

