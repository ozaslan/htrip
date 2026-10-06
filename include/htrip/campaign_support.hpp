#pragma once

#include "grid_map.hpp"
#include "minplus_oracle.hpp"
#include "types.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace htrip {
namespace campaign {

/**
 * @brief Summary statistics for benchmark latency distributions.
 */
struct TimingStats {
    double median_ms = 0.0;  ///< Sample median execution time (ms).
    double mean_ms = 0.0;    ///< Arithmetic sample mean execution time (ms).
    double min_ms = 0.0;     ///< Minimum recorded latency (ms).
    double max_ms = 0.0;     ///< Maximum recorded latency (ms).
    double q1_ms = 0.0;      ///< First quartile (ms).
    double q3_ms = 0.0;      ///< Third quartile (ms).
};

/**
 * @brief Connected component of passable space extracted from a ground-truth map.
 *
 * Used to ensure that sampled query viewpoints belong to the same topologically connected
 * traversable subgraph, preventing trivial infinite distances across disconnected rooms.
 */
struct ConnectedComponent {
    int id = -1;                 ///< Component identifier.
    size_t size = 0;             ///< Number of passable cells in this component.
    std::vector<int> cell_comp;  ///< Flattened component map (rows * cols, -1 if impassable).
};

/**
 * @brief Result descriptor for exactness verification between two all-pairs distance matrices.
 */
struct MatrixCheck {
    size_t entries = 0;          ///< Total number of matrix elements compared.
    size_t mismatches = 0;       ///< Number of elements differing by more than 0.
    int max_abs_diff = 0;        ///< Maximum absolute difference observed across all pairs.
    dist_t max_finite = 0;       ///< Maximum finite distance value observed.
};

/**
 * @brief Configuration descriptor for a standard benchmark map used in paper evaluations.
 */
struct PaperMap {
    const char* id;              ///< Short identifier string (e.g. "boston256", "maze512").
    const char* path;            ///< Relative path to MovingAI .map file.
    int b;                       ///< Calibrated leaf tile size b (e.g. 16 or 32).
    int g;                       ///< Calibrated internal node grouping factor g (e.g. 2).
};

inline constexpr PaperMap kPaperMaps[] = {
    {"boston256", "datasets/movingai/street/maps/Boston_0_256.map", 16, 2},
    {"boston1024", "datasets/movingai/street/maps/Boston_0_1024.map", 32, 2},
    {"maze512-1", "datasets/movingai/maze/maps/maze512-1-0.map", 16, 2},
    {"maze512-16", "datasets/movingai/maze/maps/maze512-16-0.map", 16, 2},
    {"paris512", "datasets/movingai/street/maps/Paris_0_512.map", 8, 2},
    {"berlin1024", "datasets/movingai/street/maps/Berlin_0_1024.map", 8, 2},
};

std::string jsonEscape(const std::string& s);
std::string shellCapture(const std::string& cmd);
std::string findRepoRoot();
std::string resolveExistingPath(const std::string& path);
std::string fileSha256(const std::string& path);
std::string gitSha(const std::string& repo);
bool gitDirty(const std::string& repo);
std::int64_t peakRssBytes();

TimingStats summarizeMs(std::vector<double> samples);
std::vector<int> parseIntList(const std::string& csv, const std::vector<int>& fallback);
std::string joinInts(const std::vector<int>& v);

bool loadFullyRevealed(const std::string& path, GridMap& grid);
size_t countGtPassable(const GridMap& grid);
ConnectedComponent largestGtComponent(const GridMap& grid);
void markObservedFreeExcept(GridMap& grid, const std::vector<Point>& withheld);

std::vector<Point> sampleSector(const GridMap& grid, const ConnectedComponent& cc, int divisions = 8,
                                int origin_r = 0, int origin_c = 0);
std::vector<Point> sampleSpread(const GridMap& grid, const ConnectedComponent& cc, std::uint64_t seed, size_t n);
std::vector<Point> sampleCluster(const GridMap& grid, const ConnectedComponent& cc, std::uint64_t seed, size_t n);
std::vector<Point> sampleGeometry(const std::string& name, const GridMap& grid, const ConnectedComponent& cc,
                                  std::uint64_t seed, size_t n, int origin_r = 0, int origin_c = 0);

MatrixCheck compareMatrices(const std::vector<dist_t>& a, const std::vector<dist_t>& b);

bool csvNeedsHeader(const std::string& path);
void csvAppendLine(const std::string& path, const std::string& header, const std::string& line);

std::string exactnessHeader();
std::string queryVsKHeader();
std::string queryVsKRepeatsHeader();
std::string updateHeader();
std::string rssHeader();
std::string fullBuildHeader();

} // namespace campaign
} // namespace htrip
