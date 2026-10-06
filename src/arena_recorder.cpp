#include "htrip/arena_recorder.hpp"
#include <cstdio>
#include <cmath>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <cstring>

namespace htrip {

namespace {

inline void setPixel(std::vector<uint8_t>& buf, int W, int H, int x, int y, uint8_t r, uint8_t g, uint8_t b) noexcept {
    if (x < 0 || x >= W || y < 0 || y >= H) return;
    size_t idx = static_cast<size_t>((y * W + x) * 3);
    buf[idx + 0] = r;
    buf[idx + 1] = g;
    buf[idx + 2] = b;
}

inline void blendPixel(std::vector<uint8_t>& buf, int W, int H, int x, int y, uint8_t r, uint8_t g, uint8_t b, float alpha) noexcept {
    if (x < 0 || x >= W || y < 0 || y >= H) return;
    size_t idx = static_cast<size_t>((y * W + x) * 3);
    buf[idx + 0] = static_cast<uint8_t>(buf[idx + 0] * (1.0f - alpha) + r * alpha);
    buf[idx + 1] = static_cast<uint8_t>(buf[idx + 1] * (1.0f - alpha) + g * alpha);
    buf[idx + 2] = static_cast<uint8_t>(buf[idx + 2] * (1.0f - alpha) + b * alpha);
}

void drawLine(std::vector<uint8_t>& buf, int W, int H, int x0, int y0, int x1, int y1,
              uint8_t r, uint8_t g, uint8_t b, int thickness = 1, float alpha = 1.0f) noexcept {
    int dx = std::abs(x1 - x0);
    int dy = std::abs(y1 - y0);
    int sx = (x0 < x1) ? 1 : -1;
    int sy = (y0 < y1) ? 1 : -1;
    int err = dx - dy;

    int rad = (thickness > 1) ? (thickness / 2) : 0;

    while (true) {
        if (rad == 0) {
            if (alpha >= 0.99f) setPixel(buf, W, H, x0, y0, r, g, b);
            else blendPixel(buf, W, H, x0, y0, r, g, b, alpha);
        } else {
            for (int dy_ = -rad; dy_ <= rad; ++dy_) {
                for (int dx_ = -rad; dx_ <= rad; ++dx_) {
                    if (alpha >= 0.99f) setPixel(buf, W, H, x0 + dx_, y0 + dy_, r, g, b);
                    else blendPixel(buf, W, H, x0 + dx_, y0 + dy_, r, g, b, alpha);
                }
            }
        }

        if (x0 == x1 && y0 == y1) break;
        int e2 = 2 * err;
        if (e2 > -dy) { err -= dy; x0 += sx; }
        if (e2 < dx) { err += dx; y0 += sy; }
    }
}

void drawCircleFilled(std::vector<uint8_t>& buf, int W, int H, int cx, int cy, int radius,
                      uint8_t r, uint8_t g, uint8_t b) noexcept {
    int r2 = radius * radius;
    for (int dy = -radius; dy <= radius; ++dy) {
        for (int dx = -radius; dx <= radius; ++dx) {
            if (dx * dx + dy * dy <= r2) {
                setPixel(buf, W, H, cx + dx, cy + dy, r, g, b);
            }
        }
    }
}

void drawCircleOutline(std::vector<uint8_t>& buf, int W, int H, int cx, int cy, int radius,
                       uint8_t r, uint8_t g, uint8_t b) noexcept {
    int x = radius;
    int y = 0;
    int err = 0;

    while (x >= y) {
        setPixel(buf, W, H, cx + x, cy + y, r, g, b);
        setPixel(buf, W, H, cx + y, cy + x, r, g, b);
        setPixel(buf, W, H, cx - y, cy + x, r, g, b);
        setPixel(buf, W, H, cx - x, cy + y, r, g, b);
        setPixel(buf, W, H, cx - x, cy - y, r, g, b);
        setPixel(buf, W, H, cx - y, cy - x, r, g, b);
        setPixel(buf, W, H, cx + y, cy - x, r, g, b);
        setPixel(buf, W, H, cx + x, cy - y, r, g, b);

        y += 1;
        err += 1 + 2 * y;
        if (2 * (err - x) + 1 > 0) {
            x -= 1;
            err += 1 - 2 * x;
        }
    }
}

// Compact 5x7 bitmap font for HUD telemetry text
static const uint8_t font5x7[96][5] = {
    {0x00, 0x00, 0x00, 0x00, 0x00}, // ' '
    {0x00, 0x00, 0x5F, 0x00, 0x00}, // '!'
    {0x00, 0x07, 0x00, 0x07, 0x00}, // '"'
    {0x14, 0x7F, 0x14, 0x7F, 0x14}, // '#'
    {0x24, 0x2A, 0x7F, 0x2A, 0x12}, // '$'
    {0x23, 0x13, 0x08, 0x64, 0x62}, // '%'
    {0x36, 0x49, 0x55, 0x22, 0x50}, // '&'
    {0x00, 0x05, 0x03, 0x00, 0x00}, // '''
    {0x00, 0x1C, 0x22, 0x41, 0x00}, // '('
    {0x00, 0x41, 0x22, 0x1C, 0x00}, // ')'
    {0x14, 0x08, 0x3E, 0x08, 0x14}, // '*'
    {0x08, 0x08, 0x3E, 0x08, 0x08}, // '+'
    {0x00, 0x50, 0x30, 0x00, 0x00}, // ','
    {0x08, 0x08, 0x08, 0x08, 0x08}, // '-'
    {0x00, 0x60, 0x60, 0x00, 0x00}, // '.'
    {0x20, 0x10, 0x08, 0x04, 0x02}, // '/'
    {0x3E, 0x51, 0x49, 0x45, 0x3E}, // '0'
    {0x00, 0x42, 0x7F, 0x40, 0x00}, // '1'
    {0x42, 0x61, 0x51, 0x49, 0x46}, // '2'
    {0x21, 0x41, 0x45, 0x4B, 0x31}, // '3'
    {0x18, 0x14, 0x12, 0x7F, 0x10}, // '4'
    {0x27, 0x45, 0x45, 0x45, 0x39}, // '5'
    {0x3C, 0x4A, 0x49, 0x49, 0x30}, // '6'
    {0x01, 0x71, 0x09, 0x05, 0x03}, // '7'
    {0x36, 0x49, 0x49, 0x49, 0x36}, // '8'
    {0x06, 0x49, 0x49, 0x29, 0x1E}, // '9'
    {0x00, 0x36, 0x36, 0x00, 0x00}, // ':'
    {0x00, 0x56, 0x36, 0x00, 0x00}, // ';'
    {0x08, 0x14, 0x22, 0x41, 0x00}, // '<'
    {0x14, 0x14, 0x14, 0x14, 0x14}, // '='
    {0x00, 0x41, 0x22, 0x14, 0x08}, // '>'
    {0x02, 0x01, 0x51, 0x09, 0x06}, // '?'
    {0x32, 0x49, 0x79, 0x41, 0x3E}, // '@'
    {0x7E, 0x11, 0x11, 0x11, 0x7E}, // 'A'
    {0x7F, 0x49, 0x49, 0x49, 0x36}, // 'B'
    {0x3E, 0x41, 0x41, 0x41, 0x22}, // 'C'
    {0x7F, 0x41, 0x41, 0x22, 0x1C}, // 'D'
    {0x7F, 0x49, 0x49, 0x49, 0x41}, // 'E'
    {0x7F, 0x09, 0x09, 0x09, 0x01}, // 'F'
    {0x3E, 0x41, 0x49, 0x49, 0x7A}, // 'G'
    {0x7F, 0x08, 0x08, 0x08, 0x7F}, // 'H'
    {0x00, 0x41, 0x7F, 0x41, 0x00}, // 'I'
    {0x20, 0x40, 0x41, 0x3F, 0x01}, // 'J'
    {0x7F, 0x08, 0x14, 0x22, 0x41}, // 'K'
    {0x7F, 0x40, 0x40, 0x40, 0x40}, // 'L'
    {0x7F, 0x02, 0x0C, 0x02, 0x7F}, // 'M'
    {0x7F, 0x04, 0x08, 0x10, 0x7F}, // 'N'
    {0x3E, 0x41, 0x41, 0x41, 0x3E}, // 'O'
    {0x7F, 0x09, 0x09, 0x09, 0x06}, // 'P'
    {0x3E, 0x41, 0x51, 0x21, 0x5E}, // 'Q'
    {0x7F, 0x09, 0x19, 0x29, 0x46}, // 'R'
    {0x46, 0x49, 0x49, 0x49, 0x31}, // 'S'
    {0x01, 0x01, 0x7F, 0x01, 0x01}, // 'T'
    {0x3F, 0x40, 0x40, 0x40, 0x3F}, // 'U'
    {0x1F, 0x20, 0x40, 0x20, 0x1F}, // 'V'
    {0x3F, 0x40, 0x38, 0x40, 0x3F}, // 'W'
    {0x63, 0x14, 0x08, 0x14, 0x63}, // 'X'
    {0x07, 0x08, 0x70, 0x08, 0x07}, // 'Y'
    {0x61, 0x51, 0x49, 0x45, 0x43}, // 'Z'
    {0x00, 0x7F, 0x41, 0x41, 0x00}, // '['
    {0x02, 0x04, 0x08, 0x10, 0x20}, // '\'
    {0x00, 0x41, 0x41, 0x7F, 0x00}, // ']'
    {0x04, 0x02, 0x01, 0x02, 0x04}, // '^'
    {0x40, 0x40, 0x40, 0x40, 0x40}, // '_'
    {0x00, 0x01, 0x02, 0x04, 0x00}, // '`'
    {0x20, 0x54, 0x54, 0x54, 0x78}, // 'a'
    {0x7F, 0x48, 0x44, 0x44, 0x38}, // 'b'
    {0x38, 0x44, 0x44, 0x44, 0x20}, // 'c'
    {0x38, 0x44, 0x44, 0x48, 0x7F}, // 'd'
    {0x38, 0x54, 0x54, 0x54, 0x18}, // 'e'
    {0x08, 0x7E, 0x09, 0x01, 0x02}, // 'f'
    {0x0C, 0x52, 0x52, 0x52, 0x3E}, // 'g'
    {0x7F, 0x08, 0x04, 0x04, 0x78}, // 'h'
    {0x00, 0x44, 0x7D, 0x40, 0x00}, // 'i'
    {0x20, 0x40, 0x44, 0x3D, 0x00}, // 'j'
    {0x7F, 0x10, 0x28, 0x44, 0x00}, // 'k'
    {0x00, 0x41, 0x7F, 0x40, 0x00}, // 'l'
    {0x7C, 0x04, 0x18, 0x04, 0x78}, // 'm'
    {0x7C, 0x08, 0x04, 0x04, 0x78}, // 'n'
    {0x38, 0x44, 0x44, 0x44, 0x38}, // 'o'
    {0x7C, 0x14, 0x14, 0x14, 0x08}, // 'p'
    {0x08, 0x14, 0x14, 0x18, 0x7C}, // 'q'
    {0x7C, 0x08, 0x04, 0x04, 0x08}, // 'r'
    {0x48, 0x54, 0x54, 0x54, 0x20}, // 's'
    {0x04, 0x3F, 0x44, 0x40, 0x20}, // 't'
    {0x3C, 0x40, 0x40, 0x20, 0x7C}, // 'u'
    {0x1C, 0x20, 0x40, 0x20, 0x1C}, // 'v'
    {0x3C, 0x40, 0x30, 0x40, 0x3C}, // 'w'
    {0x44, 0x28, 0x10, 0x28, 0x44}, // 'x'
    {0x0C, 0x50, 0x50, 0x50, 0x3C}, // 'y'
    {0x44, 0x64, 0x54, 0x4C, 0x44}, // 'z'
    {0x00, 0x08, 0x36, 0x41, 0x00}, // '{'
    {0x00, 0x00, 0x7F, 0x00, 0x00}, // '|'
    {0x00, 0x41, 0x36, 0x08, 0x00}, // '}'
    {0x08, 0x08, 0x2A, 0x1C, 0x08}  // '~'
};

void drawChar5x7(std::vector<uint8_t>& buf, int W, int H, int x, int y, char c,
                 uint8_t r, uint8_t g, uint8_t b, int scale = 1) noexcept {
    if (c < 32 || c > 126) c = '?';
    int c_idx = c - 32;
    for (int col = 0; col < 5; ++col) {
        uint8_t line = font5x7[c_idx][col];
        for (int row = 0; row < 7; ++row) {
            if (line & (1 << row)) {
                if (scale <= 1) {
                    setPixel(buf, W, H, x + col, y + row, r, g, b);
                } else {
                    for (int dy = 0; dy < scale; ++dy) {
                        for (int dx = 0; dx < scale; ++dx) {
                            setPixel(buf, W, H, x + col * scale + dx, y + row * scale + dy, r, g, b);
                        }
                    }
                }
            }
        }
    }
}

void drawString5x7(std::vector<uint8_t>& buf, int W, int H, int x, int y, const char* str,
                   uint8_t r, uint8_t g, uint8_t b, int scale = 1) noexcept {
    if (!str) return;
    int cur_x = x;
    int step = 6 * scale;
    for (const char* p = str; *p; ++p) {
        drawChar5x7(buf, W, H, cur_x, y, *p, r, g, b, scale);
        cur_x += step;
    }
}

inline void drawString5x7(std::vector<uint8_t>& buf, int W, int H, int x, int y, const std::string& str,
                          uint8_t r, uint8_t g, uint8_t b, int scale = 1) noexcept {
    drawString5x7(buf, W, H, x, y, str.c_str(), r, g, b, scale);
}

} // namespace

/**
 * @brief Rasterizes the continuous SE(2) multi-robot simulation state into an RGB image buffer.
 *
 * Implements a layered visualization pipeline designed for scientific figures and videos:
 * 1. Cell Layer:
 *    - Revealed Obstacles (dark graphite).
 *    - Explored Free Space (clean warm light gray floor).
 *    - Fog of War: Unexplored free corridors (dark navy) and unrevealed obstacle shadows (blueprint slate).
 * 2. Spatial Hierarchy Layer:
 *    - Quadtree tile grid lines.
 *    - Dynamic update glow: Base tile updates (amber border), H-TRIP hierarchical index updates (emerald border).
 * 3. LiDAR Rays:
 *    - Semi-transparent ray paths and obstacle hit endpoints.
 * 4. Exploration Frontiers:
 *    - Frontier cluster convex hulls and clearance-safe stand-off viewpoints.
 * 5. Multi-Robot Fleet:
 *    - Historical motion trails.
 *    - Planned obstacle-cleared corridor paths.
 *    - Circular robot footprints, heading vectors, and recovery status badges.
 * 6. Telemetry HUD:
 *    - Top banner displaying real-time exploration progress, H-TRIP vs BFS/A* speedups, and microsecond latencies.
 */
void rasterizeArena(const htrip::MultiRobotSimulator& sim, int width, int height, std::vector<uint8_t>& out_rgb) {
    int R = sim.observed.rows;
    int C = sim.observed.cols;
    if (R <= 0 || C <= 0) return;

    // Make width and height even numbers for H.264
    int W = (width / 2) * 2;
    int H = (height / 2) * 2;
    size_t total_bytes = static_cast<size_t>(W * H * 3);
    if (out_rgb.size() != total_bytes) {
        out_rgb.resize(total_bytes);
    }

    const int font_scale = (W >= 1024) ? 2 : 1;
    const int hud_h = (font_scale >= 2) ? 84 : 48; // Top HUD banner height (spacious 3-line layout)
    int arena_h = H - hud_h;

    float scale_x = static_cast<float>(W) / static_cast<float>(C);
    float scale_y = static_cast<float>(arena_h) / static_cast<float>(R);

    // Active leaf tiles lookup: Base tile updates vs H-TRIP index updates
    int b = sim.oracle.config.tile_size_b > 0 ? sim.oracle.config.tile_size_b : 16;
    int num_tx = (C + b - 1) / b;
    std::vector<uint8_t> base_tile_decay;
    std::vector<uint8_t> index_tile_decay;
    sim.getTileDecaySnapshots(base_tile_decay, index_tile_decay);

    // 1. Render Cells (Obstacles, Fog of War, Free Floor)
    for (int r = 0; r < R; ++r) {
        int y0 = hud_h + static_cast<int>(r * scale_y);
        int y1 = hud_h + static_cast<int>((r + 1) * scale_y);
        y1 = std::max(y1, y0 + 1);

        for (int c = 0; c < C; ++c) {
            int x0 = static_cast<int>(c * scale_x);
            int x1 = static_cast<int>((c + 1) * scale_x);
            x1 = std::max(x1, x0 + 1);

            auto state = sim.observed.getObserved(r, c);
            uint8_t cr, cg, cb;
            if (state == htrip::CellState::Obstacle) {
                cr = 32; cg = 34; cb = 40; // Revealed solid obstacle wall (dark graphite)
            } else if (state == htrip::CellState::Unknown) {
                // Environment map shadow / blueprint silhouette in the background
                bool gt_is_obstacle = (sim.ground_truth.getObserved(r, c) == htrip::CellState::Obstacle);
                if (gt_is_obstacle) {
                    // Obstacle wall shadow (blueprint steel-slate shadow under fog)
                    cr = 44; cg = 52; cb = 70;
                } else {
                    // Deep dark navy background for unexplored free corridors
                    cr = 14; cg = 16; cb = 22;
                }
            } else {
                // Explored clean floor (active highlights rendered as dedicated overlays)
                cr = 226; cg = 224; cb = 214;
            }

            for (int y = y0; y < y1 && y < H; ++y) {
                uint8_t* ptr = &out_rgb[(y * W + x0) * 3];
                for (int x = x0; x < x1 && x < W; ++x) {
                    ptr[0] = cr;
                    ptr[1] = cg;
                    ptr[2] = cb;
                    ptr += 3;
                }
            }
        }
    }

    // 2. Subtle Quadtree Tile Boundaries
    for (int r = 0; r <= R; r += b) {
        int y = hud_h + static_cast<int>(r * scale_y);
        if (y >= hud_h && y < H) {
            for (int x = 0; x < W; ++x) {
                blendPixel(out_rgb, W, H, x, y, 255, 230, 0, 0.25f);
            }
        }
    }
    for (int c = 0; c <= C; c += b) {
        int x = static_cast<int>(c * scale_x);
        if (x >= 0 && x < W) {
            for (int y = hud_h; y < H; ++y) {
                blendPixel(out_rgb, W, H, x, y, 255, 230, 0, 0.25f);
            }
        }
    }

    // Glowing highlight boundaries for active base updating leaf tiles (Warm Amber - SADECE ÇERÇEVE / FRAME ONLY)
    for (size_t tid = 0; tid < base_tile_decay.size(); ++tid) {
        if (base_tile_decay[tid] == 0) continue;
        int tx = static_cast<int>(tid % static_cast<size_t>(num_tx));
        int ty = static_cast<int>(tid / static_cast<size_t>(num_tx));
        int r0 = ty * b;
        int c0 = tx * b;
        int r1 = std::min(r0 + b, R);
        int c1 = std::min(c0 + b, C);

        int bx0 = static_cast<int>(c0 * scale_x);
        int by0 = hud_h + static_cast<int>(r0 * scale_y);
        int bx1 = static_cast<int>(c1 * scale_x);
        int by1 = hud_h + static_cast<int>(r1 * scale_y);

        float alpha = std::clamp(static_cast<float>(base_tile_decay[tid]) / 16.0f, 0.2f, 0.95f);
        for (int x = bx0; x <= bx1 && x < W; ++x) {
            blendPixel(out_rgb, W, H, x, by0, 255, 140, 0, alpha);
            if (by0 + 1 < H) blendPixel(out_rgb, W, H, x, by0 + 1, 255, 140, 0, alpha * 0.7f);
            if (by1 - 1 >= hud_h) blendPixel(out_rgb, W, H, x, by1 - 1, 255, 140, 0, alpha * 0.7f);
            blendPixel(out_rgb, W, H, x, by1, 255, 140, 0, alpha);
        }
        for (int y = by0; y <= by1 && y < H; ++y) {
            blendPixel(out_rgb, W, H, bx0, y, 255, 140, 0, alpha);
            if (bx0 + 1 < W) blendPixel(out_rgb, W, H, bx0 + 1, y, 255, 140, 0, alpha * 0.7f);
            if (bx1 - 1 >= 0) blendPixel(out_rgb, W, H, bx1 - 1, y, 255, 140, 0, alpha * 0.7f);
            blendPixel(out_rgb, W, H, bx1, y, 255, 140, 0, alpha);
        }
    }

    // Glowing highlight interiors for H-TRIP index updated leaf tiles (Electric Emerald - SADECE İÇ KISIM / INTERIOR ONLY WITH 5% MARGIN)
    for (size_t tid = 0; tid < index_tile_decay.size(); ++tid) {
        if (index_tile_decay[tid] == 0) continue;
        int tx = static_cast<int>(tid % static_cast<size_t>(num_tx));
        int ty = static_cast<int>(tid / static_cast<size_t>(num_tx));
        int r0 = ty * b;
        int c0 = tx * b;
        int r1 = std::min(r0 + b, R);
        int c1 = std::min(c0 + b, C);

        int bx0 = static_cast<int>(c0 * scale_x);
        int by0 = hud_h + static_cast<int>(r0 * scale_y);
        int bx1 = static_cast<int>(c1 * scale_x);
        int by1 = hud_h + static_cast<int>(r1 * scale_y);

        float tw = static_cast<float>(bx1 - bx0);
        float th = static_cast<float>(by1 - by0);
        int mx = std::max(1, static_cast<int>(std::round(tw * 0.05f)));
        int my = std::max(1, static_cast<int>(std::round(th * 0.05f)));

        int ix0 = bx0 + mx;
        int iy0 = by0 + my;
        int ix1 = bx1 - mx;
        int iy1 = by1 - my;

        float alpha = std::clamp((static_cast<float>(index_tile_decay[tid]) / 20.0f) * 0.45f, 0.1f, 0.55f);
        for (int y = iy0; y <= iy1 && y < H; ++y) {
            for (int x = ix0; x <= ix1 && x < W; ++x) {
                blendPixel(out_rgb, W, H, x, y, 0, 240, 160, alpha);
            }
        }
    }

    // 3. Robot Exploration Trails (Dynamic comet tail with fading alpha)
    for (const auto& bot : sim.robots) {
        if (bot.trail.size() < 2) continue;
        uint8_t tr_r = static_cast<uint8_t>(bot.color & 0xFF);
        uint8_t tr_g = static_cast<uint8_t>((bot.color >> 8) & 0xFF);
        uint8_t tr_b = static_cast<uint8_t>((bot.color >> 16) & 0xFF);

        size_t count = std::min(bot.trail.size(), size_t{80});
        size_t start_i = bot.trail.size() - count;
        for (size_t i = start_i + 1; i < bot.trail.size(); ++i) {
            float progress = static_cast<float>(i - start_i) / static_cast<float>(count);
            float trail_alpha = 0.15f + 0.65f * progress;
            int x0 = static_cast<int>((bot.trail[i - 1].first / sim.cell_size_m) * scale_x);
            int y0 = hud_h + static_cast<int>((bot.trail[i - 1].second / sim.cell_size_m) * scale_y);
            int x1 = static_cast<int>((bot.trail[i].first / sim.cell_size_m) * scale_x);
            int y1 = hud_h + static_cast<int>((bot.trail[i].second / sim.cell_size_m) * scale_y);
            drawLine(out_rgb, W, H, x0, y0, x1, y1, tr_r, tr_g, tr_b, 2, trail_alpha);
        }
    }

    // 4. LiDAR Laser Scanners (High-Visibility Sensor Fan & Obstacle Impact Sparks)
    for (const auto& bot : sim.robots) {
        int rx = static_cast<int>((bot.x / sim.cell_size_m) * scale_x);
        int ry = hud_h + static_cast<int>((bot.y / sim.cell_size_m) * scale_y);
        
        // Boost laser beam color for vivid illumination
        uint8_t bot_r = static_cast<uint8_t>(bot.color & 0xFF);
        uint8_t bot_g = static_cast<uint8_t>((bot.color >> 8) & 0xFF);
        uint8_t bot_b = static_cast<uint8_t>((bot.color >> 16) & 0xFF);
        uint8_t l_r = static_cast<uint8_t>(std::min(255, static_cast<int>(bot_r) + 70));
        uint8_t l_g = static_cast<uint8_t>(std::min(255, static_cast<int>(bot_g) + 70));
        uint8_t l_b = static_cast<uint8_t>(std::min(255, static_cast<int>(bot_b) + 70));

        // Draw outer scan boundary arc
        int prev_hx = -1, prev_hy = -1;

        // Draw all beams with prominent alpha so the laser cone/sweep is clearly visible
        for (size_t i = 0; i < bot.scan.size(); ++i) {
            const auto& beam = bot.scan[i];
            int hx = static_cast<int>((beam.hit_x / sim.cell_size_m) * scale_x);
            int hy = hud_h + static_cast<int>((beam.hit_y / sim.cell_size_m) * scale_y);

            // Laser beam ray: alpha 0.55f for strong visibility
            drawLine(out_rgb, W, H, rx, ry, hx, hy, l_r, l_g, l_b, 1, 0.55f);

            // Connect scan perimeter arc between adjacent beams
            if (prev_hx >= 0 && prev_hy >= 0) {
                drawLine(out_rgb, W, H, prev_hx, prev_hy, hx, hy, l_r, l_g, l_b, 1, 0.65f);
            }
            prev_hx = hx;
            prev_hy = hy;

            // Obstacle hit endpoints: bright laser impact points
            if (beam.hit_obstacle) {
                // Bright glowing impact spark on the obstacle wall
                drawCircleFilled(out_rgb, W, H, hx, hy, 2, 255, 230, 80); // Bright gold impact spark
                setPixel(out_rgb, W, H, hx, hy, 255, 255, 255);            // White central core
            } else {
                // Free space termination dot (max range reached)
                blendPixel(out_rgb, W, H, hx, hy, l_r, l_g, l_b, 0.75f);
            }
        }
    }

    // 5. Robot Planned Paths to Assigned Frontiers
    for (const auto& bot : sim.robots) {
        if (!bot.has_goal || bot.planned_path.empty()) continue;
        uint8_t p_r = static_cast<uint8_t>(bot.color & 0xFF);
        uint8_t p_g = static_cast<uint8_t>((bot.color >> 8) & 0xFF);
        uint8_t p_b = static_cast<uint8_t>((bot.color >> 16) & 0xFF);

        int prev_x = static_cast<int>((bot.x / sim.cell_size_m) * scale_x);
        int prev_y = hud_h + static_cast<int>((bot.y / sim.cell_size_m) * scale_y);

        size_t start_idx = std::min(bot.path_index, bot.planned_path.size() - 1);
        for (size_t i = start_idx; i < bot.planned_path.size(); ++i) {
            int px = static_cast<int>((bot.planned_path[i].c + 0.5f) * scale_x);
            int py = hud_h + static_cast<int>((bot.planned_path[i].r + 0.5f) * scale_y);
            drawLine(out_rgb, W, H, prev_x, prev_y, px, py, p_r, p_g, p_b, 2, 0.85f);
            prev_x = px;
            prev_y = py;
        }
    }

    // 5b. Frontier Cluster Boundary Glow (Highlights live exploration wavefront)
    for (const auto& f : sim.frontiers) {
        for (const auto& pt : f.cells) {
            int px = static_cast<int>((pt.c + 0.5f) * scale_x);
            int py = hud_h + static_cast<int>((pt.r + 0.5f) * scale_y);
            blendPixel(out_rgb, W, H, px, py, 0, 240, 255, 0.45f);
        }
    }

    // 6. Frontier Viewpoint Markers (Clean circular style with 2x size: fr_rad 6-12)
    for (size_t i = 0; i < sim.viewpoints.size(); ++i) {
        const auto& vp = sim.viewpoints[i];
        int cx = static_cast<int>((vp.c + 0.5f) * scale_x);
        int cy = hud_h + static_cast<int>((vp.r + 0.5f) * scale_y);
        int fr_rad = std::clamp(static_cast<int>(scale_x * 1.8f), 6, 12); // Exactly 2x previous size (was scale_x * 0.9f)

        // Check if assigned to any robot
        bool is_assigned = false;
        uint8_t fr_r = 0, fr_g = 230, fr_b = 255;
        for (const auto& bot : sim.robots) {
            if (bot.has_goal && bot.assigned_frontier == vp) {
                fr_r = static_cast<uint8_t>(bot.color & 0xFF);
                fr_g = static_cast<uint8_t>((bot.color >> 8) & 0xFF);
                fr_b = static_cast<uint8_t>((bot.color >> 16) & 0xFF);
                is_assigned = true;
                break;
            }
        }

        // Draw outer halo and filled core
        drawCircleOutline(out_rgb, W, H, cx, cy, fr_rad + 2, fr_r, fr_g, fr_b);
        drawCircleFilled(out_rgb, W, H, cx, cy, fr_rad, fr_r, fr_g, fr_b);
        drawCircleOutline(out_rgb, W, H, cx, cy, fr_rad, 20, 20, 20);
        // Inner bright white core dot
        drawCircleFilled(out_rgb, W, H, cx, cy, std::max(2, fr_rad / 3), 255, 255, 255);

        if (is_assigned) {
            // Distinct targeting reticle for assigned frontiers
            drawCircleOutline(out_rgb, W, H, cx, cy, fr_rad + 4, 255, 255, 255);
            drawLine(out_rgb, W, H, cx - fr_rad - 4, cy, cx + fr_rad + 4, cy, 255, 255, 255, 1);
            drawLine(out_rgb, W, H, cx, cy - fr_rad - 4, cx, cy + fr_rad + 4, 255, 255, 255, 1);
        }
    }

    // 7. Robots (Chassis, Heading Needle, ID Badge)
    for (const auto& bot : sim.robots) {
        int rx = static_cast<int>((bot.x / sim.cell_size_m) * scale_x);
        int ry = hud_h + static_cast<int>((bot.y / sim.cell_size_m) * scale_y);
        int bot_rad = std::clamp(static_cast<int>((sim.robot_config.radius_m / sim.cell_size_m) * scale_x), 4, 14);

        uint8_t b_r = static_cast<uint8_t>(bot.color & 0xFF);
        uint8_t b_g = static_cast<uint8_t>((bot.color >> 8) & 0xFF);
        uint8_t b_b = static_cast<uint8_t>((bot.color >> 16) & 0xFF);

        // Body
        drawCircleFilled(out_rgb, W, H, rx, ry, bot_rad, b_r, b_g, b_b);
        drawCircleOutline(out_rgb, W, H, rx, ry, bot_rad, 15, 15, 15);

        // Heading arrow
        int arrow_len = bot_rad + 6;
        int ax = rx + static_cast<int>(arrow_len * std::cos(bot.theta));
        int ay = ry + static_cast<int>(arrow_len * std::sin(bot.theta));
        drawLine(out_rgb, W, H, rx, ry, ax, ay, 255, 255, 255, 2);

        const bool recovering = (bot.recovery_phase != htrip::RobotState::RecoveryPhase::None);
        const bool controller_stuck = (bot.stuck_ticks >= 6);
        const int warn_thresh = std::max(20, sim.auto_respawn_ticks / 2);
        const bool translation_stuck = (bot.stuck_translation_ticks >= warn_thresh);
        if (recovering || controller_stuck || translation_stuck) {
            const bool critical = recovering || (bot.stuck_translation_ticks >= sim.auto_respawn_ticks * 3 / 4);
            const uint8_t wr = 255;
            const uint8_t wg = critical ? 40 : 140;
            const uint8_t wb = 20;
            const int ring = bot_rad + 4;
            drawCircleOutline(out_rgb, W, H, rx, ry, ring, wr, wg, wb);
            drawCircleOutline(out_rgb, W, H, rx, ry, ring + 2, wr, wg, wb);
            // Four warning ticks around the chassis
            for (int k = 0; k < 4; ++k) {
                double ang = bot.theta + k * (M_PI * 0.5);
                int tx0 = rx + static_cast<int>((ring + 1) * std::cos(ang));
                int ty0 = ry + static_cast<int>((ring + 1) * std::sin(ang));
                int tx1 = rx + static_cast<int>((ring + 7) * std::cos(ang));
                int ty1 = ry + static_cast<int>((ring + 7) * std::sin(ang));
                drawLine(out_rgb, W, H, tx0, ty0, tx1, ty1, wr, wg, wb, 2);
            }
            const char* tag = "STUCK";
            if (bot.recovery_phase == htrip::RobotState::RecoveryPhase::BackUp) tag = "BACKUP";
            else if (bot.recovery_phase == htrip::RobotState::RecoveryPhase::RepulsiveSpin) tag = "SPIN";
            else if (bot.recovery_phase == htrip::RobotState::RecoveryPhase::ReplanNudge) tag = "REPLAN";
            char badge[48];
            std::snprintf(badge, sizeof(badge), "%s %d/%d", tag,
                          bot.stuck_translation_ticks, sim.auto_respawn_ticks);
            const int text_x = rx - static_cast<int>(std::strlen(badge) * 3);
            const int text_y = std::max(hud_h + 2, ry - bot_rad - 16);
            // Opaque label plate
            for (int dy = -2; dy < 10; ++dy) {
                for (int dx = -3; dx < static_cast<int>(std::strlen(badge) * 6) + 3; ++dx) {
                    setPixel(out_rgb, W, H, text_x + dx, text_y + dy, 160, 16, 16);
                }
            }
            drawString5x7(out_rgb, W, H, text_x, text_y, badge, 255, 230, 80);
            if (bot.recovery_phase == htrip::RobotState::RecoveryPhase::BackUp) {
                int bx = rx - static_cast<int>((bot_rad + 10) * std::cos(bot.theta));
                int by = ry - static_cast<int>((bot_rad + 10) * std::sin(bot.theta));
                drawLine(out_rgb, W, H, rx, ry, bx, by, 255, 80, 40, 3);
            }
        }
    }

    // 8. Top HUD Telemetry Banner
    for (int y = 0; y < hud_h; ++y) {
        uint8_t* ptr = &out_rgb[(y * W) * 3];
        for (int x = 0; x < W; ++x) {
            ptr[0] = 13; ptr[1] = 16; ptr[2] = 22; // Sleek dark cyber-slate HUD background
            ptr += 3;
        }
    }
    // Bottom border of HUD with accent glow
    for (int x = 0; x < W; ++x) {
        setPixel(out_rgb, W, H, x, hud_h - 2, 28, 48, 75);
        setPixel(out_rgb, W, H, x, hud_h - 1, 0, 190, 245); // Electric cyan divider
    }

    // HUD Text formatting (3 cleanly formatted telemetry lines)
    auto formatUs = [](double us, char* out, size_t sz) {
        if (us <= 0.0) {
            std::snprintf(out, sz, "--");
        } else if (us < 1000.0) {
            std::snprintf(out, sz, "%.0fus", us);
        } else if (us < 10000.0) {
            std::snprintf(out, sz, "%.2fms", us / 1000.0);
        } else if (us < 1000000.0) {
            std::snprintf(out, sz, "%.1fms", us / 1000.0);
        } else {
            std::snprintf(out, sz, "%.2fs", us / 1000000.0);
        }
    };

    double frac = sim.explorationFraction() * 100.0;
    double sim_time_s = static_cast<double>(sim.total_ticks) * 0.05;

    char map_upd_str[32];
    formatUs(sim.last_htrip_update_us, map_upd_str, sizeof(map_upd_str));

    int stuck_n = 0;
    for (const auto& bot : sim.robots) {
        if (bot.recovery_phase != RobotState::RecoveryPhase::None ||
            bot.stuck_ticks >= 6 ||
            bot.stuck_translation_ticks >= std::max(20, sim.auto_respawn_ticks / 2)) {
            stuck_n++;
        }
    }

    const int y1 = (font_scale >= 2) ? 8 : 4;
    const int y2 = (font_scale >= 2) ? 32 : 18;
    const int y3 = (font_scale >= 2) ? 56 : 32;
    const int left_pad = (font_scale >= 2) ? 16 : 8;

    char line1[256];
    if (sim.total_robot_respawns > 0) {
        std::snprintf(line1, sizeof(line1),
                      "Step: %llu (%.1fs) | Sim: %.0fHz | Upd: %s | Fleet: %d | Stuck: %d (Resp: %u)",
                      static_cast<unsigned long long>(sim.total_ticks),
                      sim_time_s,
                      sim.sim_step_freq_hz,
                      map_upd_str,
                      static_cast<int>(sim.robots.size()),
                      stuck_n,
                      sim.total_robot_respawns);
    } else {
        std::snprintf(line1, sizeof(line1),
                      "Step: %llu (%.1fs) | Sim: %.0fHz | Upd: %s | Fleet: %d | Stuck: %d",
                      static_cast<unsigned long long>(sim.total_ticks),
                      sim_time_s,
                      sim.sim_step_freq_hz,
                      map_upd_str,
                      static_cast<int>(sim.robots.size()),
                      stuck_n);
    }

    char line2[256];
    if (!sim.enable_fog_of_war) {
        std::snprintf(line2, sizeof(line2),
                      "Map: %dx%d (%.1fm) | Targets: %d Active | Visited: %zu | Known Floor",
                      sim.ground_truth.cols, sim.ground_truth.rows,
                      sim.world_w_m,
                      static_cast<int>(sim.frontiers.size()),
                      sim.targets_visited_count);
    } else {
        std::snprintf(line2, sizeof(line2),
                      "Map: %dx%d (%.1fm) | Frontiers: %d Active | Cov: %.1f%% (%d/%d)",
                      sim.ground_truth.cols, sim.ground_truth.rows,
                      sim.world_w_m,
                      static_cast<int>(sim.frontiers.size()),
                      frac,
                      sim.explored_free_cells, sim.total_free_cells_gt);
    }

    drawString5x7(out_rgb, W, H, left_pad, y1, line1, 0, 230, 255, font_scale);
    drawString5x7(out_rgb, W, H, left_pad, y2, line2, 195, 220, 255, font_scale);

    char htrip_qry_str[32];
    formatUs(sim.last_htrip_query_us, htrip_qry_str, sizeof(htrip_qry_str));

    char bfs_qry_str[32];
    formatUs(sim.last_bfs_query_us, bfs_qry_str, sizeof(bfs_qry_str));

    char astar_qry_str[32];
    formatUs(sim.last_astar_query_us, astar_qry_str, sizeof(astar_qry_str));

    char plan_str[32];
    formatUs(sim.last_mtsp_plan_us, plan_str, sizeof(plan_str));

    double speedup_val = sim.last_paired_speedup.load();
    char speedup_str[32];
    if (speedup_val >= 10.0) {
        std::snprintf(speedup_str, sizeof(speedup_str), "%.0fx", speedup_val);
    } else if (speedup_val > 0.0) {
        std::snprintf(speedup_str, sizeof(speedup_str), "%.1fx", speedup_val);
    } else {
        std::snprintf(speedup_str, sizeof(speedup_str), "--");
    }

    double astar_speedup_val = sim.last_paired_astar_speedup.load();
    char astar_speedup_str[32];
    if (astar_speedup_val >= 10.0) {
        std::snprintf(astar_speedup_str, sizeof(astar_speedup_str), "%.0fx", astar_speedup_val);
    } else if (astar_speedup_val > 0.0) {
        std::snprintf(astar_speedup_str, sizeof(astar_speedup_str), "%.1fx", astar_speedup_val);
    } else {
        std::snprintf(astar_speedup_str, sizeof(astar_speedup_str), "--");
    }

    // Line 3: Algorithmic comparison rendered with distinct high-contrast color badges
    auto drawSegment = [&](const char* text, uint8_t r, uint8_t g, uint8_t b, int& cx) {
        drawString5x7(out_rgb, W, H, cx, y3, text, r, g, b, font_scale);
        cx += static_cast<int>(std::strlen(text)) * (6 * font_scale);
    };

    int cur_x = left_pad;

    char seg_htrip[64];
    std::snprintf(seg_htrip, sizeof(seg_htrip), "H-TRIP: %s", htrip_qry_str);
    drawSegment(seg_htrip, 50, 255, 170, cur_x); // Electric Emerald

    drawSegment(" | ", 110, 130, 155, cur_x);

    char seg_bfs[128];
    std::snprintf(seg_bfs, sizeof(seg_bfs), "BFS: %s (%s)", bfs_qry_str, speedup_str);
    drawSegment(seg_bfs, 255, 195, 60, cur_x); // Warm Amber

    if (sim.enable_astar_baseline && sim.last_astar_query_us.load() > 0.0) {
        drawSegment(" | ", 110, 130, 155, cur_x);

        char seg_astar[128];
        std::snprintf(seg_astar, sizeof(seg_astar), "A*: %s (%s)", astar_qry_str, astar_speedup_str);
        drawSegment(seg_astar, 70, 220, 255, cur_x); // Electric Cyan
    }

    drawSegment(" | ", 110, 130, 155, cur_x);

    char seg_plan[64];
    std::snprintf(seg_plan, sizeof(seg_plan), "Plan: %s", plan_str);
    drawSegment(seg_plan, 205, 215, 230, cur_x); // Soft Slate
}

bool writeP6Ppm(const std::string& path, int width, int height, const uint8_t* rgb) {
    FILE* fp = std::fopen(path.c_str(), "wb");
    if (!fp) return false;
    std::fprintf(fp, "P6\n%d %d\n255\n", width, height);
    std::fwrite(rgb, 1, static_cast<size_t>(width * height * 3), fp);
    std::fclose(fp);
    return true;
}

ArenaRecorder::~ArenaRecorder() {
    stopFfmpegStream();
    if (compile_thread.joinable()) {
        compile_thread.join();
    }
    cleanupPpmFrames();
}

void ArenaRecorder::cleanupPpmFrames() noexcept {
    if (!delete_ppm_after_compile) return;
    if (session_dir.empty()) return;
    if (!has_compiled_video) return;

    try {
        std::error_code ec;
        if (!std::filesystem::exists(session_dir, ec)) return;
        for (const auto& entry : std::filesystem::directory_iterator(session_dir, ec)) {
            if (entry.is_regular_file(ec) && entry.path().extension() == ".ppm") {
                std::filesystem::remove(entry.path(), ec);
            }
        }
        has_compiled_video = false;
    } catch (...) {}
}

void ArenaRecorder::startRecordingAt(const std::string& session_directory) {
    if (is_compiling.load()) {
        status_msg = "Cannot start new recording while ffmpeg compilation is running.";
        return;
    }

    waitForCompile();
    cleanupPpmFrames();

    std::error_code ec;
    std::filesystem::create_directories(session_directory, ec);

    session_dir = session_directory;
    frame_count = 0;
    has_compiled_video = false;
    last_video_path.clear();
    state = RecordingState::Recording;
    status_msg = "Recording active -> " + session_dir;
    if (stream_mp4) {
        if (!startFfmpegStream()) {
            status_msg = "ffmpeg stream failed to start; falling back to PPM frames";
            stream_mp4 = false;
        }
    }
}

void ArenaRecorder::startRecording(const std::string& base_dir) {
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);
    std::tm bt{};
    localtime_r(&in_time_t, &bt);

    char buf[64];
    std::strftime(buf, sizeof(buf), "rec_%Y%m%d_%H%M%S", &bt);
    startRecordingAt((std::filesystem::path(base_dir) / buf).string());
}

void ArenaRecorder::pauseRecording() noexcept {
    if (state == RecordingState::Recording) {
        state = RecordingState::Paused;
        status_msg = "Recording paused (" + std::to_string(frame_count) + " frames)";
    }
}

void ArenaRecorder::resumeRecording() noexcept {
    if (state == RecordingState::Paused) {
        state = RecordingState::Recording;
        status_msg = "Recording resumed -> " + session_dir;
    }
}

void ArenaRecorder::stopRecording() noexcept {
    if (state == RecordingState::Recording || state == RecordingState::Paused) {
        state = RecordingState::Idle;
        stopFfmpegStream();
        status_msg = "Recording finished: " + std::to_string(frame_count) + " frames in " + session_dir;
    }
}

void ArenaRecorder::captureTick(const htrip::MultiRobotSimulator& sim) {
    if (state != RecordingState::Recording || session_dir.empty()) return;

    int res = (resolution_px / 2) * 2;
    res = std::clamp(res, 384, 2048);

    rasterizeArena(sim, res, res, rgb_buffer);

    if (ffmpeg_pipe != nullptr) {
        const size_t nbytes = rgb_buffer.size();
        const size_t nw = std::fwrite(rgb_buffer.data(), 1, nbytes, ffmpeg_pipe);
        if (nw == nbytes) {
            frame_count++;
        }
        return;
    }

    char frame_name[64];
    std::snprintf(frame_name, sizeof(frame_name), "/frame_%06d.ppm", frame_count);
    std::string frame_path = session_dir + frame_name;

    if (writeP6Ppm(frame_path, res, res, rgb_buffer.data())) {
        frame_count++;
    }
}

bool ArenaRecorder::startFfmpegStream() {
    stopFfmpegStream();
    stream_res = (resolution_px / 2) * 2;
    stream_res = std::clamp(stream_res, 384, 2048);
    const std::string out_mp4 = session_dir + "/exploration_video_" + std::string(qualityTag()) + ".mp4";
    std::string cmd = "ffmpeg -y -f rawvideo -pix_fmt rgb24 -s " +
                      std::to_string(stream_res) + "x" + std::to_string(stream_res) +
                      " -r " + std::to_string(fps) +
                      " -i - -an -c:v libx264 -preset " + std::string(getPreset()) +
                      " -crf " + std::to_string(getCrf()) +
                      " -pix_fmt yuv420p -movflags +faststart " +
                      "\"" + out_mp4 + "\" > \"" + session_dir + "/ffmpeg.log\" 2>&1";
    ffmpeg_pipe = popen(cmd.c_str(), "w");
    if (ffmpeg_pipe == nullptr) return false;
    std::setvbuf(ffmpeg_pipe, nullptr, _IOFBF, 1 << 20);
    last_video_path = out_mp4;
    status_msg = "Streaming H.264 -> " + out_mp4;
    return true;
}

bool ArenaRecorder::stopFfmpegStream() noexcept {
    if (ffmpeg_pipe == nullptr) return has_compiled_video;
    const int rc = pclose(ffmpeg_pipe);
    ffmpeg_pipe = nullptr;
    if (rc == 0) {
        has_compiled_video = true;
        status_msg = "Video ready: " + last_video_path;
        return true;
    }
    status_msg = "ffmpeg stream error (code " + std::to_string(rc) + "). See " + session_dir + "/ffmpeg.log";
    return false;
}

bool ArenaRecorder::saveStill(const MultiRobotSimulator& sim, const std::string& path) {
    int res = (resolution_px / 2) * 2;
    res = std::clamp(res, 384, 2048);
    rasterizeArena(sim, res, res, rgb_buffer);
    std::error_code ec;
    std::filesystem::create_directories(std::filesystem::path(path).parent_path(), ec);
    return writeP6Ppm(path, res, res, rgb_buffer.data());
}

bool ArenaRecorder::saveLastFrame(const std::string& path) const {
    int res = (resolution_px / 2) * 2;
    res = std::clamp(res, 384, 2048);
    if (rgb_buffer.size() < static_cast<size_t>(res * res * 3)) return false;
    std::error_code ec;
    std::filesystem::create_directories(std::filesystem::path(path).parent_path(), ec);
    return writeP6Ppm(path, res, res, rgb_buffer.data());
}

void ArenaRecorder::waitForCompile() {
    if (compile_thread.joinable()) {
        compile_thread.join();
    }
}

namespace {

int runFfmpegEncode(const std::string& s_dir, int f_fps, int f_crf, const std::string& f_preset,
                    const std::string& q_tag, std::string& out_mp4, std::string& err_msg) {
    out_mp4 = s_dir + "/exploration_video_" + q_tag + ".mp4";
    std::string cmd = "ffmpeg -y -framerate " + std::to_string(f_fps) +
                      " -i \"" + s_dir + "/frame_%06d.ppm\" " +
                      "-c:v libx264 -preset " + f_preset +
                      " -crf " + std::to_string(f_crf) +
                      " -pix_fmt yuv420p " +
                      "\"" + out_mp4 + "\" > \"" + s_dir + "/ffmpeg.log\" 2>&1";
    int ret = std::system(cmd.c_str());
    if (ret != 0) {
        err_msg = "ffmpeg error (code " + std::to_string(ret) + "). See " + s_dir + "/ffmpeg.log";
    }
    return ret;
}

}  // namespace

void ArenaRecorder::compileVideoAsync() {
    if (is_compiling.load() || frame_count == 0 || session_dir.empty()) return;
    is_compiling.store(true);
    status_msg = "Compiling MP4 video with ffmpeg...";

    if (compile_thread.joinable()) {
        compile_thread.join();
    }

    std::string s_dir = session_dir;
    int f_count = frame_count;
    int f_fps = fps;
    int f_crf = getCrf();
    std::string f_preset = getPreset();
    VideoQuality f_quality = quality;

    compile_thread = std::thread([this, s_dir, f_count, f_fps, f_crf, f_preset, f_quality]() {
        (void)f_count;
        std::string q_tag = (f_quality == VideoQuality::High) ? "high" :
                            (f_quality == VideoQuality::Low) ? "low" : "medium";
        std::string out_mp4;
        std::string err;
        int ret = runFfmpegEncode(s_dir, f_fps, f_crf, f_preset, q_tag, out_mp4, err);
        if (ret == 0) {
            std::error_code ec;
            std::filesystem::copy_file(out_mp4, s_dir + "/exploration_video.mp4",
                                       std::filesystem::copy_options::overwrite_existing, ec);
            this->last_video_path = out_mp4;
            this->status_msg = "Video ready: " + out_mp4;
            this->has_compiled_video = true;
        } else {
            this->status_msg = err;
        }
        this->is_compiling.store(false);
    });
}

bool ArenaRecorder::compileVideoBlocking() {
    waitForCompile();
    if (has_compiled_video && !last_video_path.empty()) return true;
    if (frame_count == 0 || session_dir.empty()) return false;
    is_compiling.store(true);
    status_msg = "Compiling MP4 video with ffmpeg...";

    std::string q_tag = qualityTag();
    std::string out_mp4;
    std::string err;
    int ret = runFfmpegEncode(session_dir, fps, getCrf(), getPreset(), q_tag, out_mp4, err);
    if (ret == 0) {
        last_video_path = out_mp4;
        status_msg = "Video ready: " + out_mp4;
        has_compiled_video = true;
        if (delete_ppm_after_compile) {
            cleanupPpmFrames();
        }
        is_compiling.store(false);
        return true;
    }
    status_msg = err;
    is_compiling.store(false);
    return false;
}

} // namespace htrip
