#pragma once

#include "robot_simulator.hpp"
#include <string>
#include <vector>
#include <cstdint>
#include <atomic>
#include <cstdio>
#include <filesystem>
#include <thread>

namespace htrip {

enum class RecordingState {
    Idle,
    Recording,
    Paused
};

enum class VideoQuality : int {
    Low = 0,
    Medium = 1,
    High = 2
};

/**
 * @brief High-performance non-blocking arena video recorder and snapshot engine.
 *
 * Designed to capture live simulation experiments for publication videos, supplementary material,
 * and visual performance auditing.
 *
 * Supports two recording modes:
 * 1. Disk-Backed Binary PPM (P6) Streaming:
 *    Writes raw uncompressed PPM frames at each simulation tick with negligible CPU overhead (< 0.5 ms/frame).
 *    Once the mission finishes, ffmpeg is invoked asynchronously in the background to encode an H.264 MP4 video.
 * 2. Live Pipe Streaming (`stream_mp4 = true`):
 *    Spawns an ffmpeg subprocess connected via stdin pipe (`popen`), streaming RGB frames directly
 *    without intermediate disk writes.
 */
struct ArenaRecorder {
    RecordingState state = RecordingState::Idle;    ///< Current recording lifecycle state.
    VideoQuality quality = VideoQuality::Medium;    ///< Compression quality and resolution tier.
    std::string session_dir;                        ///< Filesystem directory storing frame dumps and output MP4.
    std::string last_video_path;                    ///< Absolute path to the most recently compiled MP4 video.
    std::string status_msg = "Ready";               ///< Human-readable status message for GUI telemetry.
    int frame_count = 0;                            ///< Monotonic counter of captured video frames.
    int fps = 20;                                   ///< Target playback framerate (frames per second).
    int resolution_px = 768;                        ///< Video square resolution (e.g. 512, 768, or 1024 px).
    bool delete_ppm_after_compile = true;           ///< Cleans up raw PPM frames after successful MP4 compilation.
    bool has_compiled_video = false;                ///< Set to true once an MP4 video has been successfully produced.
    bool stream_mp4 = false;                        ///< Live ffmpeg stdin pipe streaming mode.
    std::atomic<bool> is_compiling{false};          ///< Thread-safe flag indicating background ffmpeg encoding is active.

    std::vector<uint8_t> rgb_buffer;                ///< Persistent reusable RGB pixel buffer.
    std::thread compile_thread;                     ///< Dedicated worker thread for non-blocking ffmpeg compilation.
    std::FILE* ffmpeg_pipe = nullptr;               ///< Pipe handle for live ffmpeg streaming.
    int stream_res = 0;                             ///< Stream resolution width/height in pixels.

    ~ArenaRecorder();

    /**
     * @brief Retrieves the Constant Rate Factor (CRF) for x264 video compression.
     */
    int getCrf() const noexcept {
        switch (quality) {
            case VideoQuality::Low: return 24;
            case VideoQuality::Medium: return 18;
            case VideoQuality::High: return 12;
        }
        return 18;
    }

    /**
     * @brief Retrieves the ffmpeg x264 encoding speed preset.
     */
    const char* getPreset() const noexcept {
        switch (quality) {
            case VideoQuality::Low: return "fast";
            case VideoQuality::Medium: return "medium";
            case VideoQuality::High: return "slow";
        }
        return "medium";
    }

    /**
     * @brief Configures video quality and corresponding square pixel resolution.
     */
    void setQuality(VideoQuality q) noexcept {
        quality = q;
        switch (q) {
            case VideoQuality::Low:
                resolution_px = 512;
                break;
            case VideoQuality::Medium:
                resolution_px = 768;
                break;
            case VideoQuality::High:
                resolution_px = 1024;
                break;
        }
    }

    /**
     * @brief Returns a string tag identifying the active quality tier.
     */
    [[nodiscard]] const char* qualityTag() const noexcept {
        switch (quality) {
            case VideoQuality::Low: return "low";
            case VideoQuality::High: return "high";
            case VideoQuality::Medium: return "medium";
        }
        return "medium";
    }

    /**
     * @brief Initiates a new recording session under a timestamped directory inside base_dir.
     */
    void startRecording(const std::string& base_dir = "recordings");

    /**
     * @brief Initiates a new recording session in a specific target directory.
     */
    void startRecordingAt(const std::string& session_directory);

    /**
     * @brief Pauses frame capture without closing the active recording session.
     */
    void pauseRecording() noexcept;

    /**
     * @brief Resumes frame capture into the active recording session.
     */
    void resumeRecording() noexcept;

    /**
     * @brief Finalizes frame capture and closes active ffmpeg pipes.
     */
    void stopRecording() noexcept;

    /**
     * @brief Rasterizes the current simulation state and writes a video frame.
     */
    void captureTick(const MultiRobotSimulator& sim);

    /**
     * @brief Rasterizes the current simulation state and exports a single high-resolution still image.
     */
    bool saveStill(const MultiRobotSimulator& sim, const std::string& path);

    /**
     * @brief Saves the most recently captured frame buffer to a file.
     */
    bool saveLastFrame(const std::string& path) const;

    /**
     * @brief Spawns a background thread to compile recorded PPM frames into an H.264 MP4 video via ffmpeg.
     */
    void compileVideoAsync();

    /**
     * @brief Compiles recorded PPM frames into an MP4 video synchronously on the calling thread.
     */
    bool compileVideoBlocking();

    /**
     * @brief Blocks until active background ffmpeg compilation completes.
     */
    void waitForCompile();

    /**
     * @brief Deletes temporary PPM frame files after successful MP4 encoding.
     */
    void cleanupPpmFrames() noexcept;

    /**
     * @brief Opens an ffmpeg subprocess pipe for live H.264 stream encoding.
     */
    bool startFfmpegStream();

    /**
     * @brief Flushes and closes the live ffmpeg subprocess pipe.
     */
    bool stopFfmpegStream() noexcept;
};

/**
 * @brief Rasterize continuous multi-robot simulator arena into an RGB image buffer.
 */
void rasterizeArena(const htrip::MultiRobotSimulator& sim, int width, int height, std::vector<uint8_t>& out_rgb);

/**
 * @brief Fast binary P6 PPM image writer.
 */
bool writeP6Ppm(const std::string& path, int width, int height, const uint8_t* rgb);

} // namespace htrip
