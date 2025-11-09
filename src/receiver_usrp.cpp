// ============================================================================
// usrp_interface.cpp - Multi-rate capture with power-gated FIFO blocks
// ============================================================================
#include "receiver_usrp.hpp"
#include "filters.hpp" 
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <deque>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <cmath>
#include <csignal>
#include <chrono>
#include <numeric>

// Forward declarations - will be defined later
class BlockFIFO;
extern BlockFIFO g_fifo;
extern std::atomic<bool> g_stop;

// ---------------- USRP Interface implementation ----------------
ReceiverUSRP::ReceiverUSRP(const std::string& device_args,
                             double center_freq_hz,
                             double sample_rate_sps,
                             double gain_db,
                             const std::string& antenna,
                             const std::string& subdev)
: freq_(center_freq_hz), sps_(sample_rate_sps), gain_(gain_db), ant_(antenna), subdev_(subdev)
{
    std::cout << "Creating USRP interface..." << std::endl;
    usrp_ = uhd::usrp::multi_usrp::make(device_args);
    std::cout << "USRP created successfully" << std::endl;
    
    if (!subdev_.empty()) {
        std::cout << "Setting subdevice: " << subdev_ << std::endl;
        usrp_->set_rx_subdev_spec(subdev_);
    }
    
    std::cout << "Setting sample rate: " << sps_ << std::endl;
    usrp_->set_rx_rate(sps_);
    
    std::cout << "Setting frequency: " << freq_ << std::endl;
    usrp_->set_rx_freq(freq_);
    
    std::cout << "Setting gain: " << gain_ << std::endl;
    usrp_->set_rx_gain(gain_);
    
    std::cout << "Setting antenna: " << ant_ << std::endl;
    usrp_->set_rx_antenna(ant_);
    
    std::cout << "Creating streamer..." << std::endl;
    uhd::stream_args_t stream_args("fc32"); // complex float32
    rx_streamer_ = usrp_->get_rx_stream(stream_args);
    std::cout << "Streamer created successfully" << std::endl;
}

USRPInterface::~USRPInterface() { stop(); }

void USRPInterface::set_center_freq(double f) { freq_ = f; usrp_->set_rx_freq(f); }
void USRPInterface::set_sample_rate(double r) { sps_ = r; usrp_->set_rx_rate(r); }
void USRPInterface::set_gain(double g) { gain_ = g; usrp_->set_rx_gain(g); }

void USRPInterface::start() {
    std::cout << "USRPInterface::start() called" << std::endl;
    if (streaming_) {
        std::cout << "Already streaming, returning" << std::endl;
        return;
    }
    std::cout << "Creating stream command..." << std::endl;
    uhd::stream_cmd_t cmd(uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS);
    cmd.stream_now = true;
    cmd.num_samps = 0;
    std::cout << "Issuing stream command..." << std::endl;
    rx_streamer_->issue_stream_cmd(cmd);
    std::cout << "Stream command issued" << std::endl;
    streaming_ = true;
    std::cout << "Streaming flag set to true" << std::endl;
}
void USRPInterface::stop() {
    if (!streaming_) return;
    uhd::stream_cmd_t cmd(uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS);
    rx_streamer_->issue_stream_cmd(cmd);
    streaming_ = false;
}

size_t USRPInterface::read(std::complex<float>* out, size_t max_samples, double timeout_s) {
    size_t num_rx_samps = rx_streamer_->recv(out, max_samples, md_, timeout_s, false);
    
    // Handle error codes
    if (md_.error_code == uhd::rx_metadata_t::ERROR_CODE_TIMEOUT) {
        return 0;
    }
    if (md_.error_code == uhd::rx_metadata_t::ERROR_CODE_OVERFLOW) {
        return num_rx_samps;  // return what we got, caller will try again
    }
    if (md_.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE) {
        return 0;  // other errors
    }
    
    return num_rx_samps;
}

// Read exactly 10,000 samples either from one or multiple blocks.
bool USRPInterface::read_exact(std::complex<float>* out,
                               size_t n,
                               double total_timeout_s)
{
    std::complex<float>* ptr = out;
    size_t total = 0;

    auto deadline =
        std::chrono::steady_clock::now() +
        std::chrono::duration<double>(total_timeout_s);

    while (total < n) {
        double remaining =
            std::chrono::duration<double>(
                deadline - std::chrono::steady_clock::now()).count();
        if (remaining <= 0.0) return false;     // overall timeout

        size_t need = n - total;
        void* buffs[] = { ptr };
        size_t got = rx_streamer_->recv(
            buffs, need, md_, std::min(remaining, 0.2), /*one_packet=*/false);

        // handle metadata
        if (md_.error_code == uhd::rx_metadata_t::ERROR_CODE_TIMEOUT) continue;
        if (md_.error_code == uhd::rx_metadata_t::ERROR_CODE_OVERFLOW) continue;
        if (md_.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE) return false;

        if (got == 0) continue;
        total += got;
        ptr += got;  // advance write pointer
    }
    return true;                                // filled n samples
}


// ---- Baseband low-pass (passband = fs_out/4) for rational U/D ----
static std::vector<std::complex<float>>
design_baseband_lp_taps_for_resampler(int D, int taps_len = 161)
{
    // normalized cutoff at the upsampled rate: fc_norm = 1/(2D)
    // use a slight safety margin (e.g., 0.9) to keep ripple/alias small
    double fc_norm = 0.9 * (1.0 / (2.0 * D));
    const int N = (taps_len % 2 == 1) ? taps_len : (taps_len + 1); // force odd length
    const int mid = (N - 1) / 2;
    std::vector<std::complex<float>> h(N);
    double sum = 0.0;
    for (int n = 0; n < N; ++n) {
        int m = n - mid;
        double sinc = (m == 0) ? 1.0 : std::sin(2.0 * M_PI * fc_norm * m) / (M_PI * m);
        double win  = 0.54 - 0.46 * std::cos(2.0 * M_PI * n / (N - 1)); // Hamming
        double v = sinc * win;
        h[n] = { static_cast<float>(v), 0.0f };
        sum += v;
    }
    // Normalize DC gain to 1.0
    for (auto &c : h) c.real(static_cast<float>(c.real() / sum));
    return h;
}

// ---------------- Power gating and FIFO ----------------
struct PowerGate {
    float th_hi{1.0f};  // enter
    float th_lo{0.5f};  // exit
    bool  on{false};
    bool update(float p) {
        if (!on && p >= th_hi) on = true;
        else if (on && p <= th_lo) on = false;
        return on;
    }
};

inline float inst_power(const std::complex<float>& s) { return std::norm(s); }

class BlockFIFO {
    std::queue<std::vector<std::complex<float>>> q_;
    std::mutex m_;
    std::condition_variable cv_;
public:
    void push(std::vector<std::complex<float>>&& blk) {
        std::lock_guard<std::mutex> lk(m_);
        q_.push(std::move(blk));
        cv_.notify_one();
    }
    bool pop(std::vector<std::complex<float>>& out) {
        std::unique_lock<std::mutex> lk(m_);
        // Use wait_for with timeout and check g_stop
        if (!cv_.wait_for(lk, std::chrono::milliseconds(200), [&]{ return !q_.empty() || g_stop; })) {
            return false;  // timeout
        }
        if (g_stop && q_.empty()) {
            return false;  // stopping and no data
        }
        out = std::move(q_.front());
        q_.pop();
        return true;
    }
    void notify_stop() {
        cv_.notify_all();  // Wake up any waiting threads
    }
};

BlockFIFO g_fifo;  // Global instance

// ---------------- Signal handling ----------------
std::atomic<bool> g_stop{false};
static void on_signal(int){ 
    g_stop = true; 
    g_fifo.notify_stop();  // Wake up consumer thread
    std::cout << "\nStopping... (please wait)\n" << std::flush;
}


// ---- Cutoff frequency in upsampled domain (fc_up = 0.25 / D) ----
static std::vector<std::complex<float>>
design_lp_taps_updomain_for_UD(int D, int taps_len = 201, double margin = 0.0)
{
    // Normalized cyclic cutoff (cycles/sample at fs_up)
    double fc_norm_up = (0.25 / static_cast<double>(D)) * (1.0 - margin);

    const int N   = (taps_len % 2) ? taps_len : (taps_len + 1);
    const int mid = (N - 1) / 2;

    std::vector<std::complex<float>> h(N);
    double sum = 0.0;

    for (int n = 0; n < N; ++n) {
        int m = n - mid;
        double sinc = (m == 0)
            ? 1.0
            : std::sin(2.0 * M_PI * fc_norm_up * m) / (M_PI * m);
        double win = 0.54 - 0.46 * std::cos(2.0 * M_PI * n / (N - 1));
        double v = sinc * win;
        h[n] = { static_cast<float>(v), 0.0f };
        sum += v;
    }

    for (auto &c : h) c.real(static_cast<float>(c.real() / sum));
    return h;
}

// --- Smoothed (first-order IIR) power detector: y[n] = a*y[n-1] + (1-a)*|x[n]|^2
struct PowerIIR {
    float a;          // 0 <= a < 1
    float y = 0.0f;   // state
    explicit PowerIIR(float alpha) : a(alpha) {}
    inline float update(const std::complex<float>& s) {
        // |x|^2 = I^2 + Q^2
        float p = std::norm(s);
        y = a * y + (1.0f - a) * p;
        return y;
    }
};

// ============================================================================
// Producer & Consumer with Lab 3 Requirements
// ============================================================================
// RX SAMPLING RATE CHOICE: 800 kHz (matches TX symbol rate)
// 
// Rationale:
// - TX: 800 kSps symbol rate, RRC β=0.25 → 1 MHz occupied bandwidth
// - Nyquist: Need fs >= 1 MHz, but after baseband filtering
// - Raw USRP rate: 1 MHz, then resample to 800 kHz
// - Why 800 kHz?
//   1. Matches TX symbol rate (convenient for processing)
//   2. Satisfies Nyquist after LP filtering (passband = fs_out/4 = 200 kHz)
//   3. Reduces computational load vs. 1 MHz
//   4. Adequate for energy detection (not demodulating)
//
// FILTERING CHOICE: Yes - Rational resampling with LP filter
//
// Rationale:
// - Raw USRP at 1 MHz needs resampling to 800 kHz
// - Polyphase LP filter (201 taps) used during resampling
// - Passband: ~200 kHz (fs_out/4)
// - Benefits:
//   1. Anti-aliasing for resampling
//   2. Rejects out-of-band noise
//   3. Improves SNR for energy detection
// - Alternative (no filtering) considered but rejected:
//   Would need fs_out = fs_raw (no resampling)
//   More samples to process, no SNR benefit
// ============================================================================

void producer(USRPInterface& radio,
                     double fs_target,
                     float th_lin,  // linear power threshold
                     float alpha)   // IIR smoothing parameter, 0<=alpha<1
{
    std::cout << "Producer function started" << std::endl;
    const size_t RAW_READ = 10000;               // strict USRP read size
    const double fs_in = radio.actual_sample_rate();
    std::cout << "Got sample rate: " << fs_in << std::endl;

    // Rational resampling ratio U/D ~= fs_out / fs_in (reduced)
    std::cout << "Converting to integers..." << std::endl;
    auto to_int = [](double v){ return static_cast<int>(std::llround(v)); };
    int in_i  = to_int(fs_in);
    int out_i = to_int(fs_target);
    std::cout << "in_i: " << in_i << ", out_i: " << out_i << std::endl;
    
    std::cout << "Calculating GCD..." << std::endl;
    int g = std::gcd(out_i, in_i);
    std::cout << "GCD: " << g << std::endl;
    
    std::cout << "Calculating U and D..." << std::endl;
    int U = out_i / g;
    int D = in_i / g;
    std::cout << "U: " << U << ", D: " << D << std::endl;

    // Design baseband LP taps for this D (passband = fs_out/4)
    // You can tune tap length (e.g., 161, 201, 241) for stopband/CPU tradeoff.
    std::cout << "Designing filter taps..." << std::endl;
    auto taps = design_lp_taps_updomain_for_UD(D, /*taps_len*/ 201, /*margin*/ 0.0);
    std::cout << "Filter taps designed, size: " << taps.size() << std::endl;

    // Choose a processing block size for FilterPolyphase: must be divisible by D
    const size_t XLEN = RAW_READ - (RAW_READ % D);  // e.g., RAW_READ if divisible
    std::cout << "XLEN: " << XLEN << std::endl;
    std::cout << "Creating FilterPolyphase..." << std::endl;
    FilterPolyphase mr(U, D, static_cast<int>(XLEN), static_cast<int>(taps.size()), taps.data(), /*threads*/ 1);
    std::cout << "FilterPolyphase created" << std::endl;
    mr.set_head(true); // first block is "head"
    std::cout << "FilterPolyphase head set" << std::endl;

    std::cout << "Creating PowerIIR..." << std::endl;
    PowerIIR pwr(alpha);
    std::cout << "PowerIIR created" << std::endl;
    
    // ========================================================================
    // PACKET CAPTURE SIZE CALCULATION (Lab 3 Requirements)
    // ========================================================================
    // TX parameters:
    //   - 1000 bits/packet
    //   - 800 kSps symbol rate
    //   - 1.25 samples/symbol at 1 MHz
    //   - After resampling to 800 kHz: 1000 samples/packet
    //   - RRC filter transients: ~120 samples
    //   - Total TX packet: ~1120 samples at 800 kHz
    //
    // Capture strategy (to capture WHOLE packet with sync margin):
    //   - Pre-trigger: 100 samples (before detection trigger)
    //   - Packet data: 1500 samples (includes TX packet + margin)
    //   - Post-trigger: 100 samples (after packet end)
    //   - TOTAL: 1700 samples
    //
    // Rationale:
    //   - Pre-trigger captures packet start (IIR detector has delay)
    //   - Margin handles timing uncertainty and RRC transients
    //   - Post-trigger provides clean end boundary
    //   - Total provides ample margin for future synchronization
    // ========================================================================
    const size_t PRE_TRIGGER_SAMPLES = 100;
    const size_t PACKET_DATA_SAMPLES = 1500;
    const size_t POST_TRIGGER_SAMPLES = 100;
    const size_t TOTAL_CAPTURE_SIZE = PRE_TRIGGER_SAMPLES + 
                                      PACKET_DATA_SAMPLES + 
                                      POST_TRIGGER_SAMPLES;  // = 1700
    
    std::cout << "Creating vectors..." << std::endl;
    std::cout << "Packet capture size: " << TOTAL_CAPTURE_SIZE << " samples\n";
    std::cout << "  Pre-trigger:  " << PRE_TRIGGER_SAMPLES << " samples\n";
    std::cout << "  Packet data:  " << PACKET_DATA_SAMPLES << " samples\n";
    std::cout << "  Post-trigger: " << POST_TRIGGER_SAMPLES << " samples\n";
    
    // Pre-trigger circular buffer
    std::deque<std::complex<float>> pre_trigger_buffer;
    
    std::vector<std::complex<float>> pack;
    pack.reserve(TOTAL_CAPTURE_SIZE);
    bool capturing = false;
    size_t cap_count = 0;

    // staging FIFO for raw samples (to satisfy XLEN divisibility each filter call)
    std::deque<std::complex<float>> stage;
    std::vector<std::complex<float>> in_blk(XLEN);
    // Worst-case output per call: floor(XLEN/D)*U  (FilterPolyphase contract)
    std::vector<std::complex<float>> out_blk((XLEN / D) * U + 8);
    std::cout << "Vectors created" << std::endl;

    std::cout << "Starting radio..." << std::endl;
    radio.start();
    std::cout << "Radio started" << std::endl;
    while (!g_stop) {
        // 1) Read samples (may get partial reads)
        std::vector<std::complex<float>> raw(RAW_READ);
        size_t got = radio.read(raw.data(), RAW_READ, /*timeout_s=*/0.2);
        if (got == 0) continue; // timeout, try again
        
        raw.resize(got); // trim to actual size received
            
        // 2) Append to staging buffer
        stage.insert(stage.end(), raw.begin(), raw.end());

        // 3) While we have at least XLEN samples (divisible by D), filter them
        while (stage.size() >= XLEN) {
            // copy XLEN samples into a contiguous input block
            std::copy_n(stage.begin(), XLEN, in_blk.begin());
            stage.erase(stage.begin(), stage.begin() + XLEN);

            // continuous streaming after the first call
            mr.set_head(false);
            int nout = mr.filter(in_blk.data(), out_blk.data());

            // ========================================================
            // PACKET CAPTURE WITH PRE-TRIGGER (Lab 3 Requirements)
            // ========================================================
            // Captures WHOLE packet (1700 samples total):
            //   - 100 pre-trigger (before detection)
            //   - 1500 packet data (TX packet + margin)
            //   - 100 post-trigger (after packet end)
            // Handles block boundaries via deque accumulation
            // ========================================================
            for (int i = 0; i < nout; ++i) {
                const auto& s = out_blk[i];
                float yiir = pwr.update(s);       // smoothed power
                
                // Maintain pre-trigger circular buffer
                pre_trigger_buffer.push_back(s);
                if (pre_trigger_buffer.size() > PRE_TRIGGER_SAMPLES) {
                    pre_trigger_buffer.pop_front();
                }
                
                if (!capturing && yiir > th_lin) {
                    // Trigger detected! Start capturing
                    capturing = true;
                    cap_count = 0;
                    pack.clear();
                    pack.reserve(TOTAL_CAPTURE_SIZE);
                    
                    // Copy pre-trigger samples (captures packet start)
                    for (const auto& pre_samp : pre_trigger_buffer) {
                        pack.push_back(pre_samp);
                        cap_count++;
                    }
                }
                
                if (capturing) {
                    // Only add current sample if not already added from pre-trigger
                    if (cap_count >= PRE_TRIGGER_SAMPLES) {
                        pack.push_back(s);
                        ++cap_count;
                    }
                    
                    if (cap_count >= TOTAL_CAPTURE_SIZE) {
                        // ====================================================
                        // AGC: Normalize RMS magnitude to 1
                        // ====================================================
                        // For BPSK symbols with magnitude 1, RMS should be 1
                        // Compute RMS: sqrt(mean(|sample|^2))
                        // Scale factor: 1.0 / RMS
                        
                        double sum_power = 0.0;
                        for (const auto& sample : pack) {
                            sum_power += std::norm(sample);  // |sample|^2
                        }
                        double mean_power = sum_power / pack.size();
                        float rms = std::sqrt(static_cast<float>(mean_power));
                        
                        // Avoid division by zero
                        if (rms > 1e-6f) {
                            float scale = 1.0f / rms;
                            
                            // Apply AGC: Scale all samples
                            for (auto& sample : pack) {
                                sample *= scale;
                            }
                            
                            std::cout << "AGC applied: RMS=" << std::fixed 
                                      << std::setprecision(6) << rms 
                                      << ", scale=" << scale << "\n";
                        } else {
                            std::cout << "Warning: RMS too small, skipping AGC\n";
                        }
                        
                        // Push normalized packet to FIFO
                        g_fifo.push(std::move(pack));
                        pack.clear();
                        pack.reserve(TOTAL_CAPTURE_SIZE);
                        capturing = false;
                    } 
                }
                else {
                    // idle
                }
            }
        }
    }
    radio.stop();
}

static void consumer(const std::string& save_file = "") {
    std::ofstream ofs;
    bool saving = false;
    
    if (!save_file.empty()) {
        ofs.open(save_file, std::ios::binary);
        if (ofs) {
            saving = true;
            std::cout << "Saving detected samples to: " << save_file << std::endl;
        } else {
            std::cerr << "Warning: Could not open " << save_file << " for writing" << std::endl;
        }
    }
    
    // ========================================================================
    // PACKET COUNTING FOR ENERGY DETECTOR VERIFICATION
    // ========================================================================
    // Purpose: Monitor packet detection rate to verify energy detector
    // - TX sends ~1 packet/second → expect ~10 packets per 10 seconds
    // - Too few: detector threshold too high (missing packets)
    // - Too many: detector threshold too low (false alarms)
    // ========================================================================
    
    std::vector<std::complex<float>> blk;
    int total_packets = 0;
    int packets_in_interval = 0;
    
    // Start timer for 10-second intervals
    auto interval_start = std::chrono::steady_clock::now();
    const std::chrono::seconds REPORT_INTERVAL(10);
    
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "Energy Detector Performance Monitor\n";
    std::cout << std::string(70, '=') << "\n";
    std::cout << "Counting packets every 10 seconds...\n";
    std::cout << "Expected: ~10 packets/10s (1 packet/second from TX)\n";
    std::cout << std::string(70, '=') << "\n\n";
    
    while (!g_stop) {
        if (!g_fifo.pop(blk)) continue;
        
        total_packets++;
        packets_in_interval++;
        
        // Calculate and print statistics after AGC
        double sum_power = 0.0;
        for (auto &s : blk) {
            sum_power += std::norm(s);
        }
        double mean_power = sum_power / blk.size();
        double rms = std::sqrt(mean_power);
        
        std::cout << "Packet #" << total_packets << " (AGC normalized):\n";
        std::cout << "  Samples: " << blk.size() << "\n";
        std::cout << "  RMS magnitude: " << std::fixed << std::setprecision(6) << rms 
                  << " (should be ~1.0)\n";
        std::cout << "  Avg power: " << mean_power << " (should be ~1.0)\n";
        
        // Save to file if enabled
        if (saving) {
            ofs.write(reinterpret_cast<const char*>(blk.data()),
                     blk.size() * sizeof(std::complex<float>));
        }
        
        // ====================================================================
        // CHECK IF 10 SECONDS HAVE ELAPSED
        // ====================================================================
        auto now = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
            now - interval_start);
        
        if (elapsed >= REPORT_INTERVAL) {
            // ================================================================
            // PRINT 10-SECOND STATISTICS
            // ================================================================
            std::cout << "\n" << std::string(70, '=') << "\n";
            std::cout << "╔══════════════════════════════════════════════════════════════════╗\n";
            std::cout << "║          ENERGY DETECTOR PERFORMANCE (10-second report)          ║\n";
            std::cout << "╚══════════════════════════════════════════════════════════════════╝\n";
            std::cout << std::string(70, '=') << "\n";
            
            // Calculate actual elapsed time (may be slightly > 10s)
            double actual_elapsed = elapsed.count();
            double packets_per_second = packets_in_interval / actual_elapsed;
            
            std::cout << "  Elapsed time:          " << std::fixed << std::setprecision(2) 
                      << actual_elapsed << " seconds\n";
            std::cout << "  Packets captured:      " << packets_in_interval << "\n";
            std::cout << "  Packet rate:           " << std::setprecision(3) 
                      << packets_per_second << " packets/second\n";
            std::cout << "  Total packets so far:  " << total_packets << "\n";
            
            // ================================================================
            // ANALYSIS: Energy Detector Performance Assessment
            // ================================================================
            std::cout << "\n--- Energy Detector Assessment ---\n";
            
            // Expected: ~10 packets per 10 seconds (1 packet/second from TX)
            const double EXPECTED_RATE = 1.0;  // packets per second
            const double TOLERANCE = 0.2;       // ±20% tolerance
            
            if (std::abs(packets_per_second - EXPECTED_RATE) <= TOLERANCE) {
                std::cout << "✓ GOOD: Detection rate matches TX rate (~1 packet/s)\n";
                std::cout << "  Detector threshold appears properly set.\n";
            } else if (packets_per_second < EXPECTED_RATE - TOLERANCE) {
                std::cout << "⚠ WARNING: Detection rate LOW (" << packets_per_second 
                          << " < " << (EXPECTED_RATE - TOLERANCE) << " packets/s)\n";
                std::cout << "  Possible causes:\n";
                std::cout << "    • Detector threshold TOO HIGH (missing packets)\n";
                std::cout << "    • RX gain too low (weak signal)\n";
                std::cout << "    • TX not transmitting continuously\n";
                std::cout << "  Suggested fix: Lower threshold (--th_db -25 or -30)\n";
            } else if (packets_per_second > EXPECTED_RATE + TOLERANCE) {
                std::cout << "⚠ WARNING: Detection rate HIGH (" << packets_per_second 
                          << " > " << (EXPECTED_RATE + TOLERANCE) << " packets/s)\n";
                std::cout << "  Possible causes:\n";
                std::cout << "    • Detector threshold TOO LOW (false alarms)\n";
                std::cout << "    • Interference from other sources\n";
                std::cout << "    • Noise triggering detector\n";
                std::cout << "  Suggested fix: Raise threshold (--th_db -15 or -10)\n";
            }
            
            // Additional diagnostics
            if (packets_in_interval == 0) {
                std::cout << "\n✗ ERROR: NO PACKETS DETECTED in 10 seconds!\n";
                std::cout << "  Check:\n";
                std::cout << "    • Is transmitter running?\n";
                std::cout << "    • Is RX frequency correct (--fc)?\n";
                std::cout << "    • Is RX gain sufficient (try --gain 50)?\n";
                std::cout << "    • Is threshold too high (try --th_db -30)?\n";
            }
            
            std::cout << std::string(70, '=') << "\n\n";
            
            // Reset interval counter and timer
            packets_in_interval = 0;
            interval_start = now;
        }
    }
    
    // ========================================================================
    // FINAL STATISTICS
    // ========================================================================
    if (saving) {
        ofs.close();
        std::cout << "\nSaved " << total_packets << " packets to " << save_file << std::endl;
    }
    
    // Calculate total runtime
    auto final_time = std::chrono::steady_clock::now();
    auto total_runtime = std::chrono::duration_cast<std::chrono::seconds>(
        final_time - interval_start + REPORT_INTERVAL * (total_packets / 10));
    
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "FINAL STATISTICS\n";
    std::cout << std::string(70, '=') << "\n";
    std::cout << "Total packets captured: " << total_packets << "\n";
    if (total_runtime.count() > 0) {
        double avg_rate = static_cast<double>(total_packets) / total_runtime.count();
        std::cout << "Average packet rate:    " << std::fixed << std::setprecision(3) 
                  << avg_rate << " packets/second\n";
    }
    std::cout << std::string(70, '=') << "\n";
}

// ---------------- Main demo ----------------
static void usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " --dev <args> --fc <Hz> --fs_raw <sps> --fs_out <sps> --gain <dB> [--th_db <dB>] [--alpha <0-1>] [--save <file>]\n"
        << "Example: " << prog << " --dev \"type=usrp\" --fc 2.427e9 --fs_raw 1e6 --fs_out 0.8e6 --gain 40 --th_db -20 --alpha 0.5 --save samples_iq.bin\n"
        << "\nOptions:\n"
        << "  --dev <args>    : USRP device arguments (default: \"type=usrp\")\n"
        << "  --fc <Hz>       : Center frequency in Hz (default: 2.4e9)\n"
        << "  --fs_raw <sps>  : Raw USRP sampling rate in Hz (default: 1e6)\n"
        << "  --fs_out <sps>  : Output sampling rate after filtering in Hz (default: 800e3)\n"
        << "  --gain <dB>     : RX gain in dB (default: 40.0)\n"
        << "  --th_db <dB>    : Power detection threshold in dB (default: -20.0)\n"
        << "  --alpha <0-1>   : IIR smoothing parameter, 0<=alpha<1 (default: 0.5)\n"
        << "  --save <file>   : Save detected samples to binary file (optional)\n"
        << "\nLab 3 Features:\n"
        << "  • RX rate: 800 kHz (matches TX symbol rate)\n"
        << "  • Filtering: Polyphase LP filter during resampling\n"
        << "  • Packet capture: 1700 samples (100 pre + 1500 data + 100 post)\n"
        << "  • Pre-trigger buffer: Captures packet start despite detector delay\n"
        << "  • Block-boundary handling: Deque accumulation across USRP blocks\n"
        << "  • AGC: Normalizes RMS magnitude of each packet to 1 (BPSK magnitude)\n"
        << "  • Performance monitoring: Packet counting every 10 seconds\n";
}

int main(int argc, char** argv) {
    std::string dev = "type=usrp";
    // double fc = 915e6;  
    double fc = 2.4e9;
    double fs_raw = 1e6;
    double fs_out = 800e3;  // Lab 3: Match TX symbol rate (800 kSps)
    double gain = 40.0;
    double th_db = -20.0;
    double hyst_db = 3.0;
    double alpha = 0.5;    // IIR smoothing, 0<=alpha<1 (higher = smoother/longer memory)
    std::string save_file = "";  // Optional file to save detected samples

    for (int i=1; i<argc; ++i) {
        std::string a = argv[i];
        auto need = [&](int more){ if (i+more >= argc) { usage(argv[0]); return false; } return true; };
        if (a == "--dev" && need(1)) dev = argv[++i];
        else if (a == "--fc" && need(1)) fc = std::stod(argv[++i]);
        else if (a == "--fs_raw" && need(1)) fs_raw = std::stod(argv[++i]);
        else if (a == "--fs_out" && need(1)) fs_out = std::stod(argv[++i]);
        else if (a == "--gain" && need(1)) gain = std::stod(argv[++i]);
        else if (a == "--th_db" && need(1)) th_db = std::stod(argv[++i]);
        else if (a == "--hyst_db" && need(1)) hyst_db = std::stod(argv[++i]);
        else if (a == "--alpha" && need(1)) alpha   = std::stod(argv[++i]);
        else if (a == "--save" && need(1)) save_file = argv[++i];
        else if (a == "--help") { usage(argv[0]); return 0; }
    }

    try {
        std::cout << "Creating USRPInterface..." << std::endl;
        USRPInterface radio(dev, fc, fs_raw, gain, "TX/RX", "A:0");
        std::cout << "USRPInterface created successfully" << std::endl;
        
        std::cout << "Setting up signal handlers..." << std::endl;
        std::signal(SIGINT, on_signal);
        std::signal(SIGTERM, on_signal);
        std::cout << "Signal handlers set" << std::endl;

        // Convert dB thresholds to linear power
        // instantaneous power is |IQ|^2 with unity expected around 1 for full-scale
        std::cout << "Converting threshold..." << std::endl;
        float th_hi = static_cast<float>(std::pow(10.0, th_db/10.0));
        std::cout << "Threshold converted: " << th_hi << std::endl;
        
        std::cout << "Starting producer thread..." << std::endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        std::thread t_prod(producer, std::ref(radio), fs_out,
                   /*th_lin=*/th_hi,
                   /*alpha=*/static_cast<float>(alpha));
        std::cout << "Producer thread started" << std::endl;
        
        std::cout << "Starting consumer thread..." << std::endl;
        std::thread t_cons(consumer, save_file);
        std::cout << "Consumer thread started" << std::endl;

        t_prod.join();
        g_stop = true; // ensure consumer exits if producer ended
        g_fifo.notify_stop();  // Wake up consumer
        t_cons.join();
    } catch (const std::exception& ex) {
        std::cerr << "Fatal: " << ex.what() << std::endl;
        return 1;
    }
    return 0;
}
