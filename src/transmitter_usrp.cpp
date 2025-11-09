
// ============================================================================
// transmitter_alice.cpp - BPSK/GMSK Transmitter with Improved RRC
// ============================================================================


#include "transmitter_usrp.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <chrono>
#include <csignal>
#include <cmath>
#include <thread>
#include <algorithm>
#include <stdexcept>

// ============================================================================
// Global stop flag and signal handler
// ============================================================================
std::atomic<bool> g_stop{false};

void signal_handler(int) {
    g_stop = true;
    std::cout << "\nStopping transmitter...\n" << std::flush;
}

// ============================================================================
// NEW: Truncated RRC Pulse Generation Function
// ============================================================================
// PURPOSE: Generate unit-energy truncated RRC impulse response for pulse shaping
// EXAMPLE:
//   For 1 MHz sample rate, 800 kHz symbol rate:
//   - samples_per_symbol = 1.25 = 5/4
//   - U = 5, D = 4
//   - beta = (5-4)/4 = 0.25 (25% excess bandwidth)
//   - Occupied BW = 800kHz × (1+0.25) = 1 MHz
// ============================================================================
void rrc_pulse(std::complex<float>* h, int len, int U, int D, float beta)
{
    // CENTER TAP (n=0): Special case using L'Hôpital's rule
    // h[len] corresponds to n=0 in the impulse response
    h[len] = 1.0f - beta + 4.0f * beta / M_PI;
    
    // Initialize energy accumulator for normalization
    float scale = std::norm(h[len]);  // |h[len]|^2
    
    // COMPUTE OFF-CENTER TAPS (n=1 to len)
    for (int n = 1; n <= len; n++) {
        // SPECIAL CASE: Check if n is at the zero-crossing of denominator
        // This occurs at n = U/(4*beta)
        // Use tolerance to handle floating-point comparison
        if (std::abs(n - U / beta / 4.0f) < 1e-6f) {
            // Apply special formula for this singularity
            h[len + n] = beta / std::sqrt(2.0f) * 
                        ((1.0f + 2.0f / M_PI) * std::sin(M_PI / 4.0f / beta) +
                         (1.0f - 2.0f / M_PI) * std::cos(M_PI / 4.0f / beta));
        } 
        // GENERAL CASE: Standard RRC formula
        else {
            // Numerator: sin term + cos term with frequency modulation
            float numerator = std::sin(n * M_PI * (1.0f - beta) / U) + 
                             4.0f * n * beta / U * std::cos(n * M_PI * (1.0f + beta) / U);
            
            // Denominator: includes zero-crossing term (1 - 16n²β²/U²)
            float denominator = n * M_PI / U * (1.0f - 16.0f * n * n * beta * beta / U / U);
            
            // Complete formula with normalization factor U/π
            h[len + n] = numerator * U / denominator;
        }
        
        // SYMMETRY: RRC impulse response is symmetric
        h[len - n] = h[len + n];
        
        // ACCUMULATE ENERGY: Count both positive and negative taps
        scale += 2.0f * std::norm(h[len + n]);
    }
    
    // NORMALIZE TO UNIT ENERGY
    // Total energy = sum of |h[n]|²
    // Scale to make sum equal to 1
    scale = std::sqrt(scale);
    for (int n = 0; n < 2 * len + 1; n++) {
        h[n] /= scale;
    }
}

// ============================================================================
// ThreadSafeBitFIFO Implementation - Thread-safe FIFO queue for passing bit blocks between threads
// ============================================================================

void ThreadSafeBitFIFO::push(std::vector<uint8_t>&& bits) {
    std::lock_guard<std::mutex> lock(mutex_);
    queue_.push(std::move(bits));
    cv_.notify_one();
}

bool ThreadSafeBitFIFO::pop(std::vector<uint8_t>& bits, double timeout_s) {
    std::unique_lock<std::mutex> lock(mutex_);
    if (!cv_.wait_for(lock, std::chrono::duration<double>(timeout_s),
                      [this] { return !queue_.empty() || g_stop; })) {
        return false;
    }
    if (g_stop && queue_.empty()) {
        return false;
    }
    bits = std::move(queue_.front());
    queue_.pop();
    return true;
}

void ThreadSafeBitFIFO::notify_stop() {
    cv_.notify_all();
}

// ============================================================================
// RRCPulseShaper Implementation - Root Raised Cosine pulse shaping using improved truncated RRC
// Root Raised Cosine pulse shaping using improved truncated RRC
// ============================================================================

RRCPulseShaper::RRCPulseShaper(double symbol_rate, double sample_rate, 
                               double excess_bandwidth, int num_symbols)
    : symbol_rate_(symbol_rate), 
      sample_rate_(sample_rate),
      excess_bw_(excess_bandwidth) {
    
    // CALCULATE SAMPLES PER SYMBOL
    double sps = sample_rate / symbol_rate;
    samples_per_symbol_ = static_cast<int>(std::round(sps));
    
    // FIND INTEGER RATIO U/D THAT APPROXIMATES samples_per_symbol
    // For common case: 1MHz / 800kHz = 1.25 = 5/4
    // This allows exact representation without floating-point error
    find_ud_ratio(sps, U_, D_);
    
    // RECALCULATE ACTUAL SAMPLES PER SYMBOL FROM INTEGER RATIO
    // This may differ slightly from rounded value for better accuracy
    samples_per_symbol_ = U_;  // Use U as samples per symbol for upsampling
    
    // BETA (rolloff factor) is the excess bandwidth
    beta_ = excess_bandwidth;
    
    // VERIFY: Beta should equal (U-D)/D for consistency
    float beta_from_ud = static_cast<float>(U_ - D_) / D_;
    if (std::abs(beta_from_ud - beta_) > 0.01f) {
        std::cout << "Beta mismatch: " << beta_ << " != " << beta_from_ud << "\n";
    }
    
    // FILTER LENGTH: Span multiple symbol periods
    // len is half-length, total taps = 2*len+1
    int len = (U_ * num_symbols) / 2;
    num_taps_ = 2 * len + 1;
    
    // DESIGN RRC FILTER using new function
    design_rrc_taps(len);
    
    std::cout << "RRC Filter designed:\n"
              << "  Taps: " << num_taps_ << "\n"
              << "  U (samples per symbol): " << U_ << "\n"
              << "  D (downsampling factor): " << D_ << "\n"
              << "  U/D ratio: " << static_cast<float>(U_) / D_ << "\n"
              << "  Beta (rolloff factor): " << beta_ << "\n"
              << "  Filter half-length: " << len << "\n";
}

void RRCPulseShaper::find_ud_ratio(double ratio, int& U, int& D) {
    // FIND INTEGER RATIO U/D that approximates the given ratio
    // Constraint: D <= U <= 2*D (from rrc_pulse requirements)
    //
    // Common ratios:
    //   1.00 = 1/1, 2/2, 3/3, ...
    //   1.25 = 5/4, 10/8, 15/12, ...
    //   1.50 = 3/2, 6/4, 9/6, ...
    //   1.67 = 5/3, 10/6, ...
    //   2.00 = 2/1, 4/2, 6/3, ...
    
    // Try denominators from 1 to 20 and find best approximation
    double best_error = 1e10;
    int best_U = 1, best_D = 1;
    
    for (int d = 1; d <= 20; d++) {
        int u = static_cast<int>(std::round(ratio * d));
        
        // Check constraint: D <= U <= 2*D
        if (u < d || u > 2 * d) continue;
        
        double approx_ratio = static_cast<double>(u) / d;
        double error = std::abs(approx_ratio - ratio);
        
        if (error < best_error) {
            best_error = error;
            best_U = u;
            best_D = d;
        }
        
        // If exact match found, use it
        if (error < 1e-10) break;
    }
    
    U = best_U;
    D = best_D;
    
    double actual_ratio = static_cast<double>(U) / D;
    if (std::abs(actual_ratio - ratio) > 0.01) {
        std::cout << "Unable to find exact U/D ratio.\n"
    }

void RRCPulseShaper::design_rrc_taps(int len) {
    // ALLOCATE FILTER TAPS
    taps_.resize(2 * len + 1);
    
    // ALLOCATE TEMPORARY COMPLEX ARRAY for rrc_pulse function
    std::vector<std::complex<float>> h_complex(2 * len + 1);
    
    // GENERATE RRC IMPULSE RESPONSE using new function
    // Pass beta explicitly to avoid calculation issues
    rrc_pulse(h_complex.data(), len, U_, D_, beta_);
    
    // EXTRACT REAL PART (RRC for real-valued signals)
    // For BPSK, we use real-valued RRC on real symbols
    // The complex version would be used for complex modulations like QPSK
    for (int i = 0; i < 2 * len + 1; i++) {
        taps_[i] = h_complex[i].real();
    }
    
    // VERIFICATION: Check filter energy (should be ~1.0 after normalization)
    float energy = 0.0f;
    for (const auto& tap : taps_) {
        energy += tap * tap;
    }
    
    std::cout << "  Filter energy: " << energy 
              << " (should be ~1.0)\n";
    
    // OPTIONAL: Print first few taps for debugging
    if (taps_.size() <= 21) {
        std::cout << "  Filter taps: [";
        for (size_t i = 0; i < taps_.size(); i++) {
            std::cout << std::fixed << std::setprecision(4) << taps_[i];
            if (i < taps_.size() - 1) std::cout << ", ";
        }
        std::cout << "]\n";
    }
}

std::vector<std::complex<float>> RRCPulseShaper::shape(const std::vector<std::complex<float>>& symbols) {
    // PULSE SHAPING PROCESS
    // 1. Upsample symbols by inserting zeros
    // 2. Convolve with RRC filter
    // Bandwidth-limited continuous waveform is generated
    
    size_t num_symbols = symbols.size();
    size_t upsampled_len = num_symbols * samples_per_symbol_;
    
    // STEP 1: UPSAMPLING
    std::vector<std::complex<float>> upsampled(upsampled_len, {0.0f, 0.0f});
    for (size_t i = 0; i < num_symbols; i++) {
        upsampled[i * samples_per_symbol_] = symbols[i];
    }
    
    // STEP 2: Impulse response of RRC filter is convolved with the upsampled symbols
    size_t output_len = upsampled_len + num_taps_ - 1;
    std::vector<std::complex<float>> filtered(output_len, {0.0f, 0.0f});
    
    for (size_t n = 0; n < output_len; n++) {
        std::complex<float> sum = {0.0f, 0.0f};
        for (int k = 0; k < num_taps_; k++) {
            int idx = static_cast<int>(n) - k;
            if (idx >= 0 && idx < static_cast<int>(upsampled_len)) {
                sum += upsampled[idx] * taps_[k];
            }
        }
        filtered[n] = sum;
    }
    
    return filtered;
}

int RRCPulseShaper::get_samples_per_symbol() const { 
    return samples_per_symbol_; 
}

// ============================================================================
// GaussianFilter Implementation
// ============================================================================

GaussianFilter::GaussianFilter(double symbol_rate, double sample_rate,
                               double bt_product, int num_symbols)
    : symbol_rate_(symbol_rate), sample_rate_(sample_rate),
      bt_product_(bt_product) {
    
    samples_per_symbol_ = static_cast<int>(std::round(sample_rate / symbol_rate));
    num_taps_ = samples_per_symbol_ * num_symbols;
    if (num_taps_ % 2 == 0) num_taps_++;
    
    design_gaussian_taps();
    
    std::cout << "Gaussian Filter designed:\n"
              << "  Taps: " << num_taps_ << "\n"
              << "  Samples per symbol: " << samples_per_symbol_ << "\n"
              << "  BT product: " << bt_product_ << "\n";
}

void GaussianFilter::design_gaussian_taps() {
    taps_.resize(num_taps_);
    int center = num_taps_ / 2;
    double Ts = 1.0 / symbol_rate_;
    double B = bt_product_ / Ts;
    double alpha = std::sqrt(std::log(2.0)) / (M_PI * B * Ts);
    
    double sum = 0.0;
    for (int i = 0; i < num_taps_; i++) {
        double t = (i - center) / sample_rate_;
        double val = std::exp(-std::pow(t / alpha, 2));
        taps_[i] = static_cast<float>(val);
        sum += val;
    }
    
    for (auto& tap : taps_) {
        tap /= sum;
    }
}

std::vector<float> GaussianFilter::filter(const std::vector<float>& nrz_signal) {
    size_t output_len = nrz_signal.size() + num_taps_ - 1;
    std::vector<float> filtered(output_len, 0.0f);
    
    for (size_t n = 0; n < output_len; n++) {
        float sum = 0.0f;
        for (int k = 0; k < num_taps_; k++) {
            int idx = static_cast<int>(n) - k;
            if (idx >= 0 && idx < static_cast<int>(nrz_signal.size())) {
                sum += nrz_signal[idx] * taps_[k];
            }
        }
        filtered[n] = sum;
    }
    
    return filtered;
}

int GaussianFilter::get_samples_per_symbol() const {
    return samples_per_symbol_;
}

// ============================================================================
// BPSK Modulator Implementation
// ============================================================================

BPSKModulator::BPSKModulator(std::shared_ptr<RRCPulseShaper> shaper)
    : pulse_shaper_(shaper) {}

std::vector<std::complex<float>> BPSKModulator::bits_to_symbols(const std::vector<uint8_t>& bits) {
    std::vector<std::complex<float>> symbols(bits.size());
    for (size_t i = 0; i < bits.size(); i++) {
        symbols[i] = {(bits[i] == 0) ? 1.0f : -1.0f, 0.0f};
    }
    return symbols;
}

std::vector<std::complex<float>> BPSKModulator::modulate(const std::vector<uint8_t>& bits) {
    auto symbols = bits_to_symbols(bits);
    return pulse_shaper_->shape(symbols);
}

// ============================================================================
// GMSK Modulator Implementation
// ============================================================================

GMSKModulator::GMSKModulator(double symbol_rate, double sample_rate, double bt_product)
    : sample_rate_(sample_rate), symbol_rate_(symbol_rate), phase_state_(0.0f) {
    
    samples_per_symbol_ = static_cast<int>(std::round(sample_rate / symbol_rate));
    gaussian_filter_ = std::make_shared<GaussianFilter>(symbol_rate, sample_rate, bt_product, 4);
}

std::vector<std::complex<float>> GMSKModulator::modulate(const std::vector<uint8_t>& bits) {
    std::vector<float> nrz(bits.size() * samples_per_symbol_);
    for (size_t i = 0; i < bits.size(); i++) {
        float val = (bits[i] == 0) ? 1.0f : -1.0f;
        for (int j = 0; j < samples_per_symbol_; j++) {
            nrz[i * samples_per_symbol_ + j] = val;
        }
    }
    
    auto filtered = gaussian_filter_->filter(nrz);
    
    std::vector<std::complex<float>> samples(filtered.size());
    float h = 0.5f;
    
    for (size_t i = 0; i < filtered.size(); i++) {
        phase_state_ += M_PI * h * filtered[i] / samples_per_symbol_;
        samples[i] = {std::cos(phase_state_), std::sin(phase_state_)};
    }
    
    return samples;
}

// ============================================================================
// ===================== Modulator Factory Implementation =====================
// ============================================================================

std::unique_ptr<IModulator> ModulatorFactory::create(ModulationType type,
                                                     double symbol_rate,
                                                     double sample_rate,
                                                     double excess_bw,
                                                     double bt_product) {
    switch (type) {
        case ModulationType::BPSK: {
            auto shaper = std::make_shared<RRCPulseShaper>(symbol_rate, sample_rate, excess_bw);
            return std::make_unique<BPSKModulator>(shaper);
        }
        case ModulationType::GMSK: {
            return std::make_unique<GMSKModulator>(symbol_rate, sample_rate, bt_product);
        }
        default:
            throw std::runtime_error("Unknown modulation type");
    }
}

// ============================================================================
// SignalScaler Implementation with Histogram and Statistics
// ============================================================================
// PURPOSE: Prevent clipping and analyze signal distribution
//
// UHD CLIPPING PREVENTION:
// - UHD assumes sample values in [-1, 1]
// - Values outside this range get clipped → nonlinear distortions
// - BPSK symbols through RRC filter often exceed [-1, 1]
// - Solution: Scale down to target peak (default 0.7 = 70% of full scale)
//
// HISTOGRAM ANALYSIS:
// - Shows distribution of sample magnitudes
// - Helps determine optimal scaling factor
// - Reveals headroom margin
// - Detects potential clipping issues
// ============================================================================

float SignalScaler::analyze_peak(const std::vector<std::complex<float>>& samples) {
    // FIND MAXIMUM MAGNITUDE across all samples
    float peak = 0.0f;
    for (const auto& s : samples) {
        float mag = std::abs(s);
        if (mag > peak) peak = mag;
    }
    return peak;
}

SignalScaler::Statistics SignalScaler::compute_statistics(
    const std::vector<std::complex<float>>& samples) {
    // COMPUTE COMPREHENSIVE STATISTICS
    // Useful for understanding signal distribution and clipping risk
    
    Statistics stats;
    if (samples.empty()) return stats;
    
    // Initialize min/max
    stats.min_magnitude = std::numeric_limits<float>::max();
    stats.max_magnitude = 0.0f;
    
    // First pass: Compute min, max, mean
    float sum_magnitude = 0.0f;
    float sum_power = 0.0f;
    float max_power = 0.0f;
    
    for (const auto& s : samples) {
        float mag = std::abs(s);
        float power = std::norm(s);  // |s|^2
        
        if (mag < stats.min_magnitude) stats.min_magnitude = mag;
        if (mag > stats.max_magnitude) stats.max_magnitude = mag;
        
        sum_magnitude += mag;
        sum_power += power;
        if (power > max_power) max_power = power;
    }
    
    stats.mean_magnitude = sum_magnitude / samples.size();
    float mean_power = sum_power / samples.size();
    
    // Second pass: Compute standard deviation
    float sum_squared_diff = 0.0f;
    for (const auto& s : samples) {
        float mag = std::abs(s);
        float diff = mag - stats.mean_magnitude;
        sum_squared_diff += diff * diff;
    }
    stats.std_deviation = std::sqrt(sum_squared_diff / samples.size());
    
    // Compute PAPR
    if (mean_power > 1e-10f) {
        stats.papr_db = 10.0f * std::log10(max_power / mean_power);
    } else {
        stats.papr_db = 0.0f;
    }
    
    // Check for clipping (samples at or near ±1.0)
    const float clip_threshold = 0.99f;
    stats.samples_near_clipping = 0;
    for (const auto& s : samples) {
        float mag = std::abs(s);
        if (mag >= clip_threshold) {
            stats.samples_near_clipping++;
        }
    }
    
    return stats;
}

SignalScaler::Histogram SignalScaler::generate_histogram(
    const std::vector<std::complex<float>>& samples,
    int num_bins) {
    // GENERATE HISTOGRAM of sample magnitudes
    // Shows distribution: helps determine if scaling is appropriate
    
    Histogram hist;
    hist.num_bins = num_bins;
    hist.bins.resize(num_bins, 0);
    
    if (samples.empty()) return hist;
    
    // Find range
    float min_mag = std::numeric_limits<float>::max();
    float max_mag = 0.0f;
    for (const auto& s : samples) {
        float mag = std::abs(s);
        if (mag < min_mag) min_mag = mag;
        if (mag > max_mag) max_mag = mag;
    }
    
    hist.min_value = min_mag;
    hist.max_value = max_mag;
    
    // Handle edge case
    if (max_mag - min_mag < 1e-10f) {
        hist.bins[0] = samples.size();
        return hist;
    }
    
    // Compute bin edges
    hist.bin_edges.resize(num_bins + 1);
    float bin_width = (max_mag - min_mag) / num_bins;
    for (int i = 0; i <= num_bins; i++) {
        hist.bin_edges[i] = min_mag + i * bin_width;
    }
    
    // Fill histogram
    for (const auto& s : samples) {
        float mag = std::abs(s);
        int bin = static_cast<int>((mag - min_mag) / bin_width);
        if (bin >= num_bins) bin = num_bins - 1;  // Handle edge case
        if (bin < 0) bin = 0;
        hist.bins[bin]++;
    }
    
    return hist;
}

void SignalScaler::print_histogram(const Histogram& hist, 
                                   const std::string& title,
                                   int width) {
    // PRINT ASCII HISTOGRAM to console
    // Visual representation of sample distribution
    
    std::cout << "\n" << title << "\n";
    std::cout << std::string(width + 20, '=') << "\n";
    
    if (hist.bins.empty()) {
        std::cout << "No data\n";
        return;
    }
    
    // Find maximum count for scaling
    int max_count = *std::max_element(hist.bins.begin(), hist.bins.end());
    if (max_count == 0) {
        std::cout << "All bins empty\n";
        return;
    }
    
    // Print histogram bars
    for (int i = 0; i < hist.num_bins; i++) {
        float bin_start = hist.bin_edges[i];
        float bin_end = hist.bin_edges[i + 1];
        int count = hist.bins[i];
        
        // Scale bar length
        int bar_length = (count * width) / max_count;
        
        // Print range and bar
        std::cout << std::fixed << std::setprecision(3)
                  << "[" << std::setw(6) << bin_start 
                  << " - " << std::setw(6) << bin_end << "] "
                  << std::string(bar_length, '#')
                  << " " << count << "\n";
    }
    
    std::cout << std::string(width + 20, '=') << "\n";
    
    // Highlight clipping zone
    if (hist.max_value > 0.95f) {
        std::cout << "Samples approaching clipping threshold (1.0)\n";
    }
}

void SignalScaler::save_histogram(const Histogram& hist,
                                  const std::string& filename) {
    // SAVE HISTOGRAM DATA to file for external analysis
    // Format: bin_center, count
    
    std::ofstream ofs(filename);
    if (!ofs) {
        std::cerr << "Failed to open histogram file: " << filename << "\n";
        return;
    }
    
    ofs << "# Histogram of sample magnitudes\n";
    ofs << "# bin_center,count\n";
    
    for (int i = 0; i < hist.num_bins; i++) {
        float bin_center = (hist.bin_edges[i] + hist.bin_edges[i + 1]) / 2.0f;
        ofs << bin_center << "," << hist.bins[i] << "\n";
    }
    
    ofs.close();
    std::cout << "Histogram saved to: " << filename << "\n";
}

float SignalScaler::compute_scale_factor(const std::vector<std::complex<float>>& samples,
                                         float target_peak) {
    // COMPUTE SCALE FACTOR to achieve target peak
    // target_peak should be < 1.0 to prevent clipping
    // Typical values: 0.7 (30% headroom), 0.8 (20% headroom)
    
    float peak = analyze_peak(samples);
    if (peak < 1e-6f) return 1.0f;
    return target_peak / peak;
}

void SignalScaler::scale_samples(std::vector<std::complex<float>>& samples, 
                                 float scale_factor) {
    // APPLY SCALING to all samples
    for (auto& s : samples) {
        s *= scale_factor;
    }
}

float SignalScaler::compute_papr(const std::vector<std::complex<float>>& samples) {
    // COMPUTE PEAK-TO-AVERAGE POWER RATIO
    if (samples.empty()) return 0.0f;
    
    float peak_power = 0.0f;
    float avg_power = 0.0f;
    
    for (const auto& s : samples) {
        float power = std::norm(s);
        avg_power += power;
        if (power > peak_power) peak_power = power;
    }
    
    avg_power /= samples.size();
    
    if (avg_power < 1e-10f) return 0.0f;
    return 10.0f * std::log10(peak_power / avg_power);
}

int SignalScaler::count_clipped_samples(const std::vector<std::complex<float>>& samples,
                                        float clip_threshold) {
    // COUNT samples that would be clipped by UHD
    // UHD clips samples outside [-1, 1] range
    
    int count = 0;
    for (const auto& s : samples) {
        float mag = std::abs(s);
        if (mag >= clip_threshold) {
            count++;
        }
    }
    return count;
}

// ============================================================================
// ================== PSD (Power Spectral Density) Analysis ==================
// ============================================================================
// Verify that TX signal bandwidth is within 1 MHz specification
//
// METHOD:
// 1. Compute FFT of TX samples
// 2. Calculate power spectral density
// 3. Measure bandwidth (-3dB and 99% occupied)
// 4. Verify BW ≤ 1 MHz
//
// REQUIREMENT:
// "Check the PSD of the TX signal to verify that its bandwidth is within 1 MHz.
//  You may do so by capturing the TX samples sent to the USRP radio."
// ============================================================================

int SignalScaler::next_power_of_2(int n) {
    // Find next power of 2 >= n
    int p = 1;
    while (p < n) p *= 2;
    return p;
}

void SignalScaler::fft(std::vector<std::complex<float>>& x) {
    // SIMPLE FFT IMPLEMENTATION: Cooley-Tukey algorithm
    // This is a basic radix-2 FFT for power-of-2 sizes
    
    int N = x.size();
    if (N <= 1) return;
    
    // Check if power of 2
    if ((N & (N - 1)) != 0) {
        std::cerr << "FFT size must be power of 2\n";
        return;
    }
    
    // Bit-reversal permutation
    int j = 0;
    for (int i = 0; i < N - 1; i++) {
        if (i < j) std::swap(x[i], x[j]);
        
        int k = N / 2;
        while (k <= j) {
            j -= k;
            k /= 2;
        }
        j += k;
    }
    
    // Cooley-Tukey FFT
    for (int s = 1; s <= std::log2(N); s++) {
        int m = 1 << s;  // 2^s
        int m2 = m / 2;
        std::complex<float> w(1, 0);
        std::complex<float> wm = std::exp(std::complex<float>(0, -2.0f * M_PI / m));
        
        for (int j = 0; j < m2; j++) {
            for (int k = j; k < N; k += m) {
                std::complex<float> t = w * x[k + m2];
                std::complex<float> u = x[k];
                x[k] = u + t;
                x[k + m2] = u - t;
            }
            w *= wm;
        }
    }
}

SignalScaler::PSDResult SignalScaler::compute_psd(
    const std::vector<std::complex<float>>& samples,
    double sample_rate,
    int fft_size) {
    
    // COMPUTE POWER SPECTRAL DENSITY
    
    PSDResult result;
    
    if (samples.empty()) {
        std::cerr << "Empty samples for PSD computation\n";
        return result;
    }
    
    // Determine FFT size
    if (fft_size == 0) {
        // Auto: Use next power of 2 >= sample count (up to 4096)
        fft_size = std::min(next_power_of_2(samples.size()), 4096);
    } else {
        // Ensure power of 2
        fft_size = next_power_of_2(fft_size);
    }
    
    // Prepare data for FFT (zero-pad if necessary)
    std::vector<std::complex<float>> fft_data(fft_size, {0.0f, 0.0f});
    size_t copy_len = std::min(samples.size(), static_cast<size_t>(fft_size));
    std::copy(samples.begin(), samples.begin() + copy_len, fft_data.begin());
    
    // Apply window (Hamming) to reduce spectral leakage
    for (int i = 0; i < copy_len; i++) {
        float window = 0.54f - 0.46f * std::cos(2.0f * M_PI * i / (copy_len - 1));
        fft_data[i] *= window;
    }
    
    // Compute FFT
    fft(fft_data);
    
    // Compute PSD (power spectral density)
    result.psd_db.resize(fft_size);
    result.frequencies.resize(fft_size);
    
    float df = sample_rate / fft_size;  // Frequency resolution
    float normalization = 1.0f / (fft_size * fft_size);
    
    for (int i = 0; i < fft_size; i++) {
        // Power = |X[k]|^2
        float power = std::norm(fft_data[i]) * normalization;
        
        // Convert to dB
        if (power > 1e-20f) {
            result.psd_db[i] = 10.0f * std::log10(power);
        } else {
            result.psd_db[i] = -200.0f;  // Floor
        }
        
        // Frequency: FFT shift to center DC at 0
        if (i < fft_size / 2) {
            result.frequencies[i] = i * df;
        } else {
            result.frequencies[i] = (i - fft_size) * df;
        }
    }
    
    // FFT shift: Move zero frequency to center
    int mid = fft_size / 2;
    std::rotate(result.psd_db.begin(), result.psd_db.begin() + mid, result.psd_db.end());
    std::rotate(result.frequencies.begin(), result.frequencies.begin() + mid, result.frequencies.end());
    
    // Find peak power and frequency
    auto max_it = std::max_element(result.psd_db.begin(), result.psd_db.end());
    int peak_idx = std::distance(result.psd_db.begin(), max_it);
    result.peak_freq = result.frequencies[peak_idx];
    float peak_power_db = result.psd_db[peak_idx];
    
    // Calculate total power
    float total_power = 0.0f;
    for (const auto& psd_val : result.psd_db) {
        total_power += std::pow(10.0f, psd_val / 10.0f);
    }
    result.total_power_db = 10.0f * std::log10(total_power * df);
    
    // Measure -3dB bandwidth
    float threshold_3db = peak_power_db - 3.0f;
    int lower_idx = peak_idx;
    int upper_idx = peak_idx;
    
    // Find lower -3dB point
    while (lower_idx > 0 && result.psd_db[lower_idx] > threshold_3db) {
        lower_idx--;
    }
    
    // Find upper -3dB point
    while (upper_idx < fft_size - 1 && result.psd_db[upper_idx] > threshold_3db) {
        upper_idx++;
    }
    
    result.bw_3db = std::abs(result.frequencies[upper_idx] - result.frequencies[lower_idx]);
    
    // Measure occupied bandwidth (99% power)
    std::vector<std::pair<float, float>> power_freq_pairs;
    for (int i = 0; i < fft_size; i++) {
        float power = std::pow(10.0f, result.psd_db[i] / 10.0f);
        power_freq_pairs.push_back({power, result.frequencies[i]});
    }
    
    // Sort by power (descending)
    std::sort(power_freq_pairs.begin(), power_freq_pairs.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });
    
    // Find frequencies containing 99% of power
    float target_power = 0.99f * total_power;
    float accumulated_power = 0.0f;
    float min_freq = 1e9f, max_freq = -1e9f;
    
    for (const auto& pair : power_freq_pairs) {
        accumulated_power += pair.first;
        min_freq = std::min(min_freq, pair.second);
        max_freq = std::max(max_freq, pair.second);
        
        if (accumulated_power >= target_power) break;
    }
    
    result.bw_occupied_99 = max_freq - min_freq;
    
    // Check if bandwidth is within specification (≤ 1 MHz)
    result.bw_within_spec = (result.bw_occupied_99 <= 1.0e6f);
    
    return result;
}

void SignalScaler::print_psd_summary(const PSDResult& psd,
                                     double sample_rate,
                                     double symbol_rate) {
    // PRINT PSD ANALYSIS RESULTS
    
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "PSD (Power Spectral Density) Analysis\n";
    std::cout << std::string(70, '=') << "\n";
    
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  Sample rate:         " << sample_rate / 1e6 << " MHz\n";
    std::cout << "  Symbol rate:         " << symbol_rate / 1e3 << " kSps\n";
    std::cout << "  FFT size:            " << psd.frequencies.size() << "\n";
    std::cout << "  Frequency resolution:" << (sample_rate / psd.frequencies.size()) / 1e3 
              << " kHz\n";
    std::cout << "\n";
    
    std::cout << "  Peak frequency:      " << std::setprecision(3) 
              << psd.peak_freq / 1e3 << " kHz\n";
    std::cout << "  Total power:         " << std::setprecision(2) 
              << psd.total_power_db << " dB\n";
    std::cout << "\n";
    
    std::cout << "Bandwidth Measurements:\n";
    std::cout << "  -3dB bandwidth:      " << std::setprecision(1) 
              << psd.bw_3db / 1e3 << " kHz\n";
    std::cout << "  99% occupied BW:     " << psd.bw_occupied_99 / 1e3 << " kHz\n";
    std::cout << "\n";
    
    // Verification
    std::cout << "Bandwidth Specification Check:\n";
    std::cout << "  Required:            ≤ 1000 kHz (1 MHz)\n";
    std::cout << "  Measured:            " << psd.bw_occupied_99 / 1e3 << " kHz\n";
    
    if (psd.bw_within_spec) {
        std::cout << "  Status:              PASS - Bandwidth within specification\n";
        float margin = 1.0e6f - psd.bw_occupied_99;
        std::cout << "  Margin:              " << margin / 1e3 << " kHz\n";
    } else {
        std::cout << "  Status:              FAIL - Bandwidth exceeds specification!\n";
        float excess = psd.bw_occupied_99 - 1.0e6f;
        std::cout << "  Excess:              " << excess / 1e3 << " kHz\n";
        std::cout << " WARNING: TX bandwidth exceeds 1 MHz specification!\n";
        std::cout << "    Consider reducing symbol rate or excess bandwidth.\n";
    }
    
    std::cout << std::string(70, '=') << "\n";
}

void SignalScaler::save_psd(const PSDResult& psd,
                            const std::string& filename) {
    // SAVE PSD DATA to CSV for external plotting
    
    std::ofstream ofs(filename);
    if (!ofs) {
        std::cerr << "Failed to open PSD file: " << filename << "\n";
        return;
    }
    
    ofs << "# Power Spectral Density\n";
    ofs << "# frequency_Hz,psd_dB\n";
    ofs << std::fixed << std::setprecision(6);
    
    for (size_t i = 0; i < psd.frequencies.size(); i++) {
        ofs << psd.frequencies[i] << "," << psd.psd_db[i] << "\n";
    }
    
    ofs.close();
    std::cout << "PSD data saved to: " << filename << "\n";
}

USRPTransmitter::USRPTransmitter(const std::string& device_args,
                                 double center_freq,
                                 double sample_rate,
                                 double gain,
                                 const std::string& antenna,
                                 const std::string& subdev)
    : center_freq_(center_freq), sample_rate_(sample_rate), 
      gain_(gain), streaming_(false) {
    
    std::cout << "Creating USRP transmitter...\n";
    usrp_ = uhd::usrp::multi_usrp::make(device_args);
    
    if (!subdev.empty()) {
        usrp_->set_tx_subdev_spec(subdev);
    }
    
    std::cout << "Setting TX rate: " << sample_rate_ << " Hz\n";
    usrp_->set_tx_rate(sample_rate_);
    
    std::cout << "Setting TX freq: " << center_freq_ << " Hz\n";
    usrp_->set_tx_freq(center_freq_);
    
    std::cout << "Setting TX gain: " << gain_ << " dB\n";
    usrp_->set_tx_gain(gain_);
    
    std::cout << "Setting TX antenna: " << antenna << "\n";
    usrp_->set_tx_antenna(antenna);
    
    uhd::stream_args_t stream_args("fc32", "sc16");
    tx_streamer_ = usrp_->get_tx_stream(stream_args);
    
    metadata_.start_of_burst = false;
    metadata_.end_of_burst = false;
    metadata_.has_time_spec = false;
    
    std::cout << "USRP transmitter initialized successfully\n";
    std::cout << "Actual TX rate: " << usrp_->get_tx_rate() << " Hz\n";
    std::cout << "Actual TX freq: " << usrp_->get_tx_freq() << " Hz\n";
    std::cout << "Actual TX gain: " << usrp_->get_tx_gain() << " dB\n";
}

USRPTransmitter::~USRPTransmitter() {
    stop();
}

void USRPTransmitter::transmit(const std::vector<std::complex<float>>& samples) {
    if (samples.empty()) return;
    
    size_t num_sent = 0;
    const size_t total_samples = samples.size();
    
    while (num_sent < total_samples && !g_stop) {
        size_t num_to_send = total_samples - num_sent;
        size_t sent = tx_streamer_->send(
            samples.data() + num_sent,
            num_to_send,
            metadata_,
            3.0
        );
        
        if (sent == 0) {
            std::cerr << "Warning: TX timeout\n";
            break;
        }
        
        num_sent += sent;
    }
}

void USRPTransmitter::stop() {
    if (streaming_) {
        metadata_.end_of_burst = true;
        std::vector<std::complex<float>> dummy(1, {0.0f, 0.0f});
        tx_streamer_->send(dummy.data(), 1, metadata_, 1.0);
        streaming_ = false;
    }
}

double USRPTransmitter::get_sample_rate() const { 
    return usrp_->get_tx_rate(); 
}

// ============================================================================
// ===================== Thread Functions =====================
// ============================================================================


void bit_generator_thread(ThreadSafeBitFIFO& fifo, size_t bits_per_block) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, 1);
    
    std::cout << "Bit generator thread started\n";
    
    while (!g_stop) {
        std::vector<uint8_t> bits(bits_per_block);
        for (size_t i = 0; i < bits_per_block; i++) {
            bits[i] = static_cast<uint8_t>(dist(gen));
        }
        
        std::cout << "Generated " << bits_per_block << " random bits\n";
        fifo.push(std::move(bits));
        
        auto start = std::chrono::steady_clock::now();
        while (!g_stop) {
            auto elapsed = std::chrono::steady_clock::now() - start;
            if (elapsed >= std::chrono::seconds(1)) break;
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
    }
    
    fifo.notify_stop();
    std::cout << "Bit generator thread stopped\n";
}

void modulator_transmitter_thread(ThreadSafeBitFIFO& fifo,
                                  USRPTransmitter& transmitter,
                                  ModulationType mod_type,
                                  double symbol_rate,
                                  double sample_rate,
                                  const std::string& save_file,
                                  float target_peak,
                                  bool show_histogram,
                                  bool save_histogram_file,
                                  bool show_psd,
                                  bool save_psd_file) {
    std::cout << "Modulator/transmitter thread started\n";
    
    auto modulator = ModulatorFactory::create(mod_type, symbol_rate, sample_rate, 0.25, 0.3);
    std::cout << "Using " << modulator->get_name() << " modulation\n";
    std::cout << "Bits per symbol: " << modulator->get_bits_per_symbol() << "\n";
    std::cout << "Target peak for clipping prevention: " << target_peak << "\n";
    std::cout << "PSD analysis: " << (show_psd ? "Enabled" : "Disabled") << "\n\n";
    
    std::ofstream ofs;
    bool saving = false;
    if (!save_file.empty()) {
        ofs.open(save_file, std::ios::binary);
        if (ofs) {
            saving = true;
            std::cout << "Saving TX samples to: " << save_file << "\n";
        }
    }
    
    std::vector<uint8_t> bits;
    int block_count = 0;
    
    while (!g_stop) {
        if (!fifo.pop(bits, 0.5)) {
            continue;
        }
        
        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "Processing block " << ++block_count 
                  << " (" << bits.size() << " bits)\n";
        std::cout << std::string(70, '=') << "\n";
        
        // MODULATE
        auto samples = modulator->modulate(bits);
        std::cout << "  Generated " << samples.size() << " samples\n";
        
        // ====================================================================
        // PSD ANALYSIS (as per requirement) - BEFORE SCALING
        // ====================================================================
        if (show_psd && block_count == 1) {
            std::cout << "\n>>> Analyzing PSD of TX signal (BEFORE scaling)...\n";
            
            // Compute PSD
            auto psd_before = SignalScaler::compute_psd(samples, sample_rate, 2048);
            
            // Display summary
            SignalScaler::print_psd_summary(psd_before, sample_rate, symbol_rate);
            
            // Save if requested
            if (save_psd_file) {
                SignalScaler::save_psd(psd_before, "psd_before_scaling.csv");
            }
        }
        
        // ====================================================================
        // CLIPPING PREVENTION ANALYSIS
        // ====================================================================
        
        // BEFORE SCALING: Analyze original signal
        std::cout << "\n--- BEFORE SCALING ---\n";
        auto stats_before = SignalScaler::compute_statistics(samples);
        
        std::cout << "  Magnitude range: [" << std::fixed << std::setprecision(4)
                  << stats_before.min_magnitude << ", " 
                  << stats_before.max_magnitude << "]\n";
        std::cout << "  Mean magnitude:  " << stats_before.mean_magnitude << "\n";
        std::cout << "  Std deviation:   " << stats_before.std_deviation << "\n";
        std::cout << "  PAPR:            " << stats_before.papr_db << " dB\n";
        
        // Check if clipping would occur
        int would_clip = SignalScaler::count_clipped_samples(samples, 1.0f);
        if (would_clip > 0) {
            float clip_percent = 100.0f * would_clip / samples.size();
            std::cout << "WARNING: " << would_clip << " samples (" 
                      << std::setprecision(2) << clip_percent 
                      << "%) exceed UHD range [-1,1] and would be CLIPPED!\n";
        } else {
            std::cout << "No samples exceed UHD range [-1,1]\n";
        }
        
        // GENERATE HISTOGRAM (if enabled and first block)
        if (show_histogram && block_count == 1) {
            auto hist_before = SignalScaler::generate_histogram(samples, 20);
            SignalScaler::print_histogram(hist_before, 
                "Sample Magnitude Distribution (BEFORE SCALING)", 50);
            
            if (save_histogram_file) {
                SignalScaler::save_histogram(hist_before, "histogram_before_scaling.csv");
            }
        }
        
        // COMPUTE SCALING FACTOR
        float scale_factor = SignalScaler::compute_scale_factor(samples, target_peak);
        std::cout << "\n--- SCALING OPERATION ---\n";
        std::cout << "  Peak before:     " << stats_before.max_magnitude << "\n";
        std::cout << "  Target peak:     " << target_peak << "\n";
        std::cout << "  Scale factor:    " << std::setprecision(4) << scale_factor << "\n";
        std::cout << "  Power reduction: " << std::setprecision(2) 
                  << 20.0f * std::log10(scale_factor) << " dB\n";
        
        // APPLY SCALING
        SignalScaler::scale_samples(samples, scale_factor);
        
        // AFTER SCALING: Verify clipping prevention
        std::cout << "\n--- AFTER SCALING ---\n";
        auto stats_after = SignalScaler::compute_statistics(samples);
        
        std::cout << "  Magnitude range: [" << std::fixed << std::setprecision(4)
                  << stats_after.min_magnitude << ", " 
                  << stats_after.max_magnitude << "]\n";
        std::cout << "  Mean magnitude:  " << stats_after.mean_magnitude << "\n";
        std::cout << "  Std deviation:   " << stats_after.std_deviation << "\n";
        std::cout << "  PAPR:            " << stats_after.papr_db << " dB\n";
        
        // Verify no clipping
        int will_clip = SignalScaler::count_clipped_samples(samples, 1.0f);
        if (will_clip > 0) {
            std::cout << " ERROR: " << will_clip 
                      << " samples STILL exceed [-1,1]! Increase scaling!\n";
        } else {
            float headroom_db = 20.0f * std::log10(1.0f / stats_after.max_magnitude);
            std::cout << "  ✓ All samples within UHD range [-1,1]\n";
            std::cout << "  Headroom:        " << std::setprecision(2) 
                      << headroom_db << " dB\n";
        }
        
        // Check if samples are near clipping threshold (>0.95)
        if (stats_after.samples_near_clipping > 0) {
            float near_clip_percent = 100.0f * stats_after.samples_near_clipping / samples.size();
            std::cout << stats_after.samples_near_clipping << " samples ("
                      << std::setprecision(1) << near_clip_percent 
                      << "%) near clipping (>0.95)\n";
        }
        
        // HISTOGRAM AFTER SCALING
        if (show_histogram && block_count == 1) {
            auto hist_after = SignalScaler::generate_histogram(samples, 20);
            SignalScaler::print_histogram(hist_after,
                "Sample Magnitude Distribution (AFTER SCALING)", 50);
            
            if (save_histogram_file) {
                SignalScaler::save_histogram(hist_after, "histogram_after_scaling.csv");
            }
        }
        
        // ====================================================================
        // PSD ANALYSIS (as per requirement) - AFTER SCALING
        // ====================================================================
        if (show_psd && block_count == 1) {
            std::cout << "\n>>> Analyzing PSD of TX signal (AFTER scaling)...\n";
            
            // Compute PSD
            auto psd_after = SignalScaler::compute_psd(samples, sample_rate, 2048);
            
            // Display summary with bandwidth verification
            SignalScaler::print_psd_summary(psd_after, sample_rate, symbol_rate);
            
            // Save if requested
            if (save_psd_file) {
                SignalScaler::save_psd(psd_after, "psd_after_scaling.csv");
            }
            
            // Additional verification message
            if (psd_after.bw_within_spec) {
                std::cout << "\nBANDWIDTH VERIFICATION PASSED\n";
                std::cout << "  TX signal bandwidth is within 1 MHz specification.\n";
            } else {
                std::cout << "\nBANDWIDTH VERIFICATION FAILED\n";
                std::cout << "  TX signal bandwidth EXCEEDS 1 MHz specification!\n";
                std::cout << "  Action required: Reduce symbol rate or excess bandwidth.\n";
            }
        }
        
        std::cout << "\n--- TRANSMISSION ---\n";
        std::cout << "  Transmitting " << samples.size() << " samples to USRP...\n";
        
        // TRANSMIT
        transmitter.transmit(samples);
        std::cout << " Transmission complete\n";

        
        // SAVE TO FILE (optional)
        if (saving) {
            ofs.write(reinterpret_cast<const char*>(samples.data()),
                     samples.size() * sizeof(std::complex<float>));
        }
        
        std::cout << std::string(70, '=') << "\n";
    }
    
    if (saving) {
        ofs.close();
        std::cout << "Saved " << block_count << " blocks to " << save_file << "\n";
    }
    
    std::cout << "Modulator/transmitter thread stopped\n";
}

ModulationType string_to_modulation_type(const std::string& str) {
    std::string lower = str;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    
    if (lower == "bpsk") return ModulationType::BPSK;
    if (lower == "gmsk") return ModulationType::GMSK;
    
    throw std::runtime_error("Unknown modulation type: " + str + " (supported: bpsk, gmsk)");
}

std::string modulation_type_to_string(ModulationType type) {
    switch (type) {
        case ModulationType::BPSK: return "BPSK";
        case ModulationType::GMSK: return "GMSK";
        default: return "Unknown";
    }
}

// ============================================================================
// ======================== Main Function ========================
// ============================================================================

int main(int argc, char* argv[]) {
    std::string device_args = "type=b200";
    double center_freq = 2.412e9;
    double sample_rate = 1.0e6;
    double symbol_rate = 800e3;
    double tx_gain = 30.0;
    size_t bits_per_block = 1000;  // As per README.md requirement
    std::string save_file = "";
    ModulationType mod_type = ModulationType::BPSK;
    
    // NEW: Clipping prevention parameters
    float target_peak = 0.7f;           // Default: 70% of full scale (30% headroom)
    bool show_histogram = true;         // Show histogram for first block
    bool save_histogram_file = false;   // Save histogram to CSV
    
    // NEW: PSD analysis parameters
    bool show_psd = true;               // Show PSD analysis for first block
    bool save_psd_file = false;         // Save PSD data to CSV
    
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--dev" && i + 1 < argc) {
            device_args = argv[++i];
        } else if (arg == "--fc" && i + 1 < argc) {
            center_freq = std::stod(argv[++i]);
        } else if (arg == "--fs" && i + 1 < argc) {
            sample_rate = std::stod(argv[++i]);
        } else if (arg == "--fsym" && i + 1 < argc) {
            symbol_rate = std::stod(argv[++i]);
        } else if (arg == "--gain" && i + 1 < argc) {
            tx_gain = std::stod(argv[++i]);
        } else if (arg == "--nbits" && i + 1 < argc) {
            bits_per_block = std::stoull(argv[++i]);
        } else if (arg == "--save" && i + 1 < argc) {
            save_file = argv[++i];
        } else if (arg == "--mod" && i + 1 < argc) {
            try {
                mod_type = string_to_modulation_type(argv[++i]);
            } catch (const std::exception& e) {
                std::cerr << "Error: " << e.what() << "\n";
                return 1;
            }
        } else if (arg == "--peak" && i + 1 < argc) {
            // NEW: Set target peak for clipping prevention
            target_peak = std::stof(argv[++i]);
            if (target_peak <= 0.0f || target_peak >= 1.0f) {
                std::cerr << "Error: --peak must be between 0.0 and 1.0\n";
                return 1;
            }
        } else if (arg == "--no-histogram") {
            // NEW: Disable histogram display
            show_histogram = false;
        } else if (arg == "--save-histogram") {
            // NEW: Save histogram to CSV files
            save_histogram_file = true;
        } else if (arg == "--no-psd") {
            // NEW: Disable PSD analysis
            show_psd = false;
        } else if (arg == "--save-psd") {
            // NEW: Save PSD to CSV files
            save_psd_file = true;
        } else if (arg == "--help") {
            std::cout << "Usage: " << argv[0] << " [options]\n"
                      << "\nBPSK/GMSK Transmitter with Enhanced Clipping Prevention\n"
                      << "\nBasic Options:\n"
                      << "  --dev <args>     Device arguments (default: type=b200)\n"
                      << "                   For N210: type=usrp2 or addr=192.168.10.2\n"
                      << "  --fc <Hz>        Center frequency in Hz (default: 2.412e9)\n"
                      << "  --fs <Hz>        Sample rate in Hz (default: 1e6)\n"
                      << "  --fsym <Hz>      Symbol rate in Hz (default: 800e3)\n"
                      << "  --gain <dB>      TX gain in dB (default: 30.0)\n"
                      << "  --nbits <n>      Bits per block (default: 1000)\n"
                      << "  --mod <type>     Modulation: bpsk or gmsk (default: bpsk)\n"
                      << "  --save <file>    Save TX samples to binary file (optional)\n"
                      << "\nClipping Prevention Options:\n"
                      << "  --peak <value>   Target peak value (0.0-1.0, default: 0.7)\n"
                      << "                   Lower value = more headroom, less TX power\n"
                      << "                   Higher value = less headroom, more TX power\n"
                      << "                   Recommended: 0.6-0.8\n"
                      << "  --no-histogram   Disable histogram display\n"
                      << "  --save-histogram Save histogram data to CSV files\n"
                      << "\nPSD Analysis Options:\n"
                      << "  --no-psd         Disable PSD analysis\n"
                      << "  --save-psd       Save PSD data to CSV files\n"
                      << "\nImportant:\n"
                      << "  UHD assumes samples in [-1, 1] range. Values outside are clipped.\n"
                      << "  BPSK symbols through RRC filter often exceed this range.\n"
                      << "  This program automatically scales samples to prevent clipping.\n"
                      << "  Use --peak to adjust the tradeoff between headroom and TX power.\n"
                      << "\nExamples:\n"
                      << "  # Basic usage with default clipping prevention (peak=0.7)\n"
                      << "  " << argv[0] << " --mod bpsk --fc 2.412e9 --gain 30\n\n"
                      << "  # More aggressive scaling for better clipping protection\n"
                      << "  " << argv[0] << " --peak 0.6\n\n"
                      << "  # Less scaling for higher TX power (risky!)\n"
                      << "  " << argv[0] << " --peak 0.8\n\n"
                      << "  # Save histogram for analysis\n"
                      << "  " << argv[0] << " --save-histogram --save samples.dat\n\n"
                      << "  # N210 with custom peak\n"
                      << "  " << argv[0] << " --dev type=usrp2 --peak 0.65\n";
            return 0;
        }
    }
    
    try {
        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "  BPSK/GMSK Transmitter with Clipping Prevention\n";
        std::cout << std::string(70, '=') << "\n";
        std::cout << "Modulation:       " << modulation_type_to_string(mod_type) << "\n";
        std::cout << "Center frequency: " << center_freq / 1e9 << " GHz\n";
        std::cout << "Sample rate:      " << sample_rate / 1e6 << " MHz\n";
        std::cout << "Symbol rate:      " << symbol_rate / 1e3 << " kSps\n";
        std::cout << "TX gain:          " << tx_gain << " dB\n";
        std::cout << "Bits per block:   " << bits_per_block << "\n";
        std::cout << "\nClipping Prevention:\n";
        std::cout << "  Target peak:    " << target_peak << " (" 
                  << std::fixed << std::setprecision(1)
                  << (100.0f * target_peak) << "% of full scale)\n";
        std::cout << "  Headroom:       " << std::setprecision(2)
                  << 20.0f * std::log10(1.0f / target_peak) << " dB\n";
        std::cout << "  Histogram:      " << (show_histogram ? "Enabled" : "Disabled") << "\n";
        std::cout << "  Save histogram: " << (save_histogram_file ? "Yes" : "No") << "\n";
        std::cout << "\nPSD Analysis:\n";
        std::cout << "  PSD display:    " << (show_psd ? "Enabled" : "Disabled") << "\n";
        std::cout << "  Save PSD:       " << (save_psd_file ? "Yes" : "No") << "\n";
        if (show_psd) {
            std::cout << "  Note: PSD analysis verifies BW ≤ 1 MHz specification\n";
        }
        std::cout << std::string(70, '=') << "\n\n";
        
        std::signal(SIGINT, signal_handler);
        std::signal(SIGTERM, signal_handler);
        
        USRPTransmitter transmitter(device_args, center_freq, sample_rate, 
                                   tx_gain, "TX/RX", "A:0");
        
        ThreadSafeBitFIFO bit_fifo;
        
        std::thread bit_gen(bit_generator_thread, std::ref(bit_fifo), bits_per_block);
        std::thread mod_tx(modulator_transmitter_thread, std::ref(bit_fifo), 
                          std::ref(transmitter), mod_type, symbol_rate, sample_rate, 
                          save_file, target_peak, show_histogram, save_histogram_file,
                          show_psd, save_psd_file);
        
        std::cout << "Transmitter running. Press Ctrl+C to stop.\n\n";
        
        bit_gen.join();
        mod_tx.join();
        
        std::cout << "\nTransmitter stopped cleanly.\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}