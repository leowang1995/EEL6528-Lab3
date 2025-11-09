
// ============================================================================
// transmitter_alice.hpp - BPSK/GMSK Transmitter Header File
// ============================================================================
#ifndef TRANSMITTER_ALICE_HPP
#define TRANSMITTER_ALICE_HPP

#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/utils/thread.hpp>
#include <vector>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <complex>
#include <atomic>
#include <string>
#include <memory>

// ============================================================================
// Modulation Types
// ============================================================================
enum class ModulationType {
    BPSK,
    GMSK
};

// ============================================================================
// Global stop flag
// ============================================================================
extern std::atomic<bool> g_stop;

// Signal handler
void signal_handler(int);

// ============================================================================
// Thread-safe FIFO for bit blocks
// ============================================================================
class ThreadSafeBitFIFO {
private:
    std::queue<std::vector<uint8_t>> queue_;
    std::mutex mutex_;
    std::condition_variable cv_;

public:
    void push(std::vector<uint8_t>&& bits);
    bool pop(std::vector<uint8_t>& bits, double timeout_s = 0.5);
    void notify_stop();
};

// ============================================================================
// RRC Pulse Shaper - Generates RRC filter and performs pulse shaping
// ============================================================================
class RRCPulseShaper {
private:
    std::vector<float> taps_;
    double symbol_rate_;
    double sample_rate_;
    int samples_per_symbol_;
    double excess_bw_;
    int num_taps_;

    void design_rrc_taps();

public:
    RRCPulseShaper(double symbol_rate, double sample_rate, 
                   double excess_bandwidth = 0.25, int num_symbols = 8);
    
    // Upsample and filter: converts symbols to samples
    std::vector<std::complex<float>> shape(const std::vector<std::complex<float>>& symbols);
    
    int get_samples_per_symbol() const;
    const std::vector<float>& get_taps() const { return taps_; }
};

// ============================================================================
// Gaussian Filter for GMSK
// ============================================================================
class GaussianFilter {
private:
    std::vector<float> taps_;
    double symbol_rate_;
    double sample_rate_;
    int samples_per_symbol_;
    double bt_product_;  // BT product (bandwidth-time)
    int num_taps_;

    void design_gaussian_taps();

public:
    GaussianFilter(double symbol_rate, double sample_rate,
                   double bt_product = 0.3, int num_symbols = 4);
    
    // Filter NRZ bits to create smooth frequency modulation
    std::vector<float> filter(const std::vector<float>& nrz_signal);
    
    int get_samples_per_symbol() const;
};

// ============================================================================
// Base Modulator Interface
// ============================================================================
class IModulator {
public:
    virtual ~IModulator() = default;
    virtual std::vector<std::complex<float>> modulate(const std::vector<uint8_t>& bits) = 0;
    virtual int get_bits_per_symbol() const = 0;
    virtual std::string get_name() const = 0;
};

// ============================================================================
// BPSK Modulator - Binary Phase Shift Keying
// ============================================================================
class BPSKModulator : public IModulator {
private:
    std::shared_ptr<RRCPulseShaper> pulse_shaper_;

public:
    explicit BPSKModulator(std::shared_ptr<RRCPulseShaper> shaper);
    
    std::vector<std::complex<float>> modulate(const std::vector<uint8_t>& bits) override;
    int get_bits_per_symbol() const override { return 1; }
    std::string get_name() const override { return "BPSK"; }
    
    static std::vector<std::complex<float>> bits_to_symbols(const std::vector<uint8_t>& bits);
};

// ============================================================================
// GMSK Modulator - Gaussian Minimum Shift Keying
// ============================================================================
class GMSKModulator : public IModulator {
private:
    std::shared_ptr<GaussianFilter> gaussian_filter_;
    double sample_rate_;
    double symbol_rate_;
    int samples_per_symbol_;
    float phase_state_;  // Accumulated phase for continuous phase modulation

public:
    GMSKModulator(double symbol_rate, double sample_rate, double bt_product = 0.3);
    
    std::vector<std::complex<float>> modulate(const std::vector<uint8_t>& bits) override;
    int get_bits_per_symbol() const override { return 1; }
    std::string get_name() const override { return "GMSK"; }
    
    void reset_phase() { phase_state_ = 0.0f; }
};

// ============================================================================
// Modulator Factory
// ============================================================================
class ModulatorFactory {
public:
    static std::unique_ptr<IModulator> create(ModulationType type,
                                              double symbol_rate,
                                              double sample_rate,
                                              double excess_bw = 0.25,
                                              double bt_product = 0.3);
};

// ============================================================================
// Signal Scaler - Prevents clipping by analyzing and scaling samples
// ============================================================================
class SignalScaler {
public:
    static float analyze_peak(const std::vector<std::complex<float>>& samples);
    
    static float compute_scale_factor(const std::vector<std::complex<float>>& samples,
                                      float target_peak = 0.8f);
    
    static void scale_samples(std::vector<std::complex<float>>& samples, 
                             float scale_factor);
    
    // Compute PAPR (Peak to Average Power Ratio)
    static float compute_papr(const std::vector<std::complex<float>>& samples);
};

// ============================================================================
// USRP Transmitter - Interfaces with USRP hardware for transmission
// ============================================================================
class USRPTransmitter {
private:
    uhd::usrp::multi_usrp::sptr usrp_;
    uhd::tx_streamer::sptr tx_streamer_;
    uhd::tx_metadata_t metadata_;
    double center_freq_;
    double sample_rate_;
    double gain_;
    bool streaming_;

public:
    USRPTransmitter(const std::string& device_args,
                   double center_freq,
                   double sample_rate,
                   double gain,
                   const std::string& antenna = "TX/RX",
                   const std::string& subdev = "A:0");
    
    ~USRPTransmitter();

    void transmit(const std::vector<std::complex<float>>& samples);
    void stop();
    double get_sample_rate() const;
};

// ============================================================================
// Thread functions
// ============================================================================
void bit_generator_thread(ThreadSafeBitFIFO& fifo, size_t bits_per_block);

void modulator_transmitter_thread(ThreadSafeBitFIFO& fifo,
                                  USRPTransmitter& transmitter,
                                  ModulationType mod_type,
                                  double symbol_rate,
                                  double sample_rate,
                                  const std::string& save_file = "");

// ============================================================================
// Utility functions
// ============================================================================
ModulationType string_to_modulation_type(const std::string& str);
std::string modulation_type_to_string(ModulationType type);

#endif // TRANSMITTER_ALICE_HPP