// ============================================================================
// usrp_interface.hpp - USRP Radio Interface with Internal Buffer
// ============================================================================
#ifndef RECEIVER_USRP_HPP
#define RECEIVER_USRP_HPP

#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/stream.hpp>
#include <complex>
#include <vector>
#include <string>
#include <atomic>

class ReceiverUSRP {
public:
    ReceiverUSRP(const std::string& device_args,
                  double center_freq_hz,
                  double sample_rate_sps,
                  double gain_db,
                  const std::string& antenna = "TX/RX",
                  const std::string& subdev = "");
    ~ReceiverUSRP()

    // configure or reconfigure runtime parameters (optional)
    void set_center_freq(double freq_hz);
    void set_sample_rate(double rate_sps);
    void set_gain(double gain_db);

    // start/stop streaming
    void start();
    void stop();

    // blocking read for up to max_samples complex<float> IQ samples.
    // returns number of samples written to 'out' (0 on timeout).
    // thread-safe if called from one reader thread.
    size_t read(std::complex<float>* out, size_t max_samples, double timeout_s = 0.2);

    bool read_exact(std::complex<float>* out, size_t n, double total_timeout_s = 2.0);

    double actual_sample_rate() const { return sps_; }
    double actual_center_freq() const { return freq_; }
    double actual_gain() const { return gain_; }

private:
    uhd::usrp::multi_usrp::sptr usrp_;
    uhd::rx_streamer::sptr rx_streamer_;
    uhd::rx_metadata_t md_{};

    double freq_{0.0};
    double sps_{0.0};
    double gain_{0.0};
    std::string ant_{"TX/RX"};
    std::string subdev_{};

    std::atomic<bool> streaming_{false};
};

#endif // USRP_INTERFACE_HPP
