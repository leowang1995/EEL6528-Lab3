# Lab 3 Receiver Design - Energy Detector

## Design Rationale

### 1. RX Sampling Rate Selection: **1 MHz**

**Choice:** RX sampling rate = 1 MHz (same as TX)

**Rationale:**
- **TX signal characteristics:**
  - Symbol rate: 800 kSps
  - RRC excess bandwidth (β): 0.25
  - Occupied bandwidth: 800 kHz × (1 + 0.25) = **1 MHz**

- **Nyquist criterion:** 
  - Minimum sampling rate = Signal bandwidth = 1 MHz
  - Using 1 MHz satisfies Nyquist with no excess margin

- **Advantages of 1 MHz sampling:**
  1. **No resampling needed** - eliminates computational overhead and complexity
  2. **Exact match to TX rate** - simplifies testing and analysis
  3. **No aliasing** - signal bandwidth equals sample rate
  4. **Minimal processing** - fewer samples to process vs. oversampling
  5. **Power efficiency** - lower sample rate = lower processing power

- **Why not higher rates?**
  - 2 MHz or higher would require downsampling/decimation before processing
  - Increases computational cost without benefit for energy detection
  - More data to store if saving samples

**Conclusion:** 1 MHz is the optimal choice - meets Nyquist, matches TX, no resampling overhead.

---

### 2. Filtering Before Energy Detector: **Optional (Disabled by Default)**

**Choice:** No filtering by default, optional low-pass filter available with `--filter` flag

**Rationale:**

#### Energy Detection vs. Demodulation
- **Our goal:** Detect packet presence (energy detection)
- **Not needed:** Symbol recovery, matched filtering, clock recovery
- **Implication:** Matched RRC filtering is NOT necessary

#### Three Filtering Options Considered:

##### A) **No Filtering** (Default Choice) ✅
**Pros:**
- ✅ Simplest implementation
- ✅ Lowest latency
- ✅ No filter distortion
- ✅ Sufficient for clean channel (lab environment)
- ✅ Direct energy measurement of received signal

**Cons:**
- ❌ Susceptible to out-of-band interference
- ❌ Includes noise from full Nyquist bandwidth

**When to use:** Clean lab environment, no strong interferers

##### B) **Low-Pass Filter** (Available with `--filter`)
**Pros:**
- ✅ Rejects out-of-band noise
- ✅ Improves SNR in noisy environment
- ✅ Simple FIR filter (101 taps)
- ✅ Cutoff at 500 kHz (signal BW/2)

**Cons:**
- ❌ Adds processing complexity
- ❌ Introduces group delay
- ❌ May slightly distort signal edges

**When to use:** Noisy environment, other WiFi/BT devices nearby

##### C) **Matched RRC Filter** (Not Implemented)
**Pros:**
- ✅ Optimal SNR for symbol detection
- ✅ Required for demodulation

**Cons:**
- ❌ Unnecessary for energy detection
- ❌ More complex implementation
- ❌ Higher computational cost
- ❌ Requires precise timing alignment

**Why not used:** Overkill for simple packet detection

#### **Design Decision:**

**Default: No filtering** because:
1. Lab environment is relatively clean
2. Energy detection doesn't require matched filtering
3. Simplicity and low latency are priorities
4. TX signal is already bandwidth-limited by RRC at transmitter

**Optional filtering** available via `--filter` for:
- Noisy environments
- Presence of interferers
- Testing different configurations

---

### 3. Energy Detector Design

**Implementation:** IIR-based smoothing power detector

#### Algorithm:
```
P[n] = |I[n]|² + |Q[n]|²                    // Instantaneous power
y[n] = α·y[n-1] + (1-α)·P[n]                // IIR smoothing
trigger = (y[n] > threshold)                // Detection
```

**Key Parameters:**
- **α (alpha):** IIR smoothing factor
  - Range: 0 ≤ α < 1
  - Default: 0.9
  - Higher α = more smoothing, slower response
  - Lower α = less smoothing, faster response

- **Threshold:** Detection threshold in dB
  - Default: -20 dB
  - Adjust based on signal strength and noise level
  - Too high: Miss weak packets
  - Too low: False triggers from noise

#### Advantages:
- ✅ Simple, efficient implementation
- ✅ Configurable smoothing (α parameter)
- ✅ Handles varying signal strengths
- ✅ Low computational cost (2 multiplies, 1 add per sample)

---

## Usage Examples

### Basic Usage (Default Settings)
```bash
# No filtering, 1 MHz sampling
./receiver_bob --fc 2.412e9 --gain 40 --threshold -20
```

### With Low-Pass Filtering (Noisy Environment)
```bash
# Enable filtering to reject out-of-band noise
./receiver_bob --filter --threshold -25 --gain 40
```

### Save Packets for Analysis
```bash
# Save detected packets to file
./receiver_bob --save rx_packets.dat --threshold -20

# Then plot spectrum
python plot_spectrum.py rx_packets.dat 1e6
```

### Custom Parameters
```bash
# Adjust IIR smoothing and threshold
./receiver_bob --alpha 0.85 --threshold -22 --gain 45
```

### N210 USRP
```bash
# Use N210 instead of B200
./receiver_bob --dev type=usrp2 --fc 2.412e9 --gain 45
```

---

## Performance Considerations

### Computational Complexity

**Without filtering:**
- Per sample: 3 operations (|I|², |Q|², IIR update)
- At 1 MHz: 1M operations/sec
- **Very low CPU usage**

**With filtering (101 taps):**
- Per sample: 101 multiplies + 100 adds + energy detection
- At 1 MHz: ~100M operations/sec
- **Moderate CPU usage**

**Comparison to matched RRC:**
- Matched RRC would require ~2x filter taps for proper implementation
- Plus additional complexity for timing synchronization
- **Not justified for energy detection**

### Memory Requirements

**Buffer sizes:**
- Read buffer: 10,000 samples × 8 bytes = 80 KB
- Packet buffer: 1,000 samples × 8 bytes = 8 KB
- Filter state (if enabled): 101 taps × 8 bytes = 808 bytes
- **Total: < 100 KB**

### Latency

**Without filtering:**
- Detection latency ≈ 1/α samples
- At α=0.9: ~10 samples = 10 μs
- **Very low latency**

**With filtering:**
- Additional latency = filter group delay ≈ 50 samples
- Total latency ≈ 60 μs
- **Still acceptable**

---

## Summary

| Parameter | Choice | Rationale |
|-----------|--------|-----------|
| **RX Sample Rate** | 1 MHz | Matches TX, satisfies Nyquist, no resampling |
| **Filtering** | None (default) | Sufficient for clean channel, lowest latency |
| **Optional Filter** | 101-tap LP | Available for noisy environments |
| **Energy Detector** | IIR smoothing | Simple, efficient, configurable |
| **Threshold** | -20 dB | Adjustable based on conditions |
| **Alpha** | 0.9 | Good balance of smoothing and response |

---

## Expected Performance

### In Clean Channel (Lab Environment)
- **Packet detection rate:** ~100% (with proper threshold)
- **False alarm rate:** Very low (< 0.1%)
- **Latency:** ~10 μs
- **CPU usage:** < 5%

### With Interferers
- **Without filtering:** May detect interferers as packets
- **With filtering:** Better rejection of out-of-band signals
- **Recommendation:** Use `--filter` flag

---

## Troubleshooting

### Problem: No packets detected
**Solutions:**
1. Lower threshold: `--threshold -25` or `-30`
2. Increase RX gain: `--gain 50`
3. Check TX is transmitting
4. Verify center frequency matches

### Problem: Too many false triggers
**Solutions:**
1. Raise threshold: `--threshold -15` or `-10`
2. Enable filtering: `--filter`
3. Increase IIR smoothing: `--alpha 0.95`

### Problem: Packets detected but power too low
**Solutions:**
1. Increase RX gain: `--gain 50`
2. Increase TX gain at transmitter
3. Move antennas closer

---

## Comparison to Lab 2

| Feature | Lab 2 | Lab 3 |
|---------|-------|-------|
| **Sampling** | Variable (resampling) | Fixed 1 MHz (no resampling) |
| **Filtering** | Polyphase LP filter | Optional simple LP filter |
| **Complexity** | Higher (rational resampling) | Lower (direct processing) |
| **Use Case** | General-purpose receiver | Optimized for 1 MHz BPSK |
| **Flexibility** | Multiple sample rates | Single optimal rate |

**Lab 3 is simpler** because we know the exact TX characteristics and can optimize accordingly.

---

## References

- TX parameters from `transmitter_alice.cpp`
- Energy detection theory: Urkowitz (1967)
- IIR filter design: Oppenheim & Schafer

