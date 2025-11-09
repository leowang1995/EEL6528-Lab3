#!/usr/bin/env python3
"""
plot_spectrum_trans.py - Plot spectrum of transmitted samples with comprehensive analysis

Reads complex float32 samples from binary file saved by transmitter_usrp.cpp
and plots the power spectral density (PSD) to verify bandwidth ≤ 1 MHz.

Features:
- PSD analysis with bandwidth verification
- Passband/stopband analysis
- Comprehensive statistics
- Dual plot view (full + zoomed)

Usage:
    python plot_spectrum_trans.py [filename] [sample_rate]

Arguments:
    filename     : Binary file containing TX samples (default: tx_samples.bin)
    sample_rate  : TX sample rate in Hz (default: 1e6 = 1 MHz)

Example:
    python plot_spectrum_trans.py tx_samples.bin 1e6
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import sys
import os

# ==== Default Parameters ====
DEFAULT_FILENAME = 'tx_samples.bin'
DEFAULT_FS = 1e6  # 1 MHz (TX sample rate)
DEFAULT_NFFT = 2048  # FFT size for PSD calculation

# ==== Parse Command-Line Arguments ====
filename = DEFAULT_FILENAME
fs = DEFAULT_FS

if len(sys.argv) >= 2:
    filename = sys.argv[1]
if len(sys.argv) >= 3:
    fs = float(sys.argv[2])

print("=" * 70)
print("TX Spectrum Plotter")
print("=" * 70)
print(f"Reading file: {filename}")
print(f"Sample rate: {fs/1e6:.3f} MHz")

# ==== Load IQ Samples ====
if not os.path.exists(filename):
    print(f"Error: File '{filename}' not found!")
    print("\nTo generate TX samples, run transmitter with --save option:")
    print("  ./transmitter_usrp --save tx_samples.bin --mod bpsk")
    sys.exit(1)

try:
    # Read complex float32 (8 bytes per sample: 4 for I, 4 for Q)
    data = np.fromfile(filename, dtype=np.complex64)
    print(f"Loaded {len(data):,} complex samples")
    print(f"Duration: {len(data)/fs:.3f} seconds")
except Exception as e:
    print(f"Error reading file: {e}")
    sys.exit(1)

if len(data) == 0:
    print("Error: File is empty!")
    sys.exit(1)

# ==== Calculate Power Spectral Density ====
# Adjust nfft if data is too short
nfft = DEFAULT_NFFT
if len(data) < nfft:
    nfft = 2 ** int(np.log2(len(data)))  # Next lower power of 2
    print(f"Adjusted FFT size to {nfft} (data too short for original size)")

# Use Welch's method for PSD estimation
window = signal.windows.hann(nfft, sym=False)
f, pxx = signal.welch(data, fs=fs, window=window, 
                      nperseg=nfft, noverlap=nfft//2,
                      nfft=nfft, return_onesided=False,
                      scaling='density')

# Shift zero frequency to center
f = np.fft.fftshift(f)
pxx = np.fft.fftshift(pxx)

# Convert to dB
pxx_db = 10 * np.log10(pxx + 1e-20)  # Add small value to avoid log(0)

# ==== Bandwidth Analysis ====
# Find peak frequency
peak_idx = np.argmax(pxx_db)
peak_freq = f[peak_idx]
peak_power_db = pxx_db[peak_idx]

# Calculate -3dB bandwidth
threshold_3db = peak_power_db - 3.0
# Find lower -3dB point
lower_idx = peak_idx
while lower_idx > 0 and pxx_db[lower_idx] > threshold_3db:
    lower_idx -= 1
# Find upper -3dB point
upper_idx = peak_idx
while upper_idx < len(pxx_db) - 1 and pxx_db[upper_idx] > threshold_3db:
    upper_idx += 1

bw_3db = abs(f[upper_idx] - f[lower_idx])

# Calculate occupied bandwidth (99% power)
# Convert PSD to linear scale for power calculation
pxx_linear = 10**(pxx_db / 10)
total_power = np.sum(pxx_linear) * (f[1] - f[0])  # Approximate integration
target_power = 0.99 * total_power

# Sort frequencies by power (descending)
power_freq_pairs = list(zip(pxx_linear, f))
power_freq_pairs.sort(reverse=True, key=lambda x: x[0])

# Find frequencies containing 99% of power
accumulated_power = 0.0
min_freq = 1e9
max_freq = -1e9
for power, freq in power_freq_pairs:
    accumulated_power += power * (f[1] - f[0])  # Approximate integration
    min_freq = min(min_freq, freq)
    max_freq = max(max_freq, freq)
    if accumulated_power >= target_power:
        break

bw_occupied_99 = max_freq - min_freq

# Check if bandwidth is within specification (≤ 1 MHz)
bw_within_spec = bw_occupied_99 <= 1e6

# ==== Passband/Stopband Analysis (from plot_spectrum.py) ====
# Calculate power in passband vs stopband
# Passband edge is at ±fs/4 (Nyquist criterion for baseband)
passband_edge = fs / 4
passband_mask = np.abs(f) <= passband_edge
passband_power_db = np.mean(pxx_db[passband_mask])
stopband_power_db = np.mean(pxx_db[~passband_mask])
passband_stopband_ratio_db = passband_power_db - stopband_power_db

# ==== Plot ====
plt.figure(figsize=(14, 8))

# Main PSD plot
plt.subplot(2, 1, 1)
plt.plot(f/1e6, pxx_db, linewidth=1.5, color='blue', label='PSD')
plt.grid(True, alpha=0.3)
plt.xlabel('Frequency (MHz)', fontsize=12)
plt.ylabel('Power Spectral Density (dB)', fontsize=12)
plt.title(f'TX Signal Spectrum (fs={fs/1e6:.3f} MHz, {len(data):,} samples)', fontsize=14, fontweight='bold')

# Add vertical lines for bandwidth specification
plt.axvline(-0.5, color='red', linestyle='--', alpha=0.7, linewidth=2, label='±0.5 MHz (1 MHz BW limit)')
plt.axvline(0.5, color='red', linestyle='--', alpha=0.7, linewidth=2)

# Add passband edges (±fs/4) - from plot_spectrum.py
plt.axvline(-passband_edge/1e6, color='purple', linestyle=':', alpha=0.6, linewidth=1.5, 
            label=f'±{passband_edge/1e6:.3f} MHz (fs/4 passband edge)')
plt.axvline(passband_edge/1e6, color='purple', linestyle=':', alpha=0.6, linewidth=1.5)

# Mark measured bandwidth
plt.axvline(-bw_occupied_99/2e6, color='green', linestyle=':', alpha=0.7, linewidth=2, 
            label=f'99% Occupied BW: ±{bw_occupied_99/2e6:.3f} MHz')
plt.axvline(bw_occupied_99/2e6, color='green', linestyle=':', alpha=0.7, linewidth=2)

# Mark -3dB bandwidth
if bw_3db > 0:
    plt.axvline(-bw_3db/2e6, color='orange', linestyle='--', alpha=0.6, linewidth=1.5,
                label=f'-3dB BW: {bw_3db/1e3:.1f} kHz')
    plt.axvline(bw_3db/2e6, color='orange', linestyle='--', alpha=0.6, linewidth=1.5)

plt.axhline(threshold_3db, color='orange', linestyle='--', alpha=0.5, label='-3dB threshold')
plt.legend(fontsize=9, loc='upper right')
plt.xlim(-fs/2e6, fs/2e6)

# Zoomed view around center frequency
plt.subplot(2, 1, 2)
zoom_range = 2.0  # MHz
plt.plot(f/1e6, pxx_db, linewidth=1.5, color='blue')
plt.grid(True, alpha=0.3)
plt.xlabel('Frequency (MHz)', fontsize=12)
plt.ylabel('Power Spectral Density (dB)', fontsize=12)
plt.title(f'Zoomed View (±{zoom_range} MHz around center)', fontsize=12)

# Add all markers in zoomed view
plt.axvline(-0.5, color='red', linestyle='--', alpha=0.7, linewidth=2)
plt.axvline(0.5, color='red', linestyle='--', alpha=0.7, linewidth=2)
plt.axvline(-passband_edge/1e6, color='purple', linestyle=':', alpha=0.6, linewidth=1.5)
plt.axvline(passband_edge/1e6, color='purple', linestyle=':', alpha=0.6, linewidth=1.5)
plt.axvline(-bw_occupied_99/2e6, color='green', linestyle=':', alpha=0.7, linewidth=2)
plt.axvline(bw_occupied_99/2e6, color='green', linestyle=':', alpha=0.7, linewidth=2)
if bw_3db > 0:
    plt.axvline(-bw_3db/2e6, color='orange', linestyle='--', alpha=0.6, linewidth=1.5)
    plt.axvline(bw_3db/2e6, color='orange', linestyle='--', alpha=0.6, linewidth=1.5)

plt.xlim(-zoom_range, zoom_range)

plt.tight_layout()

# Save figure
output_file = filename.replace('.bin', '_spectrum.png')
if output_file == filename:  # If no .bin extension
    output_file = filename + '_spectrum.png'
plt.savefig(output_file, dpi=150, bbox_inches='tight')
print(f"\nSpectrum plot saved to: {output_file}")

plt.show()

# ==== Print Statistics ====
print("\n" + "=" * 70)
print("SPECTRUM ANALYSIS RESULTS")
print("=" * 70)
print(f"Sample rate:         {fs/1e6:.3f} MHz")
print(f"FFT size:            {nfft}")
print(f"Frequency resolution: {fs/nfft/1e3:.3f} kHz")
print()
print(f"Peak frequency:      {peak_freq/1e3:.3f} kHz")
print(f"Peak power:          {peak_power_db:.2f} dB")
print()

# ==== Basic Spectrum Statistics (from plot_spectrum.py) ====
print("Basic Spectrum Statistics:")
print(f"  Peak power:        {np.max(pxx_db):.2f} dB")
print(f"  Mean power:        {np.mean(pxx_db):.2f} dB")
print(f"  Median power:      {np.median(pxx_db):.2f} dB")
print()

# ==== Bandwidth Measurements ====
print("Bandwidth Measurements:")
print(f"  -3dB bandwidth:    {bw_3db/1e3:.3f} kHz")
print(f"  99% occupied bandwidth: {bw_occupied_99/1e3:.3f} kHz")
print()

# ==== Bandwidth Specification Check ====
print("Bandwidth Specification Check:")
print(f"  Required:          ≤ 1000 kHz (1 MHz)")
print(f"  Measured (99%):    {bw_occupied_99/1e3:.3f} kHz")
if bw_within_spec:
    margin = 1e6 - bw_occupied_99
    print(f"  Status:            ✓ PASS - Bandwidth within specification")
    print(f"  Margin:            {margin/1e3:.3f} kHz")
else:
    excess = bw_occupied_99 - 1e6
    print(f"  Status:            ✗ FAIL - Bandwidth exceeds specification!")
    print(f"  Excess:            {excess/1e3:.3f} kHz")
    print("  WARNING: TX bandwidth exceeds 1 MHz specification!")
    print("  Consider reducing symbol rate or excess bandwidth.")
print()

# ==== Passband/Stopband Analysis (from plot_spectrum.py) ====
print("Passband/Stopband Analysis:")
print(f"  Passband edge:     ±{passband_edge/1e6:.3f} MHz (fs/4)")
print(f"  Average power in passband (±{passband_edge/1e6:.3f} MHz): {passband_power_db:.2f} dB")
print(f"  Average power in stopband: {stopband_power_db:.2f} dB")
print(f"  Passband/Stopband ratio: {passband_stopband_ratio_db:.2f} dB")
print()

# ==== Additional Statistics ====
print("Additional Statistics:")
print(f"  Mean PSD:           {np.mean(pxx_db):.2f} dB")
print(f"  Median PSD:         {np.median(pxx_db):.2f} dB")
print(f"  Std dev PSD:        {np.std(pxx_db):.2f} dB")
print(f"  Min PSD:            {np.min(pxx_db):.2f} dB")
print(f"  Max PSD:            {np.max(pxx_db):.2f} dB")
print("=" * 70)
