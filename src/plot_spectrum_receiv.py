#!/usr/bin/env python3
"""
plot_spectrum_receiv.py - Plot spectrum of received samples with comprehensive analysis

Reads complex float32 samples from binary file saved by receiver_usrp.cpp
and plots the power spectral density (PSD) to verify bandwidth ≤ 1 MHz.

Features:
- PSD analysis with bandwidth verification
- ADC clipping detection
- AGC verification (RMS magnitude check)
- Time domain visualization
- Comprehensive statistics

The receiver captures packets with:
- 1700 samples per packet (100 pre-trigger + 1500 data + 100 post-trigger)
- 800 kHz sample rate (after resampling)
- AGC normalized to RMS magnitude = 1

Usage:
    python plot_spectrum_receiv.py <filename> <sample_rate> [options]

Examples:
    python plot_spectrum_receiv.py samples_iq.bin 800e3
    python plot_spectrum_receiv.py samples_iq.bin 800e3 --nfft 4096
    python plot_spectrum_receiv.py samples_iq.bin 800e3 --output my_plot.png
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import sys
import os
import argparse

def load_samples(filename):
    """Load complex float32 samples from binary file"""
    if not os.path.exists(filename):
        print(f"✗ Error: File '{filename}' not found!")
        print("\nTo generate RX samples, run receiver with --save option:")
        print("  ./receiver_usrp --save samples_iq.bin --fs_out 800e3")
        sys.exit(1)
    
    try:
        data = np.fromfile(filename, dtype=np.complex64)
        print(f"✓ Loaded {len(data):,} complex samples from {filename}")
        return data
    except Exception as e:
        print(f"✗ Error loading file: {e}")
        sys.exit(1)

def compute_psd(samples, fs, nfft=2048, nperseg=None):
    """
    Compute Power Spectral Density using Welch's method
    
    Parameters:
    - samples: Complex IQ samples
    - fs: Sample rate in Hz
    - nfft: FFT size
    - nperseg: Segment length for Welch's method
    
    Returns:
    - f: Frequency array (Hz)
    - psd_db: PSD in dB
    """
    # Adjust nfft if data is too short
    if len(samples) < nfft:
        nfft = 2 ** int(np.log2(len(samples)))  # Next lower power of 2
        print(f"  Adjusted FFT size to {nfft} (data too short for original size)")
    
    if nperseg is None:
        nperseg = min(nfft, len(samples))
    
    # Use Welch's method for PSD estimation
    f, pxx = signal.welch(
        samples,
        fs=fs,
        window='hann',
        nperseg=nperseg,
        noverlap=nperseg//2,
        nfft=nfft,
        return_onesided=False,
        scaling='density'
    )
    
    # Shift zero frequency to center
    f = np.fft.fftshift(f)
    pxx = np.fft.fftshift(pxx)
    
    # Convert to dB
    psd_db = 10 * np.log10(pxx + 1e-20)
    
    return f, psd_db, nfft

def measure_bandwidth(f, psd_db, method='3dB', percentage=0.99):
    """
    Measure signal bandwidth
    
    Parameters:
    - f: Frequency array
    - psd_db: PSD in dB
    - method: '3dB' or 'occupied'
    - percentage: For occupied bandwidth (default 0.99 = 99%)
    
    Returns:
    - bandwidth: Bandwidth in Hz
    - lower_freq: Lower frequency bound
    - upper_freq: Upper frequency bound
    """
    if method == '3dB':
        # Find peak power
        peak_power = np.max(psd_db)
        threshold = peak_power - 3.0  # -3 dB point
        
        # Find frequencies where power > threshold
        mask = psd_db > threshold
        freqs_above_threshold = f[mask]
        
        if len(freqs_above_threshold) > 0:
            lower_freq = np.min(freqs_above_threshold)
            upper_freq = np.max(freqs_above_threshold)
            bandwidth = upper_freq - lower_freq
        else:
            lower_freq = upper_freq = bandwidth = 0
            
    elif method == 'occupied':
        # Convert PSD from dB to linear
        pxx_linear = 10 ** (psd_db / 10.0)
        
        # Sort by power (descending)
        sorted_indices = np.argsort(pxx_linear)[::-1]
        sorted_freqs = f[sorted_indices]
        sorted_power = pxx_linear[sorted_indices]
        
        # Find frequencies containing specified percentage of power
        total_power = np.sum(sorted_power)
        target_power = percentage * total_power
        
        accumulated_power = 0
        freqs_in_band = []
        
        for i, power in enumerate(sorted_power):
            accumulated_power += power
            freqs_in_band.append(sorted_freqs[i])
            if accumulated_power >= target_power:
                break
        
        if len(freqs_in_band) > 0:
            lower_freq = np.min(freqs_in_band)
            upper_freq = np.max(freqs_in_band)
            bandwidth = upper_freq - lower_freq
        else:
            lower_freq = upper_freq = bandwidth = 0
    
    else:
        raise ValueError(f"Unknown method: {method}")
    
    return bandwidth, lower_freq, upper_freq

def check_for_clipping(samples):
    """
    Check for potential ADC clipping
    
    Clipping indicators:
    - Samples at or near ±1.0 (ADC full scale)
    - High Peak-to-Average Power Ratio (PAPR)
    - Distorted distribution
    
    Returns:
    - is_clipped: Boolean
    - stats: Dictionary of statistics
    """
    magnitudes = np.abs(samples)
    
    # Check for samples near full scale (±0.95)
    near_clipping = np.sum(magnitudes > 0.95)
    clipping_percentage = 100.0 * near_clipping / len(samples)
    
    # Calculate PAPR
    peak_power = np.max(magnitudes**2)
    avg_power = np.mean(magnitudes**2)
    papr_db = 10 * np.log10(peak_power / avg_power) if avg_power > 0 else 0
    
    # Calculate statistics
    stats = {
        'max_magnitude': np.max(magnitudes),
        'mean_magnitude': np.mean(magnitudes),
        'std_magnitude': np.std(magnitudes),
        'samples_near_clipping': near_clipping,
        'clipping_percentage': clipping_percentage,
        'papr_db': papr_db
    }
    
    # Clipping detected if:
    # 1. More than 1% of samples near full scale, OR
    # 2. Any samples at exactly 1.0 (hard clipping)
    is_clipped = (clipping_percentage > 1.0) or (np.max(magnitudes) >= 0.99)
    
    return is_clipped, stats

def verify_agc(samples):
    """
    Verify AGC normalization
    
    After AGC, RMS magnitude should be ~1.0 for BPSK symbols
    
    Returns:
    - rms_magnitude: RMS magnitude
    - mean_power: Mean power
    - is_normalized: Boolean indicating if AGC is correct
    """
    rms_magnitude = np.sqrt(np.mean(np.abs(samples)**2))
    mean_power = np.mean(np.abs(samples)**2)
    is_normalized = abs(rms_magnitude - 1.0) < 0.1
    
    return rms_magnitude, mean_power, is_normalized

def plot_spectrum(f, psd_db, fs, samples, bandwidth_3db, bandwidth_99, 
                  clipping_detected, output_file='rx_spectrum.png'):
    """
    Plot PSD spectrum with time domain view
    
    Creates 3 subplots:
    1. Full PSD spectrum
    2. Zoomed PSD view
    3. Time domain (first packet)
    """
    fig = plt.figure(figsize=(14, 12))
    
    # ==== Main PSD Plot ====
    ax1 = plt.subplot(3, 1, 1)
    ax1.plot(f/1e6, psd_db, linewidth=1.5, color='blue', label='PSD')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlabel('Frequency (MHz)', fontsize=12)
    ax1.set_ylabel('Power Spectral Density (dB)', fontsize=12)
    
    # Title with clipping warning if detected
    if clipping_detected:
        title = 'RX Signal Spectrum - ⚠ CLIPPING DETECTED'
        ax1.set_title(title, fontsize=14, color='red', fontweight='bold')
    else:
        title = f'RX Signal Spectrum (fs={fs/1e6:.3f} MHz, {len(samples):,} samples)'
        ax1.set_title(title, fontsize=14, fontweight='bold')
    
    # Mark 1 MHz bandwidth limit
    ax1.axvline(-0.5, color='red', linestyle='--', linewidth=2, 
                alpha=0.7, label='±0.5 MHz (1 MHz BW limit)')
    ax1.axvline(0.5, color='red', linestyle='--', linewidth=2, alpha=0.7)
    
    # Mark measured bandwidth
    ax1.axvline(-bandwidth_99/2e6, color='green', linestyle=':', linewidth=2,
                alpha=0.7, label=f'99% Occupied BW: ±{bandwidth_99/2e6:.3f} MHz')
    ax1.axvline(bandwidth_99/2e6, color='green', linestyle=':', linewidth=2, alpha=0.7)
    
    # Mark -3dB bandwidth
    if bandwidth_3db > 0:
        ax1.axvline(-bandwidth_3db/2e6, color='orange', linestyle='--', linewidth=1.5,
                    alpha=0.6, label=f'-3dB BW: {bandwidth_3db/1e3:.1f} kHz')
        ax1.axvline(bandwidth_3db/2e6, color='orange', linestyle='--', linewidth=1.5, alpha=0.6)
    
    ax1.legend(loc='upper right', fontsize=10)
    ax1.set_xlim([-fs/2e6, fs/2e6])
    
    # ==== Zoomed PSD View ====
    ax2 = plt.subplot(3, 1, 2)
    zoom_range = 1.0  # MHz
    mask = (f >= -zoom_range*1e6) & (f <= zoom_range*1e6)
    
    ax2.plot(f[mask]/1e6, psd_db[mask], linewidth=1.5, color='blue')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlabel('Frequency (MHz)', fontsize=12)
    ax2.set_ylabel('Power Spectral Density (dB)', fontsize=12)
    ax2.set_title(f'Zoomed View (±{zoom_range} MHz around center)', fontsize=12)
    
    # Mark bandwidth limits in zoomed view
    ax2.axvline(-0.5, color='red', linestyle='--', linewidth=2, alpha=0.7)
    ax2.axvline(0.5, color='red', linestyle='--', linewidth=2, alpha=0.7)
    ax2.axvline(-bandwidth_99/2e6, color='green', linestyle=':', linewidth=2, alpha=0.7)
    ax2.axvline(bandwidth_99/2e6, color='green', linestyle=':', linewidth=2, alpha=0.7)
    if bandwidth_3db > 0:
        ax2.axvline(-bandwidth_3db/2e6, color='orange', linestyle='--', linewidth=1.5, alpha=0.6)
        ax2.axvline(bandwidth_3db/2e6, color='orange', linestyle='--', linewidth=1.5, alpha=0.6)
    
    ax2.set_xlim([-zoom_range, zoom_range])
    
    # ==== Time Domain Plot ====
    ax3 = plt.subplot(3, 1, 3)
    samples_per_packet = 1700
    if len(samples) >= samples_per_packet:
        # Plot first packet
        first_packet = samples[:samples_per_packet]
        time_axis = np.arange(len(first_packet)) / fs * 1e6  # microseconds
        ax3.plot(time_axis, np.real(first_packet), linewidth=1, alpha=0.7, 
                label='I (real)', color='blue')
        ax3.plot(time_axis, np.imag(first_packet), linewidth=1, alpha=0.7, 
                label='Q (imag)', color='red')
        ax3.set_title(f'Time Domain: First Packet ({samples_per_packet} samples)', fontsize=12)
    else:
        # Plot all available samples
        time_axis = np.arange(len(samples)) / fs * 1e6  # microseconds
        ax3.plot(time_axis, np.real(samples), linewidth=1, alpha=0.7, 
                label='I (real)', color='blue')
        ax3.plot(time_axis, np.imag(samples), linewidth=1, alpha=0.7, 
                label='Q (imag)', color='red')
        ax3.set_title(f'Time Domain: All Samples ({len(samples)} samples)', fontsize=12)
    
    ax3.set_xlabel('Time (μs)', fontsize=12)
    ax3.set_ylabel('Amplitude', fontsize=12)
    ax3.grid(True, alpha=0.3)
    ax3.legend(fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=200, bbox_inches='tight')
    print(f"✓ Spectrum plot saved to: {output_file}")
    
    return fig

def print_analysis(samples, fs, nfft, bandwidth_3db, bandwidth_99, 
                  clipping_stats, agc_stats):
    """
    Print comprehensive analysis results
    """
    print("\n" + "="*70)
    print("SPECTRUM ANALYSIS RESULTS")
    print("="*70)
    
    print(f"\nSample Statistics:")
    print(f"  Total samples:     {len(samples):,}")
    print(f"  Sample rate:       {fs/1e6:.3f} MHz")
    print(f"  Duration:          {len(samples)/fs*1e3:.2f} ms")
    print(f"  FFT size:          {nfft}")
    print(f"  Frequency resolution: {fs/nfft/1e3:.3f} kHz")
    
    # Estimate number of packets
    samples_per_packet = 1700
    num_packets = len(samples) // samples_per_packet
    if num_packets > 0:
        print(f"  Estimated packets:  {num_packets} (assuming {samples_per_packet} samples/packet)")
    
    print(f"\nSignal Magnitude:")
    print(f"  Max magnitude:     {clipping_stats['max_magnitude']:.6f}")
    print(f"  Mean magnitude:    {clipping_stats['mean_magnitude']:.6f}")
    print(f"  Std deviation:     {clipping_stats['std_magnitude']:.6f}")
    print(f"  PAPR:              {clipping_stats['papr_db']:.2f} dB")
    
    print(f"\nBandwidth Measurements:")
    print(f"  -3dB bandwidth:     {bandwidth_3db/1e3:.1f} kHz")
    print(f"  99% occupied BW:  {bandwidth_99/1e3:.1f} kHz")
    
    print(f"\nBandwidth Verification:")
    print(f"  Requirement:       ≤ 1000 kHz (1 MHz)")
    if bandwidth_99 <= 1.0e6:
        print(f"  -3dB result:       {bandwidth_3db/1e3:.1f} kHz  ✓ PASS")
        print(f"  99% result:       {bandwidth_99/1e3:.1f} kHz  ✓ PASS")
        margin = 1.0e6 - bandwidth_99
        print(f"  Margin:           {margin/1e3:.1f} kHz")
    else:
        print(f"  -3dB result:       {bandwidth_3db/1e3:.1f} kHz")
        print(f"  99% result:        {bandwidth_99/1e3:.1f} kHz  ✗ FAIL")
        excess = bandwidth_99 - 1.0e6
        print(f"  Excess:           {excess/1e3:.1f} kHz")
        print(f"  ⚠ WARNING: Bandwidth exceeds 1 MHz specification!")
    
    print(f"\nAGC Verification:")
    print(f"  RMS magnitude:     {agc_stats['rms_magnitude']:.6f} (should be ~1.0)")
    print(f"  Mean power:        {agc_stats['mean_power']:.6f} (should be ~1.0)")
    if agc_stats['is_normalized']:
        print(f"  Status:            ✓ AGC normalization appears correct")
    else:
        print(f"  Status:            ⚠ AGC normalization may need adjustment")
        print(f"                     Expected RMS ≈ 1.0 for BPSK symbols")
    
    print(f"\nClipping Analysis:")
    print(f"  Samples near clipping: {clipping_stats['samples_near_clipping']} "
          f"({clipping_stats['clipping_percentage']:.2f}%)")
    
    if clipping_stats['clipping_percentage'] > 1.0:
        print(f"  ✗ CLIPPING DETECTED!")
        print(f"     {clipping_stats['clipping_percentage']:.1f}% of samples near ADC full scale")
        print(f"     This creates spurious spectral components!")
        print(f"\n  Recommended actions:")
        print(f"     1. Lower RX gain (try --gain 30 or less)")
        print(f"     2. Lower TX gain (try --gain 20 or less)")
        print(f"     3. Increase distance between radios (>3 feet)")
        print(f"     4. Check for proper AGC normalization")
    else:
        print(f"  ✓ No significant clipping detected")
        print(f"     Signal levels appropriate for ADC range")
    
    print("="*70 + "\n")

def main():
    parser = argparse.ArgumentParser(
        description='Plot spectrum of RX samples with comprehensive analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage (800 kHz sample rate)
  python plot_spectrum_receiv.py samples_iq.bin 800e3
  
  # With 1 MHz sample rate
  python plot_spectrum_receiv.py samples_iq.bin 1e6
  
  # Custom FFT size
  python plot_spectrum_receiv.py samples_iq.bin 800e3 --nfft 4096
  
  # Custom output filename
  python plot_spectrum_receiv.py samples_iq.bin 800e3 --output my_spectrum.png
  
  # Analysis only (no plot)
  python plot_spectrum_receiv.py samples_iq.bin 800e3 --no-plot

Notes:
  - Sample file should contain complex64 IQ samples
  - Sample rate should match receiver output rate (--fs_out)
  - Check for clipping if spurious spectral components present
  - Lower RX/TX gains if clipping detected
  - AGC should normalize RMS magnitude to ~1.0
        """
    )
    
    parser.add_argument('filename', nargs='?', default='samples_iq.bin',
                       help='Binary file with complex64 samples (default: samples_iq.bin)')
    parser.add_argument('sample_rate', nargs='?', type=float, default=800e3,
                       help='Sample rate in Hz (default: 800e3 = 800 kHz)')
    parser.add_argument('--nfft', type=int, default=2048,
                       help='FFT size (default: 2048)')
    parser.add_argument('--output', default=None,
                       help='Output plot filename (default: <filename>_spectrum.png)')
    parser.add_argument('--no-plot', action='store_true',
                       help='Skip plotting (analysis only)')
    
    args = parser.parse_args()
    
    # Load samples
    print(f"\n{'='*70}")
    print(f"RX Spectrum Plotter - Loading samples...")
    print(f"{'='*70}")
    samples = load_samples(args.filename)
    
    if len(samples) == 0:
        print("✗ Error: File is empty!")
        sys.exit(1)
    
    # Check for clipping
    print(f"\nChecking for ADC clipping...")
    clipping_detected, clipping_stats = check_for_clipping(samples)
    
    # Verify AGC
    print(f"Verifying AGC normalization...")
    rms_magnitude, mean_power, is_normalized = verify_agc(samples)
    agc_stats = {
        'rms_magnitude': rms_magnitude,
        'mean_power': mean_power,
        'is_normalized': is_normalized
    }
    
    # Compute PSD
    print(f"Computing Power Spectral Density...")
    f, psd_db, nfft_actual = compute_psd(samples, args.sample_rate, nfft=args.nfft)
    
    # Measure bandwidth
    print(f"Measuring bandwidth...")
    bandwidth_3db, _, _ = measure_bandwidth(f, psd_db, method='3dB')
    bandwidth_99, _, _ = measure_bandwidth(f, psd_db, method='occupied', percentage=0.99)
    
    # Print analysis
    print_analysis(samples, args.sample_rate, nfft_actual, bandwidth_3db, 
                  bandwidth_99, clipping_stats, agc_stats)
    
    # Plot
    if not args.no_plot:
        print(f"Generating spectrum plot...")
        if args.output is None:
            output_file = args.filename.replace('.bin', '_spectrum.png')
            if output_file == args.filename:  # If no .bin extension
                output_file = args.filename + '_spectrum.png'
        else:
            output_file = args.output
        
        plot_spectrum(f, psd_db, args.sample_rate, samples, bandwidth_3db, 
                     bandwidth_99, clipping_detected, output_file=output_file)
        print(f"\n✓ Analysis complete! Check {output_file} for spectrum plot.")
    else:
        print(f"\n✓ Analysis complete!")
    
    # Return exit code based on verification
    exit_code = 0
    if bandwidth_99 > 1.0e6:
        print(f"\n⚠ WARNING: Bandwidth exceeds 1 MHz specification")
        exit_code = 1
    elif clipping_detected:
        print(f"\n⚠ WARNING: Clipping detected - check gains")
        exit_code = 2
    elif not is_normalized:
        print(f"\n⚠ WARNING: AGC normalization may need adjustment")
        exit_code = 3
    else:
        print(f"\n✓ All checks passed!")
    
    sys.exit(exit_code)

if __name__ == '__main__':
    main()
