[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/6TkJhj8O)
# EEL6528 Lab 3: Let's shape the spectrum

## Goals:
- You will practice generating a linearly modulated TX signal by pulse shaping as discussed in class.
- You will use the simple energy detector in Lab 2 to capture TX packets at the RX as well as implement an AGC to normalize the signal level of the captured packets as discussed in class.
- Completing this lab, you will get a taste of what are needed to continuously generate (TX) and capture (RX) packets using a pair of USRP radios.


## What you need to do:
- ***You will need to use two USRP radios in this lab: one for TX and one for RX.*** You may write two programs: one for TX and one for RX. You may also write a single program with a PO switch to select between running as TX or RX. You may organize your program(s) into any classes, functions, and threads as you see fit. 

- You may want to set both the TX and RX center frequency $f_c$ to $2.412$ GHz. It seems that that the WiFi channel at $2.412$ GHz is relatively clean where the radios are placed. You may change the value of $f_c$ to match the [frequencies of the WiFi channels](https://en.wikipedia.org/wiki/List_of_WLAN_channels) to find a relatively unused WiFi channel in the vacinity of the radios.

### Signal generation:
- Use 1 MHz TX sampling rate to achieve a symbol rate of 800k symbols per second (Sps) while making sure that the bandwidth of the TX signal is limited to 1 MHz:
  1. Generate a sequence/array/vector/block of 1000 random bits and insert this sequence of bits into a FIFO queue every second.
  2. Implement a modulator (or multiple modulators in multiple threads) that grabs a block of bits from the FIFO queue, uses the bits to BPSK-modulate a sequence of RRC pulses with a symbol rate of 800 kSps.
     > BPSK symbols are in the set {+1, -1}. The mapping from bit to symbol may be arbitrary, e.g., $0 \mapsto +1$ and $1 \mapsto -1$.
  3. Send the samples of the BPSK-modulated signal to the TX USRP for transmission.

- You may use the TX pulse shaping design in class to select the excess bandwidth of the RRC pulse. You may also use the code provided to generate the impulse response of the pulse shaping filter.

- Recall that UHD assumes the sample values (the real and imaginary parts of each sample) lie within the range $[-1,1]$, which is mapped to the full dynamic range of the DAC. UHD clips any value outside of the range to the closest boundary. Such clipping if applied often will cause nonlinear distortions on the TX signal. Hence, it is important to prevent excessive clipping. In our case, passing the BPSK symbols through the RRC pulse-shaping filter will likely to generate sample values outside $[-1,1]$. A simple method to prevent clipping is to scale down the filter output signal by either reducing the filter gain or explicitly multiplying the output samples with a scalar constant smaller than one. The tradeoff is clearly a reduction in TX power (and perhaps also increasing quantization noise) if too much scaling down is applied. You may determine how much scaling is required by generating a histogram of the sample values at the output of the pulse-shaping filter.

- You should check the PSD of the TX signal in order to verify that its bandwidth is within 1 MHz. You may do so by capturing the TX samples sent to the USRP radio. Alternatively, you may repeat the steps in Lab 2 capturing the TX signal using the RX radio and then saving the captured samples to calculate and plot the spectrum (see below). 

## Signal capture:
- Implement a simple energy detector as you did in Lab 2 to capture the TX packets. You may choose an appropriate value for the RX sampling rate. You may apply any filter to the RX samples as you deem appropriate before feeding them into the energy detector. **Explain your choice of the RX sampling rate and your choice of any filtering (or not filtering).**
- The detection threshold should be chosen to detect most of the TX packets with a minimal false alarm rate. For each TX packet, you should capture all the samples corresponding to the whole packet. In fact, the number of samples to capture should be a bit larger than the number corresponding to a whole packet, with extra samples captured at both the beginning and end of the packet, to allow for proper synchronization later (no need to implement synchronization that in this lab). Remember to handle the case when the TX packet lies across the boundary between two adjacent blocks of RX samples that you grab from the USRP radio.
- Implement an AGC to normalize the magnitude level of the captured packet samples to 1 (the BPSK symbol signal magnitude). There is no need to implement a closed-loop AGC. It suffices to simply normalize the RMS magnitude of the whole block of samples to 1. Push the AGCed packet to a FIFO queue.
- Implement a function that grabs the packet samples from the FIFO queue and counts and prints the number of packets captured every 10 seconds. You should be able to use the number of packets captured every 10 seconds to check whether your energy detector is designed properly.
- Save the samples of the captured packets. Use them to calculate and plot the PSD of the TX signal. Check to make sure that the TX signal's bandwidth is no larger than 1 MHz. Note that the two USRP radios are placed only 3 feet apart in the lab. You may need to lower your RX gain if you see spurious spectral components in your PSD plot. If the TX and RX gains are too high, the RX signal may be clipped at the input of the ADC at the RX radio as we discussed in class. The clipping (nonlinear distortion) may create spurious spectral components.
 

## More experiments:
- Increase the packet generation rate at the TX to beyond 1 packet per second. Since the symbol rate is 800 kSps and each packet contains 1000 bits (symbols), the maximum possible packet generation rate (without overflowing the FIFO queue) is 800 packets per second or 1 packet per 1.25ms.
- Observe and report the length of the TX FIFO queue(s) as the packet generation rate increases.
- Observe and report the number of packets captured at the RX as the packet generation rate increases.
- Increase the TX sampling rate up to 12.5 MHz (check Lab 2 for the list of "good" sampling rates to use) while keeping the symbol rate to be 80% of the sampling rate and the bandwidth of the TX signal no larger than the sampling rate. Repeat the previous experiment by increasing the packet generation rate. What's the maximum symbol rate you can attain without underflowing at the TX, overflowing at the RX, or overflowing any queue you implement at the TX or RX? 
