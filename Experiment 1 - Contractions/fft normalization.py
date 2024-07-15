#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 13:01:28 2024

@author: ale
"""

import numpy as np


def aggregate_ffts(signals, sampling_rate):
    n_samples = len(signals)
    fft_length = None

    # Compute FFT for each signal
    ffts = []
    for signal in signals:
        signal -= np.mean(signal)  # Remove DC component (mean)
        signal /= np.std(signal)
        fft = np.fft.fft(signal)
        fft = fft / np.sqrt(len(fft))
        if fft_length is None:
            fft_length = len(fft)
        elif len(fft) != fft_length:
            raise ValueError("All signals must have the same length")
        ffts.append(fft)

    # Convert list to numpy array for easy averaging
    ffts = np.array(ffts)

    # Compute the average FFT
    average_fft = np.mean(ffts, axis=0)

    # Get the frequency bins
    freqs = np.fft.fftfreq(fft_length, d=1 / sampling_rate)

    return average_fft, freqs


# Example usage
sampling_rate = 1000  # in Hz

sr = sampling_rate = 500
# sampling interval
ts = 1 / sr
t = np.arange(0, 10, ts)

freq = 0.5
x = 1.5 * np.sin(2 * np.pi * freq * t)

freq = 2
x += 3 * np.sin(2 * np.pi * freq * t)


signals = [
    2 * np.sin(2 * np.pi * 0.5 * t),  # Replace with your actual signals
    1 * np.sin(2 * np.pi * 2 * t),
    6 * np.sin(2 * np.pi * 5 * t),
]

average_fft, freqs = aggregate_ffts(signals, sampling_rate)

# Optional: Compute the inverse FFT of the average spectrum
average_signal = np.fft.ifft(average_fft)

# Plotting the results (if needed)
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 5))

# Plot the magnitude of the average FFT
plt.subplot(121)
plt.plot(freqs, np.abs(average_fft))
plt.title("Average FFT Magnitude")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude")
plt.xlim(-0.1, 11)

# Plot the reconstructed average signal
plt.subplot(122)
plt.plot(np.real(average_signal))
plt.title("Reconstructed Signal from Average FFT")
plt.xlabel("Time")
plt.ylabel("Amplitude")

plt.tight_layout()
plt.show()
