# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 12:15:59 2024

@author: Ale Campoy
"""
import numpy as np
import matplotlib.pyplot as plt
import pywt

# Parameters
frequency = 1  # Frequency of the pulses (in Hz)
pulse_width = 0.2  # Width of each pulse (in seconds)
amplitude = 1  # Amplitude of the pulses
duration = 5  # Duration of the signal (in seconds)
sampling_rate = 1000  # Sampling rate (in Hz)

# Time array
t = np.linspace(0, duration, int(duration * sampling_rate), endpoint=False)

# Generating the original signal
original_signal = np.sin(2 * np.pi * frequency * t)

# Perform wavelet decomposition (using Daubechies wavelet)
coefficients = pywt.wavedec(original_signal, "db4")

# Reconstruct the signal using all coefficients
reconstructed_signal = pywt.waverec(coefficients, "db4")[: len(t)]

# Threshold the reconstructed signal to obtain pulsed representation
threshold = 0.1  # Adjust the threshold as needed
pulsed_signal = np.where(reconstructed_signal > threshold, amplitude, 0)

# Plotting the original and pulsed signals
plt.figure(figsize=(10, 6))
plt.plot(t, original_signal, label="Original Signal")
plt.plot(t, pulsed_signal, label="Pulsed Signal", linestyle="--")
plt.title("Original and Pulsed Signal")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.legend()
plt.grid(True)
plt.show()

# %% example 2
import numpy as np
import matplotlib.pyplot as plt
import pywt

# Parameters
frequency = 1  # Frequency of the pulses (in Hz)
pulse_width = 0.2  # Width of each pulse (in seconds)
amplitude = 1  # Amplitude of the pulses
duration = 5  # Duration of the signal (in seconds)
sampling_rate = 1000  # Sampling rate (in Hz)

# Time array
t = np.linspace(0, duration, int(duration * sampling_rate), endpoint=False)

# Generating the original signal
original_signal = np.sin(2 * np.pi * frequency * t)
original_signal += 0.2 * np.sin(2 * np.pi * 2.2 * t)
original_signal += 0.5 * np.sin(2 * np.pi * 3 * t)

# Perform Discrete Wavelet Transform (DWT)
level = 5  # Adjust the level of decomposition as needed
coefficients = pywt.wavedec(original_signal, "db4", level=level)

# Reconstruct the signal using all coefficients
reconstructed_signal = pywt.waverec(coefficients, "db4")[: len(t)]

# Threshold the reconstructed signal to obtain pulsed representation
threshold = 0.1  # Adjust the threshold as needed
pulsed_signal = np.where(reconstructed_signal > threshold, amplitude, 0)

# Plotting the original and pulsed signals
plt.figure(figsize=(10, 6))
plt.plot(t, original_signal, label="Original Signal")
plt.plot(t, pulsed_signal, label="Pulsed Signal", linestyle="--")
plt.title("Original and Pulsed Signal")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.legend()
plt.grid(True)
plt.show()
