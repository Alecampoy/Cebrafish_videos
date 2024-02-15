# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 19:03:59 2024

@author: Ale Campoy
"""
# https://pythonnumericalmethods.berkeley.edu/notebooks/chapter24.04-FFT-in-Python.html

from IPython import get_ipython

get_ipython().magic("reset -sf")

import matplotlib.pyplot as plt
import numpy as np

# %% FFT Example of sinousoidal
# sampling rate
sr = 200
# sampling interval
ts = 1 / sr
t = np.arange(0, 10, ts)

freq = 0.5
x = 1.5 * np.sin(2 * np.pi * freq * t)

freq = 2
x += 3 * np.sin(2 * np.pi * freq * t)

freq = 7
x += 2 * np.sin(2 * np.pi * freq * t)

x += 5
x += np.random.normal(1, 1, sr * 10)

x = x - np.mean(x)

plt.figure(figsize=(8, 6))
plt.plot(t, x, "r")
plt.ylabel("Amplitude")

plt.show()

# %%%
from scipy import fft

X = fft.fft(x, norm="backward")
N = len(x)
freq = fft.fftfreq(N, ts)

psd = np.abs(X) ** 2 / (sr * N)
# aquii esta la clave de la normalización. Por ello estudiar esto. Tambien usar la fft completa y no rfft para evitar problemas con parserval

plt.figure(figsize=(12, 6))
plt.subplot(121)

plt.stem(freq, psd, "b", markerfmt=" ", basefmt="-b")
plt.xlabel("Freq (Hz)")
plt.ylabel("FFT Amplitude |X(freq)|**2")
plt.xlim(-0.1, 11)

plt.subplot(122)
plt.plot(t, fft.ifft(X), "r")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.tight_layout()
plt.show()


# %%% Calculation of the power
# Parserval's theorem
V = np.sum(x**2) / sr
V
P = np.sum(psd)
P

plt.figure(figsize=(12, 6))
plt.subplot(121)

plt.stem(freq, psd / P, "b", markerfmt=" ", basefmt="-b")
plt.xlabel("Freq (Hz)")
plt.ylabel("FFT Amplitude |X(freq)|**2")
plt.xlim(-0.1, 10)

plt.subplot(122)
plt.plot(t, fft.ifft(X), "r")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.tight_layout()
plt.show()

# normalized PSD
sum(psd) / P

# %% References
# https://appliedacousticschalmers.github.io/scaling-of-the-dft/AES2020_eBrief/

# %% Squared signal

import numpy as np
import matplotlib.pyplot as plt

# Parameters
frequency = 1  # Frequency of the signal (in Hz)
amplitude = 2  # Amplitude of the signal
duration = 5  # Duration of the signal (in seconds)
sampling_rate = 400  # Sampling rate (in Hz)

# Time array
t = np.arange(0, duration, 1 / sampling_rate)

# Generating the squared periodic signal
signal = amplitude * np.sign(np.sin(2 * np.pi * frequency * t))
signal += 1 * (np.sin(2 * np.pi * 2 * t))
signal += np.random.normal(1, 0.2, duration * sampling_rate)

# Plotting the signal
plt.plot(t, signal)
plt.title("Squared Periodic Signal")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.grid(True)
plt.show()

# %%% FFT
from scipy import fft

X = fft.fft(signal, norm="backward")
N = len(signal)
freq = fft.fftfreq(N, 1 / sampling_rate)

psd = np.abs(X) ** 2 / (sampling_rate * N)
# aquii esta la clave de la normalización. Por ello estudiar esto. Tambien usar la fft completa y no rfft para evitar problemas con parserval

# plt.figure(figsize=(12, 6))
plt.subplot(121)

plt.stem(freq, psd, "b", markerfmt=" ", basefmt="-b")
plt.xlabel("Freq (Hz)")
plt.ylabel("FFT Amplitude |X(freq)|**2")
plt.xlim(-0.1, 11)

plt.subplot(122)
plt.plot(t, fft.ifft(X), "r")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.tight_layout()
plt.show()

# %% Pulsed signal

# Parameters
frequency = 2.5  # Frequency of the pulses (in Hz)
pulse_width = 0.3  # Width of each pulse (in seconds)
amplitude = 2  # Amplitude of the pulses
duration = 5  # Duration of the signal (in seconds)
sampling_rate = 1000  # Sampling rate (in Hz)

# Time array
t = np.arange(0, duration, 1 / sampling_rate)

# Generating the pulsed signal
signal = amplitude * np.where(np.sin(2 * np.pi * frequency * t) > 0, 1, 0)
signal = np.where((t % (1 / frequency)) < pulse_width, signal, 0)
signal = signal.astype("float64") + 1 * (np.sin(2 * np.pi * 1.3 * t))
signal += np.random.normal(1, 0.2, duration * sampling_rate)
signal = signal - np.mean(signal)


# Plotting the signal
plt.plot(t, signal)
plt.title("Pulsed Signal")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.grid(True)
plt.show()

"# %%% FFT"

X = fft.fft(signal, norm="backward")
N = len(signal)
freq = fft.fftfreq(N, 1 / sampling_rate)

psd = np.abs(X) ** 2 / (sampling_rate * N)
# aquii esta la clave de la normalización. Por ello estudiar esto. Tambien usar la fft completa y no rfft para evitar problemas con parserval

# plt.figure(figsize=(12, 6))
plt.subplot(121)

plt.stem(freq, psd, "b", markerfmt=" ", basefmt="-b")
plt.xlabel("Freq (Hz)")
plt.ylabel("FFT Amplitude |X(freq)|**2")
plt.xlim(-0.1, 6)

plt.subplot(122)
plt.plot(t, fft.ifft(X), "r")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.tight_layout()
plt.show()

# %% same but convolved

# Parameters
frequency = 7  # Frequency of the pulses (in Hz)
pulse_width = 0.4  # Width of each pulse (in seconds)
amplitude = 2  # Amplitude of the pulses
duration = 5  # Duration of the signal (in seconds)
sampling_rate = 1000  # Sampling rate (in Hz)

# Time array
t = np.arange(0, duration, 1 / sampling_rate)

# Generating the pulsed signal
signal = amplitude * np.where(np.sin(2 * np.pi * frequency * t) > 0, 1, 0)
signal = np.where((t % (1 / frequency)) < pulse_width, signal, 0)
signal2 = 0.2 * (np.sin(2 * np.pi * 29 * t))
signal = np.convolve(signal2, signal, mode="same")

signal += np.random.normal(1, 0.2, duration * sampling_rate)
signal = signal - np.mean(signal)


# Plotting the signal
plt.plot(t, signal)
plt.title("Pulsed Signal")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.grid(True)
plt.show()


"# %%% FFT"

X = fft.fft(signal, norm="backward")
N = len(signal)
freq = fft.fftfreq(N, 1 / sampling_rate)

psd = np.abs(X) ** 2 / (sampling_rate * N)
# aquii esta la clave de la normalización. Por ello estudiar esto. Tambien usar la fft completa y no rfft para evitar problemas con parserval

# plt.figure(figsize=(12, 6))
plt.subplot(121)

plt.stem(freq, psd, "b", markerfmt=" ", basefmt="-b")
plt.xlabel("Freq (Hz)")
plt.ylabel("FFT Amplitude |X(freq)|**2")
plt.xlim(-0.1, 50)

plt.subplot(122)
plt.plot(t, fft.ifft(X), "r")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.tight_layout()
plt.show()
