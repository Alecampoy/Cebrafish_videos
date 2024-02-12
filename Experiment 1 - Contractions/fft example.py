# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 19:03:59 2024

@author: Ale Campoy
"""
# https://pythonnumericalmethods.berkeley.edu/notebooks/chapter24.04-FFT-in-Python.html
# https://pythonnumericalmethods.berkeley.edu/notebooks/chapter24.04-FFT-in-Python.html
import matplotlib.pyplot as plt
import numpy as np


# %%
# sampling rate
sr = 2000
# sampling interval
ts = 1.0 / sr
t = np.arange(0, 1, ts)

freq = 2.0
x = 1.5 * np.sin(2 * np.pi * freq * t)

freq = 4
x += 3 * np.sin(2 * np.pi * freq * t)

freq = 7
x += 2 * np.sin(2 * np.pi * freq * t)

x += 2
plt.figure(figsize=(8, 6))
plt.plot(t, x, "r")
plt.ylabel("Amplitude")

plt.show()

# %%
from scipy import fft

X = fft.fft(x, norm="backward")
N = len(X)
n = np.arange(N)
T = N / sr
freq = fft.fftfreq(N, ts)

psd = np.abs(X) ** 2 / (sr * N)

plt.figure(figsize=(12, 6))
plt.subplot(121)

plt.stem(freq, psd, "b", markerfmt=" ", basefmt="-b")
plt.xlabel("Freq (Hz)")
plt.ylabel("FFT Amplitude |X(freq)|**2")
plt.xlim(-10, 10)

plt.subplot(122)
plt.plot(t, fft.ifft(X), "r")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.tight_layout()
plt.show()

freqs = fft.fftfreq(N, ts)

# %% Calculation of the power
# Parserval's theorem
V = np.sum(x**2) / sr

P = np.sum(psd)  # [0:int(sr/2)]


plt.figure(figsize=(12, 6))
plt.subplot(121)

plt.stem(freq, psd / P, "b", markerfmt=" ", basefmt="-b")
plt.xlabel("Freq (Hz)")
plt.ylabel("FFT Amplitude |X(freq)|**2")
plt.xlim(0, 10)

plt.subplot(122)
plt.plot(t, fft.ifft(X), "r")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.tight_layout()
plt.show()

# normalized PSD

np.sum(psd / P)

# %% Welsh
from scipy.signal import welch

frequencies_welch, psd_welch = welch(X, sr, nperseg=512)

plt.figure(figsize=(12, 6))
plt.subplot(121)

plt.stem(freq, psd / P, "b", markerfmt=" ", basefmt="-b")
plt.xlabel("Freq (Hz)")
plt.ylabel("FFT Amplitude |X(freq)|**2")
plt.xlim(0, 10)

plt.subplot(122)
plt.plot(frequencies_welch, psd_welch, "r")
plt.xlim(0, 10)
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.tight_layout()
plt.show()
