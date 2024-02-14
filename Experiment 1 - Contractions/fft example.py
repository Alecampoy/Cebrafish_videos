# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 19:03:59 2024

@author: Ale Campoy
"""
# https://pythonnumericalmethods.berkeley.edu/notebooks/chapter24.04-FFT-in-Python.html
# https://pythonnumericalmethods.berkeley.edu/notebooks/chapter24.04-FFT-in-Python.html
import matplotlib.pyplot as plt
import numpy as np

# sampling rate
sr = 2500
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
x += np.random.normal(1, 1.3, sr)

plt.figure(figsize=(8, 6))
plt.plot(t, x, "r")
plt.ylabel("Amplitude")

plt.show()

# %%
from scipy import fft

X = fft.rfft(x, norm="backward")
N = len(x)
freq = fft.rfftfreq(N, ts)

psd = np.abs(X) ** 2 / (
    sr * N
)  # aquii esta la clave de la normalización. He encontrado que con diferentes sr se obtienen diferntes valores. Por ello estudiar esto. Tambien usar la fft completa y no rfft para evitar problemas con parserval

plt.figure(figsize=(12, 6))
plt.subplot(121)

plt.stem(freq, psd, "b", markerfmt=" ", basefmt="-b")
plt.xlabel("Freq (Hz)")
plt.ylabel("FFT Amplitude |X(freq)|**2")
plt.xlim(0.1, 11)

plt.subplot(122)
plt.plot(t, fft.irfft(X), "r")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.tight_layout()
plt.show()


# %% Calculation of the power
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
plt.xlim(0, 10)

plt.subplot(122)
plt.plot(t, fft.irfft(X), "r")
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
