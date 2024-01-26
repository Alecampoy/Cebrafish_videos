# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 19:03:59 2024

@author: Ale Campoy
"""
# https://pythonnumericalmethods.berkeley.edu/notebooks/chapter24.04-FFT-in-Python.html
# https://pythonnumericalmethods.berkeley.edu/notebooks/chapter24.04-FFT-in-Python.html
import matplotlib.pyplot as plt
import numpy as np

plt.style.use("seaborn-poster")


# %%
# sampling rate
sr = 2000
# sampling interval
ts = 1.0 / sr
t = np.arange(0, 1, ts)

freq = 1.0
x = 3 * np.sin(2 * np.pi * freq * t)

freq = 4
x += np.sin(2 * np.pi * freq * t)

freq = 7
x += 0.5 * np.sin(2 * np.pi * freq * t)

plt.figure(figsize=(8, 6))
plt.plot(t, x, "r")
plt.ylabel("Amplitude")

plt.show()

# %%
from scipy.fftpack import fft, ifft

X = fft(x)

plt.figure(figsize=(12, 6))
plt.subplot(121)

plt.stem(freq, np.abs(X), "b", markerfmt=" ", basefmt="-b")
plt.xlabel("Freq (Hz)")
plt.ylabel("FFT Amplitude |X(freq)|")
plt.xlim(0, 10)

plt.subplot(122)
plt.plot(t, ifft(X), "r")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.tight_layout()
plt.show()
