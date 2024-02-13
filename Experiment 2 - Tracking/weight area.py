# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 12:27:09 2024

@author: Ale Campoy
"""
import numpy as np
from random import randrange, uniform

R = 1

r = R * np.sqrt(np.random.uniform(0, 1, 1000))
weights = np.arange(0, 1 + 1 / 10, 1 / 10)

bins = np.arange(0, 1 + 1 / 12, 1 / 12)
weights = np.pi * (bins**2)

weights = np.diff(weights)  # [:len(weights)-1]

weights = np.repeat(weights, np.diff(np.searchsorted(bins, r)))
a = np.diff(np.searchsorted(bins, r))
# weights = 1/(r)
import matplotlib.pyplot as plt

count, bins, ignored = plt.hist(r, 12, weights=weights, density=True)

plt.plot(bins, np.ones_like(bins), linewidth=2, color="r")

plt.show()
