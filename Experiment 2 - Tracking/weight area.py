# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 12:27:09 2024

@author: Ale Campoy
"""
import numpy as np
import matplotlib.pyplot as plt

R = 1
r = R * np.sqrt(
    np.random.uniform(0, 1, 5000)
)  # distribución homogenea sobre una circunferencia

nbins = 10
bins = np.arange(0, 1 + (1 / (nbins)), 1 / (nbins))
weights_a = np.pi * np.arange(0, 1 + (1 / (nbins)), 1 / (nbins)) ** 2
weights_a = np.diff(weights_a)  # [:len(weights)-1]

bin_of_r = np.searchsorted(bins, r) - 1  #
weights_a_ind = 1 / weights_a[bin_of_r]
# b = np.diff(np.searchsorted(bins, r))
weights_r = 1 / (np.pi * r)  # dada la naturaleza radial de los datos

sns.histplot(r, stat="density", binrange=[0, 1], bins=10)

plt.show()

# %%


count, bins = np.histogram(r, bins, density=True)
plt.hist(r, bins, weights=weights_r, density=True, alpha=0.5, color="b")
plt.hist(r, bins, weights=weights_a_ind, density=True, alpha=0.5, color="g")

plt.plot(bins, np.ones_like(bins), linewidth=2, color="r")
plt.plot(
    bins[: len(bins) - 1] + 0.05, count / (np.pi * weights_a), linewidth=2, color="g"
)

plt.show()
