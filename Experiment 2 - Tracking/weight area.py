# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 12:27:09 2024

@author: Ale Campoy
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

R = 1
r = 1 - R * np.sqrt(
    np.random.uniform(0, 1, 6666)
)  # distribución homogenea sobre una circunferencia
## Si la distribución va de 1 a 0 por un cambio de variable, hay que especificar 1-r en los pesos
nbins = 10
bins = np.arange(0, 1 + (1 / (nbins)), 1 / (nbins))

# metodo 1 para los pesos
weights_a = (np.pi * np.arange(0, 1 + (1 / (nbins)), 1 / (nbins)) ** 2)[::-1]
weights_a = np.diff(weights_a)  # Diferencias del area
bin_of_r = np.searchsorted(bins, r) - 1  # bin al que pertenece cada observación
weights_a_ind = 1 / weights_a[bin_of_r]
# b = np.diff(np.searchsorted(bins, r))

# método 2 para los pesos
weights_r = 1 / (np.pi * (1 - r))  # dada la naturaleza radial de los datos

sns.histplot(r, stat="density", binrange=[0, 1], bins=10)

plt.show()
sns.histplot(x=r, stat="density", weights=weights_r, binrange=[0, 1], bins=10)


count, bins = np.histogram(r, bins, density=True)
plt.hist(r, bins, weights=weights_a_ind, density=True, alpha=0.5, color="g")

# plt.hist(r, bins, weights=weights_r, density=True, alpha=0.5, color="b")
plt.plot(bins, np.ones_like(bins), linewidth=2, color="r")  # uniform line
# plt.plot(
#     bins[: len(bins) - 1] + 0.05, count / (np.pi * weights_a), linewidth=2, color="g"
# )
plt.show()
