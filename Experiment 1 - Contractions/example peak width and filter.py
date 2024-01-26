# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 11:51:04 2023

@author: Ale Campoy
"""


import numpy as np

import pandas as pd

import matplotlib.pyplot as plt

from scipy.datasets import electrocardiogram

from scipy.signal import find_peaks

x = pd.read_csv("k.txt")["345"]
peaks = find_peaks(x, prominence=1, width=40)[1]
peaks, properties = find_peaks(x, prominence=1, width=40)

properties["prominences"], properties["widths"]

plt.plot(x)

plt.plot(peaks, x[peaks], "x")

plt.vlines(
    x=peaks, ymin=x[peaks] - properties["prominences"], ymax=x[peaks], color="C1"
)

plt.hlines(
    y=properties["width_heights"],
    xmin=properties["left_ips"],
    xmax=properties["right_ips"],
    color="C1",
)

plt.show()
