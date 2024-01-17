# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 14:36:10 2024

@author: Ale Campoy
"""

import numpy as np

from scipy import stats

stats.ttest_ind(X1, X2).confidence_interval(confidence_level=0.95)
