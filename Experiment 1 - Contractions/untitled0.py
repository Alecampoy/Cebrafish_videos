# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 14:24:01 2024

@author: Ale Campoy
"""

import numpy as np
from scipy import stats
from pingouin import ttest
import statsmodels.stats.api as sms


def confidence_interval_diff_means(sample1, sample2, confidence_level=0.95):
    """
    Calculate the confidence interval for the difference of means between two samples.

    Parameters:
    - sample1: First sample data
    - sample2: Second sample data
    - confidence_level: Confidence level for the interval (default is 0.95)

    Returns:
    - Tuple containing the lower and upper bounds of the confidence interval
    """

    # Perform a two-sample t-test
    t_stat, p_value = stats.ttest_ind(sample1, sample2)

    # Degrees of freedom
    df = len(sample1) + len(sample2) - 2

    # Calculate the standard error of the difference of means
    std_error_diff = np.sqrt(
        (np.var(sample1) / len(sample1)) + (np.var(sample2) / len(sample2))
    )

    # Calculate the margin of error
    margin_of_error = (
        stats.t.ppf(confidence_level + (1 - confidence_level) / 2, df) * std_error_diff
    )

    # Calculate the confidence interval
    lower_bound = np.mean(sample1) - np.mean(sample2) - margin_of_error
    upper_bound = np.mean(sample1) - np.mean(sample2) + margin_of_error

    return lower_bound, upper_bound


# Example usage:
# sample1 = np.random.normal(loc=30, scale=10, size=100)
# sample2 = np.random.normal(loc=25, scale=12, size=120)

# confidence_interval = confidence_interval_diff_means(sample1, sample2)


def mean_diff_CI(X1=np.nan, X2=np.nan):
    if X1 is np.nan or X2 is np.nan:
        return np.nan
    else:
        cm = sms.CompareMeans(sms.DescrStatsW(X1), sms.DescrStatsW(X2))
        return cm.tconfint_diff(alpha=0.05, usevar="unequal")


# %%

X1 = [2, 3, 1, 23, 3, 4, 5, 2]
X2 = [4, 5, 6, 7, 8, 8, 7, 6, 7, 7]

confidence_interval_diff_means(X1, X2, 0.95)

mean_diff_CI(X1, X2)

stats.ttest_ind(X1, X2).confidence_interval(confidence_level=0.90)

ttest(X1, X2, confidence=0.90, correction=False)["CI90%"]
