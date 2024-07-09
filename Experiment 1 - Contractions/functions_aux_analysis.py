#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 19 19:15:39 2022

@author: Usuario

Funciones para el analisis de gusanos


"""

# from scipy.stats import sem

import os
import glob
import numpy as np


def get_files_in_folder(folder_path):
    file_list = glob.glob(folder_path + "/**", recursive=True)
    files = [file for file in file_list if not os.path.isdir(file)]
    return files


def agrupamiento_gusanos_fft(df, condicion):
    n_gus = sum([condicion in g for g in df.columns])
    df.insert(len(df.columns), "Suma", df.iloc[:, :n_gus].sum(axis=1))
    df.insert(len(df.columns), "Media", df.iloc[:, :n_gus].mean(axis=1))
    df.insert(len(df.columns), "Max", df.iloc[:, :n_gus].max(axis=1))
    df.insert(len(df.columns), "STD", df.iloc[:, :n_gus].std(axis=1))
    df.insert(len(df.columns), "SEM", df.iloc[:, :n_gus].sem(axis=1))


def remove_outliers_iqr(data):
    Q1 = np.percentile(data, 25)
    Q3 = np.percentile(data, 75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR
    return data[(data >= lower_bound) & (data <= upper_bound)]


# # change names of Batch 8

# folder_path = (
#       "/home/ale/pCloudDrive/CABD/Lab Ozren/Marta Fernandez/Experimento Coletazos/Batch 6 Results/"
# )
# file_list = os.listdir(folder_path)

# file_list
# for f in file_list:
#     os.rename(folder_path+f, folder_path + '150923_' + f[8:])
