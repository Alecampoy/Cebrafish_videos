#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 19 19:15:39 2022

@author: Usuario

Funciones para el analisis de gusanos


"""

# from scipy.stats import sem


def agrupamiento_gusanos_fft(df, condicion):
    n_gus = sum([condicion in g for g in df.columns])
    df.insert(len(df.columns), "Suma", df.iloc[:, :n_gus].sum(axis=1))
    df.insert(len(df.columns), "Media", df.iloc[:, :n_gus].mean(axis=1))
    df.insert(len(df.columns), "Max", df.iloc[:, :n_gus].max(axis=1))
    df.insert(len(df.columns), "STD", df.iloc[:, :n_gus].std(axis=1))
    df.insert(len(df.columns), "SEM", df.iloc[:, :n_gus].sem(axis=1))