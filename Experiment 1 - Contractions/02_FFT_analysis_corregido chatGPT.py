# -*- coding: utf-8 -*-
# Spyder Editor

# %% Intro [md]
"""
# **Análisis de Movimiento Zebrafish**
### Autor: Alejandro Campoy López  
"""

# %% Librerías
import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.signal import find_peaks, find_peaks_cwt, detrend, periodogram, lombscargle
from scipy.fft import fft, rfft, fftfreq, rfftfreq
from scipy.linalg import dft
from functions_aux_analysis import *

plt.rcParams["figure.figsize"] = (15, 8)

# %% Lectura Archivos [md]
'''
Lectura de todos los archivos csv con los resultados de los diferentes batches.
Se añade una columna representando el gusano y el batch mediante el uso de regex
'''

# %%% Load Files

windows = False  # Comentario: Variable no utilizada
if windows:
    folder_path = "p:\\CABD\\Lab Ozren\\Marta Fernandez\\Experimento Coletazos\\"
else:
    folder_path = (
        "/home/ale/pCloudDrive/CABD/Lab Ozren/Marta Fernandez/Experimento Coletazos/"
    )
files = get_files_in_folder(folder_path)

df = []
for f in files:
    if ".csv" in f:
        csv = pd.read_csv(f, sep=";").drop(
            [
                "FeretAngle",
                "NBranches",
                "AvgBranchLen",
                "MaxBranchLen",
                "BranchLen",
                "EuclideanDist",
                "AR",  # this is the inverse of Roundness
            ],
            axis=1,
        )
        csv.insert(0, "Batch", re.search("Batch \d+", f).group(0))
        csv.insert(1, "Phenotype", re.search("KO\d*|WT", f.upper()).group(0))
        csv.insert(2, "Fish", "ZebraF_" + re.search("(\d+)(.lif)", f).group(1))
        df.append(csv)
        del (csv, f)

df = pd.concat(df)
# renombro la columna con nombre repetido
df = df.rename(columns={"XM.1": "YM"})

# renombro KO a KO44 para el batch 6 y 7
df.loc[df.Phenotype == "KO", "Phenotype"] = "KO44"

# %% NAs [md]
'''
## Número de NAs por ZebraF
Visualizamos el número de Nas por pez. Se imputan mediante interpolación Lineal.
'''

# %%% NAs Plot
NAs = (
    df.groupby(["Batch", "Phenotype", "Fish"])
    .apply(lambda x: x.isnull().sum())[["area"]]
    .rename(columns={"area": "NAs"})
).reset_index()  # me quedo solo con una columna ya que el numero de NAN es el mismo en todas

NAs["Batch_Phenotype"] = NAs.Batch + "_" + NAs.Phenotype
# NAs_barplot = sns.barplot(x="Fish", y="NAs", hue="Batch", data=NAs.reset_index())
# plt.xticks(rotation=90)
# plt.show()

NAs_barplot = sns.catplot(
    kind="bar", data=NAs.reset_index(), x="Fish", y="NAs", col="Phenotype", row="Batch"
)
NAs_barplot.set_xticklabels(rotation=90)
plt.show()

# %%% NA Impute

df = df.interpolate(method="linear")

# %%% [md]
'''El Zebra 10 WT del batch 7 se ha eliminado por contener > 500 NAs'''

# %% Distancia Recorrida [md]
'''
## Distancia Recorrida
Se calcula la distancia que recorre el pez a lo largo del video y se grafica por batch
'''

# %%% Calculo de la distancia

df.insert(7, "X_diff", df.groupby(["Batch", "Phenotype", "Fish"]).XM.diff())
df.insert(8, "Y_diff", df.groupby(["Batch", "Phenotype", "Fish"]).YM.diff())
df.insert(9, "dist", np.sqrt((df.X_diff**2) + (df.Y_diff**2)))

# dataframe con la distancia recorrida por el  gusano
Dist = df.groupby(["Batch", "Phenotype", "Fish"])[["dist"]].sum().round().reset_index()


# %%% Box-plot por batch

grped_bplot = sns.catplot(
    x="Batch",
    y="dist",
    kind="box",
    data=Dist,
    col="Phenotype",
    showfliers=False,
    notch=True,
)
grped_bplot.set_xticklabels(rotation=45)
plt.show()

# %%%  Barplot

Dist_bplot = sns.catplot(
    x="Phenotype",
    y="dist",
    kind="bar",
    data=Dist,
    col="Batch",
    capsize=0.2,
    errwidth=2.5,
    dodge=True,
)
Dist_bplot.set_xticklabels(rotation=45)
plt.show()

# %% Media Velocidades por minuto [md]
'''
## Media Velocidades por minuto
Se calcula la media de velocidades por minuto para cada batch
'''

# %%% Add seconds column
df.insert(4, "secs", (df.MM * 60) + df.SS)

# %%% Speed mean per min
Speed_per_min = df.groupby(["Batch", "Phenotype", "Fish", "MM"])[["Speed"]].mean()
Speed_per_min.reset_index(inplace=True)
Speed_per_min = Speed_per_min.rename(columns={"Speed": "Speed_mean_min"})

# %%% Speed boxplot
grped_bplot = sns.catplot(
    x="Batch",
    y="Speed_mean_min",
    kind="box",
    data=Speed_per_min,
    col="Phenotype",
    showfliers=False,
    notch=True,
)
grped_bplot.set_xticklabels(rotation=45)
plt.show()

# %% Velocidades Medias por ZebraF y minuto [md]
'''
## Velocidades Medias por ZebraF y minuto
Se calcula la media de velocidades para cada pez a lo largo de todos los minutos y se grafica por batch
'''

# %%% Fish speeds

Fish_speeds = (
    df.groupby(["Batch", "Phenotype", "Fish"])[["Speed"]].mean().round().reset_index()
)

# %%% Fish speeds boxplot
grped_bplot = sns.catplot(
    x="Batch",
    y="Speed",
    kind="box",
    data=Fish_speeds,
    col="Phenotype",
    showfliers=False,
    notch=True,
)
grped_bplot.set_xticklabels(rotation=45)
plt.show()

# %%% Velocities per fish
