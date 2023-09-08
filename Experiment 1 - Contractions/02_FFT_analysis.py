# -*- coding: utf-8 -*-
# """
# Spyder Editor

# This is a temporary script file.
# """

# %% librerias

import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.signal import find_peaks, detrend, periodogram, lombscargle
from scipy.fft import fft, rfft, fftfreq, rfftfreq
from scipy.linalg import dft
from functions_aux_analysis import *

# %% Lectura de todos los archivos csv con los resultados de los diferentes batches.
# Se añade una columna representando el gusano y el batch mediante el uso de regex

windows = False
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
        csv = pd.read_csv(f, sep=";").drop(["FeretAngle"], axis=1)
        csv.insert(0, "Batch", re.search("Batch \d*", f).group(0))
        csv.insert(1, "Fish", "ZebraF " + re.search("(?<=hpf_)(\d*)(?=.)", f).group(0))
        df.append(csv)
        del (csv, f)

df = pd.concat(df)
# renombro la columna con nombre repetido
df = df.rename(columns={"XM.1": "YM"})


# %% Calculo de las magnitudes derivadas y añado la condición al dataset

df["Curvatura"] = df.EuclideanDist / df.BranchLen

fenotype = pd.read_excel(folder_path+"Fenotype.ods")
# añado fenotipo
df = pd.merge(df, fenotype, on=["Batch", "Fish"])


# %% Número de NAs que hay por frame. Imputo mediante interpolación


NAs = (
    df.groupby(["Batch", "Fish"])
    .apply(lambda x: x.isnull().sum())[["area"]]
    .rename(columns={"area": "NAs"})
)  # me quedo solo con una columna ya que el numero de NAN es el mismo en todas

NAs_barplot = sns.barplot(x="Fish", y="NAs", hue="Batch", data=NAs.reset_index())
plt.xticks(rotation=90)
plt.show()

NAs_barplot = sns.catplot(
    kind="bar", data=NAs.reset_index(), x="Fish", y="NAs", col="Batch"
)
NAs_barplot.set_xticklabels(rotation=90)
plt.show()

df = df.interpolate(method="linear")


# %% Distancia Recorrida

# distancia en cada paso

df.insert(6, "X_diff", df.groupby(["Batch", "Fish"]).XM.diff())
df.insert(7, "Y_diff", df.groupby(["Batch", "Fish"]).YM.diff())
df.insert(8, "dist", np.sqrt((df.X_diff**2) + (df.Y_diff**2)))

# dataframe con la distancia recorrida por el  gusano
Dist = df.groupby(["Batch", "Fish"])[["dist"]].sum().round().reset_index()
Dist = pd.merge(Dist, fenotype, on=["Batch", "Fish"], how="left").dropna()

# %%% Box-Plot de la distancia
a = sns.boxplot(x="Fenotype", y="dist", data=Dist)
a.set_title("Distancia Recorrida por el gusano en pixeles")
b = sns.stripplot(x="Fenotype", y="dist", data=Dist, color="grey", size=8)
plt.show()

# %%% Box-plot por batch

grped_bplot = sns.catplot(
    x="Batch",
    y="dist",
    hue="Fenotype",
    kind="box",
    legend=False,
    height=6,
    aspect=1.9,
    data=Dist,
)
# make grouped stripplot
grped_bplot = sns.stripplot(
    x="Batch",
    y="dist",
    hue="Fenotype",
    jitter=True,
    dodge=True,
    marker="o",
    color="black",
    # palette="Set2",
    data=Dist,
)
handles, labels = grped_bplot.get_legend_handles_labels()

l = plt.legend(handles[0:3], labels[0:3])
plt.show()


# %% Repliegamientos


# #%%% Mediante threshold (No apropiada)
# # esta medida no es apropiada, ya que indica el tiempo que pasa con algo de replegamiento y no el numero de veces que se repliega.
# Para el tiempo que pasa replegado, mejor usar solidity
# #%%%% Roundness
# # Para contar los repliegamientos, defino un threshold y simplemente cuento

# threshold = 0.5
# replieg_round = (
#     df[df.Round > threshold]
#     .groupby("Gusano")
#     .count()[["area"]]
#     .rename(columns={"area": "Repliegamientos"})
#     .reset_index()
# )
# replieg_round.insert(
#     0,
#     "Condicion",
#     replieg_round.Gusano.apply(
#         lambda x: "WT" if x[0 : x.index(" ")] == "CONTROL" else "MUT"
#     ),
# )

# a = sns.boxplot(x="Condicion", y="Repliegamientos", data=replieg_round)
# a.set_title("Numero de repliegamientos con Roundness y Threshold " + str(threshold))
# b = sns.stripplot(
#     x="Condicion", y="Repliegamientos", data=replieg_round, color="grey", size=8
# )
# plt.show()

# #%%%% Circularity
# # Para contar los repliegamientos, defino un threshold y simplemente cuento

# threshold = 0.7
# replieg_circ = (
#     df[df.Circ > threshold]
#     .groupby("Gusano")
#     .count()[["area"]]
#     .rename(columns={"area": "Repliegamientos"})
#     .reset_index()
# )
# replieg_circ.insert(
#     0,
#     "Condicion",
#     replieg_circ.Gusano.apply(
#         lambda x: "WT" if x[0 : x.index(" ")] == "CONTROL" else "MUT"
#     ),
# )

# a = sns.boxplot(x="Condicion", y="Repliegamientos", data=replieg_circ)
# a.set_title("Numero de repliegamientos con circularity y Threshold " + str(threshold))
# b = sns.stripplot(
#     x="Condicion", y="Repliegamientos", data=replieg_circ, color="grey", size=8
# )
# plt.show()

# %%% Mediante peak finder
# %%%% Función peak finder y plot de los peaks detectados
#
# Para plotear los peaks encontrados usar :
g_temp = df.Round[df.Gusano == "Series007"].to_numpy()  # reset_index(drop=True)
peaks, _ = find_peaks(g_temp, height=0.45, prominence=0.0, threshold=0.0, distance=5)
plt.plot(g_temp)
plt.plot(peaks, g_temp[peaks], "2")
plt.show()

# %%%% Repliegamientos dados por Roundness

# Creo un dataframe con el número de peaks de cada gusano
peaks_round = pd.DataFrame(columns=["Gusano", "N_peaks"])
for g in set(df.Gusano):
    g_temp = df.Round[df.Gusano == g].to_numpy()  # reset_index(drop=True)
    peaks, _ = find_peaks(
        g_temp, height=0.45, prominence=0.0, threshold=0.0, distance=5
    )
    peaks_round = pd.concat(
        [peaks_round, pd.DataFrame({"Gusano": g, "N_peaks": [len(peaks)]})],
    ).reset_index(drop=True)
peaks_round.insert(
    0,
    "Condicion",
    peaks_round.Gusano.apply(
        lambda x: "WT" if x[0 : x.index(" ")] == "CONTROL" else "MUT"
    ),
)

a = sns.boxplot(x="Condicion", y="N_peaks", data=peaks_round, order=["WT", "MUT"])
a.set_title("Numero de repliegamientos con Find_Peaks y Roundness")
b = sns.stripplot(
    x="Condicion",
    y="N_peaks",
    data=peaks_round,
    color="grey",
    size=8,
    order=["WT", "MUT"],
)
plt.show()


# #%%%% Circularity
# # La circularity NO es adecuada ya que solo se vuelve pequeña cuando el gusano se toca

# # Creo un dataframe con el número de peaks de cada gusano
# peaks_circ = pd.DataFrame(columns=["Gusano", "N_peaks"])
# for g in set(df.Gusano):
#     g_temp = df.Circ[df.Gusano == g].to_numpy()  # reset_index(drop=True)
#     peaks, _ = find_peaks(
#         g_temp, height=0.45, prominence=0.0, threshold=0.0, distance=5
#     )
#     peaks_circ = pd.concat(
#         [peaks_circ, pd.DataFrame({"Gusano": g, "N_peaks": [len(peaks)]})],
#     ).reset_index(drop=True)
# peaks_circ.insert(
#     0,
#     "Condicion",
#     peaks_circ.Gusano.apply(
#         lambda x: "WT" if x[0 : x.index(" ")] == "CONTROL" else "MUT"
#     ),
# )

# a = sns.boxplot(x="Condicion", y="N_peaks", data=peaks_circ, order=["WT", "MUT"])
# a.set_title("Numero de repliegamientos con Find_Peaks y Circularity")
# b = sns.stripplot(
#     x="Condicion",
#     y="N_peaks",
#     data=peaks_circ,
#     color="grey",
#     size=8,
#     order=["WT", "MUT"],
# )
# plt.show()

# %% Euclidean distance # Maximos Mediante peak finder

# Creo un dataframe con el número de peaks de cada gusano
peaks_euclidean = pd.DataFrame(columns=["Gusano", "N_peaks"])
for g in set(df.Gusano):
    g_temp = df.EuclideanDist[df.Gusano == g].to_numpy()  # reset_index(drop=True)
    peaks, _ = find_peaks(
        g_temp, height=0.45, prominence=0.0, threshold=0.0, distance=5
    )
    peaks_euclidean = pd.concat(
        [peaks_euclidean, pd.DataFrame({"Gusano": g, "N_peaks": [len(peaks)]})],
    ).reset_index(drop=True)
peaks_euclidean.insert(
    0,
    "Condicion",
    peaks_euclidean.Gusano.apply(
        lambda x: "WT" if x[0 : x.index(" ")] == "CONTROL" else "MUT"
    ),
)

a = sns.boxplot(x="Condicion", y="N_peaks", data=peaks_euclidean, order=["WT", "MUT"])
a.set_title("Numero de peaks en euclidean distance")
b = sns.stripplot(
    x="Condicion",
    y="N_peaks",
    data=peaks_euclidean,
    color="grey",
    size=8,
    order=["WT", "MUT"],
)
plt.show()


# %% Solidity
# Espero con esta medida comprobar que el WT esta mas tiempo en una posición contraida

# %%% media

solidity_mean = df.groupby("Gusano").mean().reset_index()
solidity_mean.insert(
    0,
    "Condicion",
    solidity_mean["Gusano"].apply(
        lambda x: "WT" if x[0 : x.index(" ")] == "CONTROL" else "MUT"
    ),
)

a = sns.boxplot(x="Condicion", y="Solidity", data=solidity_mean, order=["WT", "MUT"])
a.set_title("Solidity media")
b = sns.stripplot(
    x="Condicion",
    y="Solidity",
    data=solidity_mean,
    color="grey",
    size=8,
    order=["WT", "MUT"],
)
plt.show()

# %%% threshold

threshold = 0.7
solidity_thr = (
    df[df.Solidity > threshold]
    .groupby("Gusano")
    .count()[["Solidity"]]  # .rename(columns={"Solidity": "solidity_thr"})
    .reset_index()
)
solidity_thr.insert(
    0,
    "Condicion",
    solidity_thr.Gusano.apply(
        lambda x: "WT" if x[0 : x.index(" ")] == "CONTROL" else "MUT"
    ),
)

a = sns.boxplot(x="Condicion", y="Solidity", data=solidity_thr)
a.set_title("Tiempo total con una solidity mayor que " + str(threshold))
b = sns.stripplot(x="Condicion", y="Solidity", data=solidity_thr, color="grey", size=8)
plt.show()


# %% Curvatura
# %%% Mediante threshold

# Esta medida indica el numero de frames en los que el gusano pasa estirado
# Curvatura cercana a 1 indica que el gusano esta estirado
threshold = 0.9
curvatura_frames = (
    df[df.Curvatura > threshold].groupby("Gusano").count()[["Curvatura"]].reset_index()
)
curvatura_frames.insert(
    0,
    "Condicion",
    curvatura_frames.Gusano.apply(
        lambda x: "WT" if x[0 : x.index(" ")] == "CONTROL" else "MUT"
    ),
)

a = sns.boxplot(x="Condicion", y="Curvatura", data=curvatura_frames)
a.set_title("Numero de frames con curvatura > " + str(threshold))
b = sns.stripplot(
    x="Condicion", y="Curvatura", data=curvatura_frames, color="grey", size=8
)
plt.show()


# %%% Mediante peak finder

# Creo un dataframe con el número de peaks de cada gusano
peaks_curvatura = pd.DataFrame(columns=["Gusano", "N_peaks"])
for g in set(df.Gusano):
    g_temp = df.Curvatura[df.Gusano == g].to_numpy()  # reset_index(drop=True)
    peaks, _ = find_peaks(
        g_temp, height=0.45, prominence=0.0, threshold=0.0, distance=5
    )
    peaks_curvatura = pd.concat(
        [peaks_curvatura, pd.DataFrame({"Gusano": g, "N_peaks": [len(peaks)]})],
    ).reset_index(drop=True)
peaks_curvatura.insert(
    0,
    "Condicion",
    peaks_curvatura.Gusano.apply(
        lambda x: "WT" if x[0 : x.index(" ")] == "CONTROL" else "MUT"
    ),
)

a = sns.boxplot(x="Condicion", y="N_peaks", data=peaks_curvatura, order=["WT", "MUT"])
a.set_title("Numero de peaks con Find_Peaks para la curvatura")
b = sns.stripplot(
    x="Condicion",
    y="N_peaks",
    data=peaks_curvatura,
    color="grey",
    size=8,
    order=["WT", "MUT"],
)
plt.show()

# %% Analisis de Frecuencias del movimiento para un gusano

# %%% Detrend
# Muestro como se normaliza mediante la función detrend
g = "CONTROL 1"
x = df[["Curvatura"]][df.Gusano == g]
t = df["T_seg"][df.Gusano == g]

plt.figure(figsize=(12, 6))
plt.subplot(121)

plt.plot(t, x.Curvatura)

plt.subplot(122)
plt.plot(t, detrend(x, axis=0), "r")
plt.title("Serie Temporal sin tendencia")
plt.tight_layout()
plt.show()

# %%% FFT
# la duracion esta ajustada para cada sample por haber eliminado los NA
sample_rate = 900 / 60
duration = np.int(len(x) / sample_rate)
N_points = len(x)

# usando la función FFT
g_fft = rfft(detrend(x, axis=0), axis=0, norm="forward")
# Calculo de las frecuencias pare representar en el eje X
x_fft = rfftfreq(N_points, 1 / sample_rate)


plt.plot(x_fft, np.abs(g_fft))
plt.title("FFT")
plt.xlabel("Frecuency (Hz)")
plt.show()


# %%% DFT
# def DFT(x):
#     """
#     Function to calculate the
#     discrete Fourier Transform
#     of a 1D real-valued signal x
#     """

#     N = len(x)
#     n = np.arange(N)
#     k = n.reshape((N, 1))
#     e = np.exp(-2j * np.pi * k * n / N)

#     X = np.dot(e, x)

#     return X


# #%%% Prueba DFT
# g_dft = DFT(detrend(x, axis=0))

# plt.plot(x_fft, np.abs(g_dft)[: int(len(g_dft) / 2) + 1])
# plt.xlabel("Frecuency (Hz)")
# plt.show()


# #%%% mismo con la funcion DFT
# dft_matrix = dft(len(x))

# g_dft = dft_matrix @ detrend(x, axis=0)

# plt.plot(x_fft, np.abs(g_dft)[: int(len(g_dft) / 2) + 1])
# plt.xlabel("Frecuency (Hz)")
# plt.show()


# %% Analisis de Frecuencias (FFT) por condición
# %%% Interpolando FFT
# Se consigue así que todos tengan la misma longitud
# %%%% ajuste de un gusano (ejemplo)

x_min = 0
x_max = max(x_fft)
n = len(x_fft)
x_new = np.linspace(0.00, 6, 100)  # parámetro que ajusta
y = np.abs(g_fft).reshape(-1)
g_fft_interp = np.interp(x_new, x_fft, y)

plt.plot(x_fft, y, linewidth=0.6, color="blue")
plt.plot(x_new, g_fft_interp, linewidth=0.8, color="orange")
plt.xlabel("Frecuency (Hz)")
plt.show()

# ejemplo escalado
# x = np.linspace(0, 2 * np.pi, 10)
# y = np.sin(x)
# xvals = np.linspace(0, 2 * np.pi, 50)
# g_fft_interp = np.interp(xvals, x, y)

# %%%% DF para todos los gusanos

sample_rate = 900 / 60
max_freq = 6
rate_N = 500
x_scaled = np.linspace(0.00, max_freq, rate_N)

# WT
fft_scaled_WT = pd.DataFrame()
fft_scaled_WT.insert(0, "Freq", x_scaled)
fft_scaled_WT = fft_scaled_WT.set_index("Freq")
for g in set(df.Gusano):
    if g[0 : g.index(" ")] == "CONTROL":
        x_temp = df[["Curvatura"]][df.Gusano == g]
        duration_temp = np.int(len(x_temp) / sample_rate)
        N_points = len(x_temp)
        g_fft = rfft(detrend(x_temp, axis=0), axis=0, norm="forward")
        x_fft = rfftfreq(N_points, 1 / sample_rate)
        g_fft_interp = np.interp(x_scaled, x_fft, np.abs(g_fft).reshape(-1), right=99)
        fft_scaled_WT.insert(len(fft_scaled_WT.columns), g, g_fft_interp)
# MUT
fft_scaled_MUT = pd.DataFrame()
fft_scaled_MUT.insert(0, "Freq", x_scaled)
fft_scaled_MUT = fft_scaled_MUT.set_index("Freq")
for g in set(df.Gusano):
    if g[0 : g.index(" ")] == "MUT":
        x_temp = df[["Curvatura"]][df.Gusano == g]
        duration_temp = np.int(len(x_temp) / sample_rate)
        N_points = len(x_temp)
        g_fft = rfft(detrend(x_temp, axis=0), axis=0, norm="forward")
        x_fft = rfftfreq(N_points, 1 / sample_rate)
        g_fft_interp = np.interp(x_scaled, x_fft, np.abs(g_fft).reshape(-1), right=99)
        fft_scaled_MUT.insert(len(fft_scaled_MUT.columns), g, g_fft_interp)
funciones_analisis.agrupamiento_gusanos_fft(fft_scaled_MUT, "MUT")
funciones_analisis.agrupamiento_gusanos_fft(fft_scaled_WT, "CONTROL")


# %%%% plot Medias

plt.plot(x_scaled, fft_scaled_WT.Media, linewidth=0.8, color="blue")
plt.fill_between(
    x_scaled,
    (fft_scaled_WT.Media - fft_scaled_WT.SEM),
    (fft_scaled_WT.Media + fft_scaled_WT.SEM),
    color="blue",
    linewidth=0.1,
    alpha=0.2,
)
plt.plot(x_scaled, fft_scaled_MUT.Media, linewidth=0.8, color="orange")
plt.fill_between(
    x_scaled,
    (fft_scaled_MUT.Media - fft_scaled_MUT.SEM),
    (fft_scaled_MUT.Media + fft_scaled_MUT.SEM),
    color="orange",
    linewidth=0.1,
    alpha=0.2,
)

plt.xlabel("Frecuency (Hz)")
plt.show()

# %%%% plot Suma

plt.plot(x_scaled, fft_scaled_WT.Suma, linewidth=0.8, color="blue")
plt.plot(x_scaled, fft_scaled_MUT.Suma, linewidth=0.8, color="orange")
plt.xlabel("Frecuency (Hz)")
plt.show()

# %%%% plot Max

plt.plot(x_scaled, fft_scaled_WT.Max, linewidth=0.8, color="blue")
plt.plot(x_scaled, fft_scaled_MUT.Max, linewidth=0.8, color="orange")
plt.xlabel("Frecuency (Hz)")
plt.show()

# %%% Ajustando la longitud de entrada a la FFT a un número comun de frames

# %%%% DF generado
limite_t = N_points = 300
x_fft = rfftfreq(N_points, 1 / sample_rate)

# WT
fft_WT = pd.DataFrame()
fft_WT.insert(0, "Freq", x_fft)
fft_WT = fft_WT.set_index("Freq")
for g in set(df.Gusano):
    if g[0 : g.index(" ")] == "CONTROL":
        x_temp = df[["Curvatura"]][df.Gusano == g][:limite_t]
        g_fft = rfft(detrend(x_temp, axis=0), axis=0, norm="forward")
        fft_WT.insert(len(fft_WT.columns), g, np.abs(g_fft))
# MUT
fft_MUT = pd.DataFrame()
fft_MUT.insert(0, "Freq", x_fft)
fft_MUT = fft_MUT.set_index("Freq")
for g in set(df.Gusano):
    if g[0 : g.index(" ")] == "MUT":
        x_temp = df[["Curvatura"]][df.Gusano == g][:limite_t]
        g_fft = rfft(detrend(x_temp, axis=0), axis=0, norm="forward")
        fft_MUT.insert(len(fft_MUT.columns), g, np.abs(g_fft))
funciones_analisis.agrupamiento_gusanos_fft(fft_WT, "CONTROL")
funciones_analisis.agrupamiento_gusanos_fft(fft_MUT, "MUT")

# %%%% plot Medias

plt.plot(x_fft, fft_WT.Media, linewidth=0.8, color="blue")
plt.fill_between(
    x_fft,
    (fft_WT.Media - fft_WT.SEM),
    (fft_WT.Media + fft_WT.SEM),
    color="blue",
    linewidth=0.1,
    alpha=0.2,
)
plt.plot(x_fft, fft_MUT.Media, linewidth=0.8, color="orange")
plt.fill_between(
    x_fft,
    (fft_MUT.Media - fft_MUT.SEM),
    (fft_MUT.Media + fft_MUT.SEM),
    color="orange",
    linewidth=0.1,
    alpha=0.2,
)

plt.xlabel("Frecuency (Hz)")
plt.show()

# %%%% plot Suma

plt.plot(x_fft, fft_WT.Suma, linewidth=0.8, color="blue")
plt.plot(x_fft, fft_MUT.Suma, linewidth=0.8, color="orange")
plt.xlabel("Frecuency (Hz)")
plt.show()

# %%%% plot Max

plt.plot(x_fft, fft_WT.Max, linewidth=0.8, color="blue")
plt.plot(x_fft, fft_MUT.Max, linewidth=0.8, color="orange")
plt.xlabel("Frecuency (Hz)")
plt.show()

# %% Periodograma gusano y 'Curvatura'
# %%% funcion periodogram (usa FFT)
# Creo que no interesa de este modo ya que tenemos la FFT calculada antes

ventanas = ["boxcar", "blackman", "hann", "flattop", "nuttall"]


for ventana in ventanas:
    x_period, y_period = periodogram(
        detrend(x, axis=0),
        fs=sample_rate,
        nfft=300,
        axis=0,
        window=ventana,
        scaling="spectrum",
    )

    plt.plot(x_period, y_period)
    plt.xlabel("Frecuency (Hz)")
    plt.title("Función periodograma (usa FFT) y una ventana " + ventana)
    plt.show()
# %%% Lomb-Scargle Periodogram
periods = np.linspace(0.00001, 10, 1000)
f_periodogram = 2 * np.pi / periods
y_periodogram = lombscargle(
    t, detrend(x.Curvatura, axis=0), f_periodogram, normalize=True
)

plt.plot(periods, y_periodogram)
plt.xlabel("Period (Seg)")
plt.title("Función lombscargle -> ajusta por LSE")
plt.show()

# hacer para todos los gusanos y sumar

# %% Periodogramas por condición
# %%% Funcion periodogram conteniendo todos los gusano

# ventanas = ["boxcar","blackman","hann","flattop","nuttall"]

ventana = "flattop"
N_points = 500
x_fft = rfftfreq(N_points, 1 / sample_rate)

# generación de un df para todos los WT y otro para todos los mutantes
# en el que cada columna es periodograma de ese gusano en concreto
# WT
periodograms_F_WT = pd.DataFrame()
periodograms_F_WT.insert(0, "Freq", x_fft)
periodograms_F_WT = periodograms_F_WT.set_index("Freq")
for g in set(df.Gusano):
    if g[0 : g.index(" ")] == "CONTROL":
        g_temp = detrend(df.Curvatura[df.Gusano == g].to_numpy(), axis=0)
        x_period, y_period = periodogram(
            g_temp,
            fs=sample_rate,
            nfft=N_points,
            axis=0,
            window=ventana,
            scaling="density",
        )
        periodograms_F_WT.insert(len(periodograms_F_WT.columns), g, y_period)
# MUT
periodograms_F_MUT = pd.DataFrame()
periodograms_F_MUT.insert(0, "Freq", x_fft)
periodograms_F_MUT = periodograms_F_MUT.set_index("Freq")
for g in set(df.Gusano):
    if g[0 : g.index(" ")] == "MUT":
        g_temp = detrend(df.Curvatura[df.Gusano == g].to_numpy(), axis=0)
        x_period, y_period = periodogram(
            g_temp,
            fs=sample_rate,
            nfft=N_points,
            axis=0,
            window=ventana,
            scaling="density",
        )
        periodograms_F_MUT.insert(len(periodograms_F_MUT.columns), g, y_period)
# agrupamiento de las las funciones según
funciones_analisis.agrupamiento_gusanos_fft(periodograms_F_WT, "CONTROL")
funciones_analisis.agrupamiento_gusanos_fft(periodograms_F_MUT, "MUT")


# %%%% plot Medias

plt.plot(x_fft, periodograms_F_WT.Media, linewidth=0.8, color="blue")
plt.fill_between(
    x_fft,
    (periodograms_F_WT.Media - periodograms_F_WT.SEM),
    (periodograms_F_WT.Media + periodograms_F_WT.SEM),
    color="blue",
    linewidth=0.1,
    alpha=0.2,
)
plt.plot(x_fft, periodograms_F_MUT.Media, linewidth=0.8, color="orange")
plt.fill_between(
    x_fft,
    (periodograms_F_MUT.Media - periodograms_F_MUT.SEM),
    (periodograms_F_MUT.Media + periodograms_F_MUT.SEM),
    color="orange",
    linewidth=0.1,
    alpha=0.2,
)

plt.xlabel("Period (Seg)")
plt.show()

# %%%% plot Suma

plt.plot(x_fft, periodograms_F_WT.Suma, linewidth=0.8, color="blue")
plt.plot(x_fft, periodograms_F_MUT.Suma, linewidth=0.8, color="orange")
plt.xlabel("Period (Seg)")
plt.show()

# %%%% plot Max

plt.plot(x_fft, periodograms_F_WT.Max, linewidth=0.8, color="blue")
plt.plot(x_fft, periodograms_F_MUT.Max, linewidth=0.8, color="orange")
plt.xlabel("Period (Seg)")
plt.show()


# %%% Lomb-Scargle Periodograma conteniendo todos los gusanos

# Periodos a explorar con lombscargle, el segundo parámetro es el valor máximo
periods = np.linspace(
    2 / sample_rate, 10, 5000
)  # el valor minimo viene dado de modo que
# el valor maximo de la frecuencia cumpla el teorema de Nyquist -> Sample_rate = 2*freq_max
f_periodogram_LS = 2 * np.pi / periods

# generación de un df para todos los WT y otro para todos los mutantes
# en el que cada columna es periodograma de ese gusano en concreto
# WT
periodograms_LS_WT = pd.DataFrame()
periodograms_LS_WT.insert(0, "W", f_periodogram_LS)
periodograms_LS_WT = periodograms_LS_WT.set_index("W")
for g in set(df.Gusano):
    if g[0 : g.index(" ")] == "CONTROL":
        g_temp = detrend(df.Curvatura[df.Gusano == g].to_numpy(), axis=0)
        t_temp = df["T_seg"][df.Gusano == g]
        temp_LS = lombscargle(t_temp, g_temp, f_periodogram_LS, normalize=True)
        periodograms_LS_WT.insert(len(periodograms_LS_WT.columns), g, temp_LS)
# MUT
periodograms_LS_MUT = pd.DataFrame()
periodograms_LS_MUT.insert(0, "W", f_periodogram_LS)
periodograms_LS_MUT = periodograms_LS_MUT.set_index("W")
for g in set(df.Gusano):
    if g[0 : g.index(" ")] == "MUT":
        g_temp = detrend(df.Curvatura[df.Gusano == g].to_numpy(), axis=0)
        t_temp = df["T_seg"][df.Gusano == g]
        temp_LS = lombscargle(t_temp, g_temp, f_periodogram_LS, normalize=True)
        periodograms_LS_MUT.insert(len(periodograms_LS_MUT.columns), g, temp_LS)
# agrupamiento de las las funciones según

funciones_analisis.agrupamiento_gusanos_fft(periodograms_LS_WT, "CONTROL")
funciones_analisis.agrupamiento_gusanos_fft(periodograms_LS_MUT, "MUT")


# %%%% plot Medias

plt.plot(periods, periodograms_LS_WT.Media, linewidth=0.8, color="blue")
plt.fill_between(
    periods,
    (periodograms_LS_WT.Media - periodograms_LS_WT.SEM),
    (periodograms_LS_WT.Media + periodograms_LS_WT.SEM),
    color="blue",
    linewidth=0.1,
    alpha=0.2,
)
plt.plot(periods, periodograms_LS_MUT.Media, linewidth=0.8, color="orange")
plt.fill_between(
    periods,
    (periodograms_LS_MUT.Media - periodograms_LS_MUT.SEM),
    (periodograms_LS_MUT.Media + periodograms_LS_MUT.SEM),
    color="orange",
    linewidth=0.1,
    alpha=0.2,
)

plt.xlabel("Period (Seg)")
plt.show()

# %%%% plot Suma

plt.plot(periods, periodograms_LS_WT.Suma, linewidth=0.8, color="blue")
plt.plot(periods, periodograms_LS_MUT.Suma, linewidth=0.8, color="orange")
plt.xlabel("Period (Seg)")
plt.show()

# %%%% plot Max

plt.plot(periods, periodograms_LS_WT.Max, linewidth=0.8, color="blue")
plt.plot(periods, periodograms_LS_MUT.Max, linewidth=0.8, color="orange")
plt.xlabel("Period (Seg)")
plt.show()
