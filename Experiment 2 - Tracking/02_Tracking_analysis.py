# -*- coding: utf-8 -*-
# %%
# Spyder Editor

# %% Intro [md]
"""
# Tracking Experiment
# **Análisis del Tracking del Movimiento de Zebrafish en una placa**
### Author: Alejandro Campoy Lopez  
"""

# %% Librerias
import warnings
import pandas as pd
import re
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def get_files_in_folder(folder_path):
    file_list = glob.glob(folder_path + "/**", recursive=True)
    files = [file for file in file_list if not os.path.isdir(file)]
    return files


plt.rcParams["figure.figsize"] = (15, 8)

# %matplotlib inline

warnings.filterwarnings("ignore")

# %% Lectura Archivos [md]
"""
# Lectura de Archivos
Lectura de todos los archivos csv con los resultados de los diferentes batches.
Se añade una columna representando el gusano y el batch mediante el uso de regex
"""

# %%% Load Files

windows = True
if windows:
    folder_path = "p:\\CABD\\Lab Ozren\\Marta Fernandez\\Experimento Tracking\\resultados sucios\\"
else:
    folder_path = "/home/ale/pCloudDrive/CABD/Lab Ozren/Marta Fernandez/Experimento Tracking/resultados sucios/"
files = get_files_in_folder(folder_path)

df = []
for f in files:
    if ".csv" in f:
        csv = pd.read_csv(f, sep=";")
        csv.insert(0, "Batch", re.search("batch \d+", f.lower()).group(0))
        csv.insert(1, "Fenotype", re.search("(KO\d*|WT)", f.upper()).group(1))
        if csv.Batch.iloc[1] == "batch 7":
            csv.insert(2, "Fish", "ZebraF_" + re.search("(Ko |Wt )(\d+)", f).group(2))
        else:
            csv.insert(
                2, "Fish", "ZebraF_" + re.search("(\d+)(.mp4)", f.lower()).group(1)
            )
        df.append(csv)
        del (csv, f)

df = pd.concat(df)
df = df.drop(df.columns[-1], axis=1)


# renombro KO a KO44 para el batch 6 y 7
df.loc[df.Fenotype == "KO", "Fenotype"] = "KO44"
df.loc[df.Fenotype == "KO 44", "Fenotype"] = "KO44"
df.loc[df.Fenotype == "KO 179", "Fenotype"] = "KO179"

df["Batch"] = pd.Categorical(
    df["Batch"],
    categories=["batch 6", "batch 7", "batch 8", "batch 11"],
    ordered=True,
)

elements = round(
    pd.crosstab(index=df.Batch, columns=df.Fenotype) / 4202
)  # divided by lengh of the video
print(str(elements).replace(".0", "").replace("],", "]\n"))

# %% NAs [md]
"""
## Número de NAs por ZebraF
Visualizamos el número de Frames no segmentados apropiadamente por pez. Dado que no son demasiados, los imputo mediante interpolación Lineal.
"""

# %%% NAs Plot
NAs = (
    df.groupby(["Batch", "Fenotype", "Fish"])
    .apply(lambda x: x.isnull().sum())[["Mean-Distance"]]
    .rename(columns={"Mean-Distance": "NAs"})
).reset_index()  # me quedo solo con una columna ya que el numero de NAN es el mismo en todas

NAs_barplot = sns.catplot(
    kind="bar",
    data=NAs.reset_index(),
    x="Fish",
    y="NAs",
    col="Fenotype",
    row="Batch",
    legend=True,
)

NAs_barplot.set_xticklabels(rotation=45, size=333)
NAs_barplot.set_xlabels("Fish", fontsize=15)
plt.show()

# %%% NA Impute

df[['X', 'Y', 'Mean-Distance']] = df[['X', 'Y',
                                      'Mean-Distance']].interpolate(method="linear")

# %%% [md]
"""
No hay NAs
"""
# %% Distancia Recorrida [md]
"""
## Distancia Recorrida
Se calcula la distancia que recorre el pez a lo largo del video y se gráfica por batch
"""

# %%% Calculo de la distancia recorrida

df.insert(7, "X_diff", df.groupby(["Batch", "Fenotype", "Fish"]).X.diff())
df.insert(8, "Y_diff", df.groupby(["Batch", "Fenotype", "Fish"]).Y.diff())
df.insert(9, "dist", np.sqrt((df.X_diff**2) + (df.Y_diff**2)))

# dataframe con la distancia recorrida por el  gusano
Dist = (
    df.dropna()
    .groupby(["Batch", "Fenotype", "Fish"], as_index=False)
    .dist.sum(min_count=1)
    .round()
)  # .reset_index()

Dist = Dist.loc[Dist.dist != 0]
# %%% Box-plot por batch

grped_bplot = sns.catplot(
    x="Batch",
    y="dist",
    hue="Fenotype",
    kind="box",
    showfliers=False,
    legend=False,
    height=6,
    aspect=1.9,
    data=Dist,
    hue_order=["WT", "KO44", "KO179"],
)
# make grouped stripplot
grped_bplot = sns.stripplot(
    x="Batch",
    y="dist",
    hue="Fenotype",
    jitter=0.18,
    dodge=True,
    marker="o",
    color="black",
    # palette="Set2",
    data=Dist,
    hue_order=["WT", "KO44", "KO179"],
)
handles, labels = grped_bplot.get_legend_handles_labels()


grped_bplot.set_title("Distancia Total Recorrida por el Zebrafish (px)")
plt.legend(handles[0:3], labels[0:3])
plt.show()

# %% Distribución de la posición del Pez [md]
""" 
# Distribución de la posición del Zebra
A lo largo del video, el pez se posiciona en algún lugar de la placa. Se estima que la posición que mantiene el Zebra es comportamental, por lo que vamos a estudiar que posición mantiene con respecto al borde d=0 hasta el centro d=255. (Dado que existe simetria radial, solo usamos la distancia respecto al borde)
"""

# %%% Histograma 1 Zebra [md]
"""
## Histograma de 1 solo Pez
Vamos a evaluar el histograma de 1 Zebra. Este nos va a indicar donde se posiciona el Zebra a lo largo del tiempo del video.
Hay que normalizar el histograma
"""
# %%% Grafico Histograma 1 Zebra
df_temp = df[
    (df.Batch == "batch 7") & (df.Fenotype == "WT") & (df.Fish == "ZebraF_1")
]

sns.histplot(data=df_temp, x="Mean-Distance")
plt.show()

# %%% [md]
'''
Se observa como el Zebra se posiciona mayormente a lo largo de la franja centra, como esperable probabilisticamente, sin posicionarse apenas cerca del borde. 
'''

# %%% [md]
"""
Los picos en las señales correlacionan altamente, se ejecuta el análisis solamente sobre una única.
"""

# %% Análisis Overview[md]
"""
# Análisis

Para modelar el comportamiento del pez voy a realizar 3 aproximaciones:

## - Tiempo que pasa replegado
Dado que observamos picos, se puede evaluar el tiempo total que pasa replegado usando un threshold. 
Se asocia el tiempo total replegado a todo el tiempo que el valor esta por encima del threshold.
  
## - Contado de picos
Usando la función Peak Finder detectamos picos en las señales, que se pueden tanto contar como cuantificar con parámetros
como altura (no de interes) o anchura

## - Periodograma
Usando la FFT ver las frecuencias intrinsicas de cada uno de los peces. (puede ser interesante buscar la baseline)
"""

# %% Tiempo replegado [md]
"""
# Tiempo Replegado

Usando la Solidity = area/convex area, si su valor es superior al Threshold, 
indica que el gusano esta replegado. Se pueden usar otras magnitudes acotadas entre 0-1 
"""

# %%% Plot Ejemplo threshold

df_temp = df[
    (df.Batch == "batch 7") & (df.Fenotype == "KO44") & (df.Fish == "ZebraF_1")
]

g = sns.lineplot(data=df_temp, x="Time", y="Solidity")
g.axhline(0.88, color="red")
g.set_title("Threshold on Solidity", size=25)
plt.show()

# %%% [md]
"""
Contando para cada gusano el total del tiempo que pasa sobre el Threshold, obtenemos
"""


# %%% Comparación usando un threshold fijo

Variable_plot = "Solidity"
threshold = 0.88
time_over_Thr = (
    df.groupby(["Batch", "Fenotype", "Fish"])[Variable_plot]
    .apply(lambda x: (x > threshold).sum())
    .reset_index()
    .rename(columns={Variable_plot: "contracted"})
)

time_over_Thr["contracted_perc"] = 100 * time_over_Thr.contracted / 1550

# a = sns.boxplot(x="Fenotype", y="contracted_perc", data=time_over_Thr)
# a.set_title("Numero de tiempo replegado con Thr " + str(threshold))
# b = sns.stripplot(
#     x="Fenotype", y="contracted_perc", data=time_over_Thr, color="grey", size=8
# )
# plt.show()

grped_bplot = sns.catplot(
    x="Batch",
    y="contracted_perc",
    data=time_over_Thr,
    hue="Fenotype",
    kind="box",
    legend=False,
    showfliers=False,
    height=6,
    aspect=1.9,
    hue_order=["WT", "KO44", "KO179"],
)
# make grouped stripplot
grped_bplot = sns.stripplot(
    x="Batch",
    y="contracted_perc",
    data=time_over_Thr,
    hue="Fenotype",
    jitter=True,
    dodge=True,
    marker="o",
    color="black",
    # palette="Set2",
    hue_order=["WT", "KO44", "KO179"],
)
handles, labels = grped_bplot.get_legend_handles_labels()
grped_bplot.set_title(
    "Porcentaje del tiempo que pasa el gusano replegado - sobre el Threshold = "
    + str(threshold),
    size=20,
)
plt.legend(handles[0:3], labels[0:3])
plt.show()


# %%% Evolución del resultado con el threshold [md]
"""
## Evolución del resultado con el threshold
Dado que este resultado es sensible al Threshold, vamos a ver como cambia el resultado con el Threshold
elegido. Se representa la diferencia de la mediana por batch del tiempo que pasa replegado el KO con respecto a su mutante. 
(Este resultado puede ser sensible a la normalización de la señal que queda aún pendiente.)

### Solidity Plot
"""

# %%% plot del threshold para Solidity
threshold_result = pd.DataFrame(columns=["Threshold", "Batch", "KO44", "KO179"])

i = 0
Variable_plot = "Solidity"
for thr in np.arange(0.1, 1.01, 0.01):
    time_over_Thr = (
        df.groupby(["Batch", "Fenotype", "Fish"])[Variable_plot]
        .apply(lambda x: (x > thr).sum())
        .reset_index()
        .rename(columns={Variable_plot: "contracted"})
    )
    time_over_Thr["contracted_perc"] = 100 * time_over_Thr.contracted / 1550

    df_temp_median = time_over_Thr.groupby(["Batch", "Fenotype"])[
        "contracted_perc"
    ].median()
    threshold_result.loc[i] = [
        thr,
        "batch 6",
        df_temp_median["batch 6", "WT"] - df_temp_median["batch 6", "KO44"],
        np.nan,
    ]
    threshold_result.loc[i + 1] = [
        thr,
        "batch 7",
        df_temp_median["batch 7", "WT"] - df_temp_median["batch 7", "KO44"],
        np.nan,
    ]
    threshold_result.loc[i + 2] = [
        thr,
        "batch 8",
        df_temp_median["batch 8", "WT"] - df_temp_median["batch 8", "KO44"],
        df_temp_median["batch 8", "WT"] - df_temp_median["batch 8", "KO179"],
    ]
    i = i + 3


df_temp = threshold_result.melt(id_vars=["Threshold", "Batch"]).dropna()
df_temp["hue"] = df_temp.Batch + " - " + df_temp.variable
g = sns.lineplot(data=df_temp, x="Threshold", y="value", hue="hue")
g.set_title("Diference of Solidity Batch Median values with Threshold")
plt.show()

# %%% plot del threshold para Solidity, más batches
threshold_result = pd.DataFrame(columns=["Threshold", "Batch", "KO44", "KO179"])

i = 0
Variable_plot = "Solidity"
for thr in np.arange(0.1, 1.01, 0.01):
    time_over_Thr = (
        df.groupby(["Batch", "Fenotype", "Fish"])[Variable_plot]
        .apply(lambda x: (x > thr).sum())
        .reset_index()
        .rename(columns={Variable_plot: "contracted"})
    )
    time_over_Thr["contracted_perc"] = 100 * time_over_Thr.contracted / 1550

    df_temp_median = time_over_Thr.groupby(["Batch", "Fenotype"])[
        "contracted_perc"
    ].median()
    threshold_result.loc[i] = [
        thr,
        "batch 9",
        df_temp_median["batch 9", "WT"] - df_temp_median["batch 9", "KO44"],
        df_temp_median["batch 9", "WT"] - df_temp_median["batch 9", "KO179"],
    ]
    threshold_result.loc[i + 1] = [
        thr,
        "batch 10",
        np.nan,
        df_temp_median["batch 10", "WT"] - df_temp_median["batch 10", "KO179"],
    ]
    threshold_result.loc[i + 2] = [
        thr,
        "batch 11",
        df_temp_median["batch 11", "WT"] - df_temp_median["batch 11", "KO44"],
        df_temp_median["batch 11", "WT"] - df_temp_median["batch 11", "KO179"],
    ]
    i = i + 3


df_temp = threshold_result.melt(id_vars=["Threshold", "Batch"]).dropna()
df_temp["hue"] = df_temp.Batch + " - " + df_temp.variable
g = sns.lineplot(data=df_temp, x="Threshold", y="value", hue="hue")
g.set_title("Diference of Solidity Batch Median values with Threshold 2")
plt.show()

# %%% [md]
"""
Parece que el KO44 y el KO179 se comportan igual en el batch 8, pero los batches 6 y 7 tiene el KO44 un comportamiento opuesto

### Circularity plot
Lo mismo para la circularity
"""

# %%% plot del threshold para Circularity
threshold_result = pd.DataFrame(columns=["Threshold", "Batch", "KO44", "KO179"])

i = 0
Variable_plot = "Circ"
for thr in np.arange(0.2, 1.01, 0.01):
    time_over_Thr = (
        df.groupby(["Batch", "Fenotype", "Fish"])[Variable_plot]
        .apply(lambda x: (x > thr).sum())
        .reset_index()
        .rename(columns={Variable_plot: "contracted"})
    )
    time_over_Thr["contracted_perc"] = 100 * time_over_Thr.contracted / 1550

    df_temp_median = time_over_Thr.groupby(["Batch", "Fenotype"])[
        "contracted_perc"
    ].median()
    threshold_result.loc[i] = [
        thr,
        "batch 6",
        df_temp_median["batch 6", "WT"] - df_temp_median["batch 6", "KO44"],
        np.nan,
    ]
    threshold_result.loc[i + 1] = [
        thr,
        "batch 7",
        df_temp_median["batch 7", "WT"] - df_temp_median["batch 7", "KO44"],
        np.nan,
    ]
    threshold_result.loc[i + 2] = [
        thr,
        "batch 8",
        df_temp_median["batch 8", "WT"] - df_temp_median["batch 8", "KO44"],
        df_temp_median["batch 8", "WT"] - df_temp_median["batch 8", "KO179"],
    ]
    i = i + 3


df_temp = threshold_result.melt(id_vars=["Threshold", "Batch"]).dropna()
df_temp["hue"] = df_temp.Batch + " - " + df_temp.variable
g = sns.lineplot(data=df_temp, x="Threshold", y="value", hue="hue")
g.set_title("Diference of " + Variable_plot + "Batch Median values with Threshold")
plt.show()

# %%% plot del threshold para Circularity 2
threshold_result = pd.DataFrame(columns=["Threshold", "Batch", "KO44", "KO179"])

i = 0
Variable_plot = "Circ"
for thr in np.arange(0.2, 1.01, 0.01):
    time_over_Thr = (
        df.groupby(["Batch", "Fenotype", "Fish"])[Variable_plot]
        .apply(lambda x: (x > thr).sum())
        .reset_index()
        .rename(columns={Variable_plot: "contracted"})
    )
    time_over_Thr["contracted_perc"] = 100 * time_over_Thr.contracted / 1550

    df_temp_median = time_over_Thr.groupby(["Batch", "Fenotype"])[
        "contracted_perc"
    ].median()
    threshold_result.loc[i] = [
        thr,
        "batch 9",
        df_temp_median["batch 9", "WT"] - df_temp_median["batch 9", "KO44"],
        df_temp_median["batch 9", "WT"] - df_temp_median["batch 9", "KO179"],
    ]
    threshold_result.loc[i + 1] = [
        thr,
        "batch 10",
        np.nan,
        df_temp_median["batch 10", "WT"] - df_temp_median["batch 10", "KO179"],
    ]
    threshold_result.loc[i + 2] = [
        thr,
        "batch 11",
        df_temp_median["batch 11", "WT"] - df_temp_median["batch 11", "KO44"],
        df_temp_median["batch 11", "WT"] - df_temp_median["batch 11", "KO179"],
    ]
    i = i + 3


df_temp = threshold_result.melt(id_vars=["Threshold", "Batch"]).dropna()
df_temp["hue"] = df_temp.Batch + " - " + df_temp.variable
g = sns.lineplot(data=df_temp, x="Threshold", y="value", hue="hue")
g.set_title("Diference of " + Variable_plot + "Batch Median values with Threshold 2")
plt.show()

# %%% [md]
"""
Como es de esperar, y debido a la alta correlación entre variables, el efecto es el mismo

"""

# %% Peaks - Numero de coletazos [md]
"""
# Peaks - Numero de coletazos

Contando el número de picos de las señales anteriores se evalua el número de coletazos que ejecuta el pez en el tiempo del video.
Para ello usamos la funcion peak finder sobre la magnitud que parece que muestra un mayor rango o SNR, Roundness

### Ejemplo de Peak Finder sobre un pez
"""

# %%% Mediante Peak Finder
df_temp = df[
    (df.Batch == "batch 7") & (df.Fenotype == "KO44") & (df.Fish == "ZebraF_1")
].Round

peaks, _ = find_peaks(
    df_temp, height=0.4, prominence=0.08, threshold=0.0, distance=2, width=1
)  # ajustar estos parametros
# he comprobado que algun gusano realiza contracciones y extensiones en 2 frames, por lo que distance = 1 & width = 1

plt.plot(df_temp)
plt.plot(peaks, df_temp[peaks], "2", markersize=24)
plt.title("Picos encontrados sobre la Roundness", size=20)
plt.show()

# %%% [md]
"""
Es interesante ver como funciona sobre todos los gusanos
"""

# %%% Peak finder en todos los gusanos

df["unique_fish"] = (
    df.Batch.astype(str) + "_" + df.Fenotype.astype(str) + "_" + df.Fish.astype(str)
)

# Filtro para ver por batch
# dfa = df[df.Batch == "batch 6"]

for f in sorted(set(df.unique_fish)):
    fish_temp = df[df.unique_fish == f].Round  # ajustar estos parametros
    peaks, _ = find_peaks(
        fish_temp, height=0.4, prominence=0.08, threshold=0.0, distance=2, width=1
    )

    plt.plot(fish_temp)
    plt.plot(peaks, fish_temp[peaks], "2", markersize=24)
    plt.title(f)
    plt.show()

# %%% Aplicando un filtro

# Aplicar un filtro No es necesario, pues los peaks estan bien encontrados

# %%% [md]
"""
## Número de picos por condición

Represento el número de picos por condición y batch. Dado que todos los videos duran el mismo tiempo, se puede asociar a la frecuencia.

### Circularity

Usando la circularity
"""

# %%% Por condición

Variable_plot = "LongestShortestPath_inv"
peaks_df = (
    df.groupby(["Batch", "Fenotype", "Fish"])[Variable_plot]
    .apply(
        lambda x: len(
            find_peaks(
                # x, height=0.4, prominence=0.08, threshold=0.0, distance=2, width=1 # para magnitudes 0-1
                x,
                height=0.005,
                prominence=0.002,
                threshold=0.0,
                distance=2,
                width=1,  # para perimetro_inv
            )[0]
        )
    )
    .reset_index()
    .rename(columns={Variable_plot: "N_peaks"})
)

grped_bplot = sns.catplot(
    x="Batch",
    y="N_peaks",
    hue="Fenotype",
    kind="box",
    legend=False,
    showfliers=False,
    height=6,
    aspect=1.9,
    data=peaks_df,
    hue_order=["WT", "KO44", "KO179"],
)
# make grouped stripplot
grped_bplot = sns.stripplot(
    x="Batch",
    y="N_peaks",
    hue="Fenotype",
    jitter=0.18,
    dodge=True,
    marker="o",
    color="black",
    # palette="Set2",
    data=peaks_df,
    hue_order=["WT", "KO44", "KO179"],
)
handles, labels = grped_bplot.get_legend_handles_labels()
l = plt.legend(handles[0:3], labels[0:3])
grped_bplot.set_title("Number of peaks using " + Variable_plot)
# plt.figure(figsize=(19,8))
plt.show()


# %%% [md]
"""
He visualizado la gráfica anterior cambiando el valor de `height` en el intervalo 0.1-0.9 y 
no cambia hasta 0.8, el cual es ya un valor extremo para Roundness. Lo mismo con `prominence`en el 
intervalo 0.01-0,6 y es invariante hasta valores extremos a partir de 0.4 (altura del pico)

### Conclusión

El método funciona correctamente, pero hay mucha variabilidad interbatch. 
Recomiendo repasar las gráficas de los peces comparandolas con las fotos de microscopia y ver que videos contienen defectos. También es necesario aumentar la N
"""

# %% Intro FFT y Periodograma [md]
"""
# Fourier Transform & Periodogram

En lugar de contar el número de picos, voy a usar transformar las señales al dominios de las frecuencias. Con esto busco encontrar frecuencias más fuertes para alguna condición

"""
# %% FFT [md]
"""
# Fourier Transform & Periodogram

En lugar de contar el número de picos, voy a usar transformar las señales al dominios de las frecuencias. Con esto busco encontrar frecuencias más fuertes para alguna condición

"""

# %%% FFT Analisis de Frecuencias del movimiento para un Pez Zebra
# Ajustar al cebra
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

# %% CODIGO GUSANOS


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
    if g[0: g.index(" ")] == "CONTROL":
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
    if g[0: g.index(" ")] == "MUT":
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
    if g[0: g.index(" ")] == "CONTROL":
        x_temp = df[["Curvatura"]][df.Gusano == g][:limite_t]
        g_fft = rfft(detrend(x_temp, axis=0), axis=0, norm="forward")
        fft_WT.insert(len(fft_WT.columns), g, np.abs(g_fft))
# MUT
fft_MUT = pd.DataFrame()
fft_MUT.insert(0, "Freq", x_fft)
fft_MUT = fft_MUT.set_index("Freq")
for g in set(df.Gusano):
    if g[0: g.index(" ")] == "MUT":
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
    if g[0: g.index(" ")] == "CONTROL":
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
    if g[0: g.index(" ")] == "MUT":
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
    if g[0: g.index(" ")] == "CONTROL":
        g_temp = detrend(df.Curvatura[df.Gusano == g].to_numpy(), axis=0)
        t_temp = df["T_seg"][df.Gusano == g]
        temp_LS = lombscargle(t_temp, g_temp, f_periodogram_LS, normalize=True)
        periodograms_LS_WT.insert(len(periodograms_LS_WT.columns), g, temp_LS)
# MUT
periodograms_LS_MUT = pd.DataFrame()
periodograms_LS_MUT.insert(0, "W", f_periodogram_LS)
periodograms_LS_MUT = periodograms_LS_MUT.set_index("W")
for g in set(df.Gusano):
    if g[0: g.index(" ")] == "MUT":
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