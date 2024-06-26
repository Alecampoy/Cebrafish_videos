# -*- coding: utf-8 -*-
# Spyder Editor

# %% Intro [md]
"""
# **Análisis de Movimiento Zebrafish**
### Author: Alejandro Campoy Lopez
"""


# %% Librerias
from IPython import get_ipython

# limpia variables y librerias antiguas
get_ipython().magic("reset -sf")

import warnings
import pandas as pd
import numpy as np
import re
import os
import platform
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats, fft
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks, find_peaks_cwt, detrend, periodogram, lombscargle
from scipy.linalg import dft

# os.chdir("/home/ale/Documentos/GitHub/Cebrafish_videos/Experiment 1 - Contractions")
from functions_aux_analysis import *

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

if platform.system() == "Windows":
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
        csv.insert(0, "Batch", re.search("batch \d+", f.lower()).group(0))
        csv.insert(1, "Fenotype", re.search("_(KO\s?\d*|WT)", f.upper()).group(1))
        csv.insert(2, "Fish", "ZebraF_" + re.search("(\d+)(.lif)", f.lower()).group(1))
        df.append(csv)
        del (csv, f)

# raw strings (r' ') to avoid escape character issues.

# convierte la lista en un dataframe
df = pd.concat(df)
# renombro la columna con nombre repetido
df = df.rename(columns={"XM.1": "YM"})

# renombro KO a KO44 para el batch 6 y 7
df.loc[df.Fenotype == "KO", "Fenotype"] = "KO44"
df.loc[df.Fenotype == "KO 44", "Fenotype"] = "KO44"
df.loc[df.Fenotype == "KO 179", "Fenotype"] = "KO179"

df["Batch"] = pd.Categorical(
    df["Batch"],
    categories=["batch 6", "batch 7", "batch 8", "batch 9", "batch 10", "batch 11"],
    ordered=True,
)

pd.crosstab(index=df.Batch, columns=df.Fenotype)


# %% NAs [md]
"""
## Número de NAs por ZebraF
Visualizamos el número de Frames no segmentados apropiadamente por pez. Dado que no son demasiados, los imputo mediante interpolación Lineal.
"""

# %%% NAs Plot
NAs = (
    df.groupby(["Batch", "Fenotype", "Fish"])
    .apply(lambda x: x.isnull().sum())[["area"]]
    .rename(columns={"area": "NAs"})
).reset_index()  # me quedo solo con una columna ya que el numero de NAN es el mismo en todas

# NAs["Batch_Feno"] = NAs.Batch.astype(str) + "_" + NAs.Fenotype.astype(str)
# NAs_barplot = sns.barplot(x="Fish", y="NAs", hue="Batch", data=NAs.reset_index())
# plt.xticks(rotation=90)
# plt.show()

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
# df._get_numeric_data().columns
# Only numeric columns
df[df.select_dtypes(include=np.number).columns] = df[
    df.select_dtypes(include=np.number).columns
].interpolate(method="linear")

# %%% [md]
"""
El Zebra 10 WT del batch 7 se ha eliminado por contener > 500 NAs
El Zebra WT 1 tiene un video de la mitad de frames. Como esta totalmente quieto, se duplican sus datos para que se ajuste a la misma longitud de los demas"""

# %% Distancia Recorrida [md]
"""
## Distancia Recorrida
Se calcula la distancia que recorre el pez a lo largo del video y se gráfica por batch
"""

# %%% Calculo de la distancia

df.insert(7, "X_diff", df.groupby(["Batch", "Fenotype", "Fish"]).XM.diff())
df.insert(8, "Y_diff", df.groupby(["Batch", "Fenotype", "Fish"]).YM.diff())
df.insert(9, "dist", np.sqrt((df.X_diff**2) + (df.Y_diff**2)))

# dataframe con la distancia recorrida por el  gusano
Dist = (
    df.dropna()
    .groupby(["Batch", "Fenotype", "Fish"], as_index=False)
    .dist.sum(min_count=1)
    .round()
)
Dist = Dist.loc[
    Dist.dist != 0
]  # Aparecen como 0 las categorias para las que no hay datos, pero tambien hay zebra que no se mueven

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


grped_bplot.set_title("Distancia Recorrida por el Zebrafish durante el video (px)")
plt.legend(handles[0:3], labels[0:3])
plt.show()

# %%% Suavizado de las columnas [md]
"""
### Suavizado de las columnas
He observado que hay peces que vibran mucho, por lo que su medición muestra que se mueven mucho cuando apenas se han movido realmente. Para solucionarlo voy a aplicar un filtro de ventana gaussiana y recalcular el ultimo gráfico
"""

# %%%% Ventana Gausiana & calculo distancia suavizada


def apply_gaussian_filter(group_df, column, new_column_name, sigma=3.0):
    group_df[new_column_name] = gaussian_filter1d(group_df[column], sigma)
    return group_df


df = df.groupby(["Batch", "Fenotype", "Fish"], group_keys=False).apply(
    apply_gaussian_filter, column="XM", new_column_name="XM_filt"
)
df = df.groupby(["Batch", "Fenotype", "Fish"], group_keys=False).apply(
    apply_gaussian_filter, column="YM", new_column_name="YM_filt"
)

df["X_filt_diff"] = df.groupby(["Batch", "Fenotype", "Fish"]).XM_filt.diff()
df["Y_filt_diff"] = df.groupby(["Batch", "Fenotype", "Fish"]).YM_filt.diff()
df["dist_filt"] = np.sqrt((df.X_filt_diff**2) + (df.Y_filt_diff**2))

# %%%% Boxplot por batch Filtrado
# dataframe con la distancia recorrida por el  gusano
Dist_filt = (
    df.dropna()
    .groupby(["Batch", "Fenotype", "Fish"], as_index=False)
    .dist_filt.sum(min_count=1)
    .round()
)
Dist_filt = Dist_filt.loc[
    Dist_filt.dist_filt != 0
]  # Aparecen como 0 las categorias para las que no hay datos, pero tambien hay zebra que no se mueven


grped_bplot = sns.catplot(
    x="Batch",
    y="dist_filt",
    hue="Fenotype",
    kind="box",
    showfliers=False,
    legend=False,
    height=6,
    aspect=1.9,
    data=Dist_filt,
    hue_order=["WT", "KO44", "KO179"],
)
# make grouped stripplot
grped_bplot = sns.stripplot(
    x="Batch",
    y="dist_filt",
    hue="Fenotype",
    jitter=0.18,
    dodge=True,
    marker="o",
    color="black",
    # palette="Set2",
    data=Dist_filt,
    hue_order=["WT", "KO44", "KO179"],
)
handles, labels = grped_bplot.get_legend_handles_labels()


grped_bplot.set_title(
    "Distancia Suavizada Recorrida por el Zebrafish durante el video (px)"
)
plt.legend(handles[0:3], labels[0:3])
plt.show()

# %%% [md]
"""
La aplicación del filtro no cambia los resultados
"""
# %% Evolución temporal de todas las variables [md]
"""
## Evolución temporal de las variables

Dado que buscamos evaluar el comportamiento de los peces, vamos a representar
las variables que hemos extraido del análisis de los videos con el tiempo.


El estado basal del pez se considera estirado, así, cuando el pez realiza alguna acción, se reflejará como un cambios en las variables.
Las contracciones que esperamos, se reflejaran como picos en la evolución temporal.


Busco picos, por lo que las magnitudes son interesantes si presentan la linea basal baja, así que calculo la inversa de las que no la tienen.
"""
# %%% Inversa de algunas magnitudes
# Inversa de algunas magnitudes
df["area_inv"] = 1 / df.area
df["Perim_inv"] = 1 / df.Perim
df["LongestShortestPath_inv"] = 1 / df.LongestShortestPath
df["Feret_inv"] = 1 / df.Feret

# df['Prueba'] = (df.Circ + df.Round)*df.Solidity # composición de varias magnitudes

# %%% Grafico todas las magnitudes temporales
df_temp = df[
    (df.Batch == "batch 11") & (df.Fenotype == "KO44") & (df.Fish == "ZebraF_1")
].melt(
    id_vars=["Frame"],
    value_vars=[
        # "area",
        "area_inv",
        # "Perim",
        "Perim_inv",
        "Circ",
        "Round",
        # "AR",
        # "LongestShortestPath",
        "LongestShortestPath_inv",
        # "Feret",
        "Feret_inv",
        # "MinFeret",
        "Solidity",
    ],
)

g = sns.FacetGrid(
    df_temp,
    hue="variable",
    row="variable",
    sharex="col",
    sharey=False,
    height=5,
    aspect=4,
    margin_titles=True,
)
g.map(sns.lineplot, "Frame", "value")
g.set_axis_labels(fontsize=20)
g.fig.suptitle(
    "Evolución temporal de todas las variables para un pez de ejemplo",
    fontsize=24,
    fontdict={"weight": "bold"},
)
g.fig.subplots_adjust(top=0.97)
# sns.set(font_scale=2)

plt.show()

# %%% Correlación entre variables

#  Correlación entre variables

# sns.scatterplot(data=df[(df.Batch == "batch 11") & (df.Fenotype == "KO44") & (df.Fish == "ZebraF_1")], x="Round", y="Solidity", hue = "Frame")
# plt.show()

# Todas las señales correlacionan altamente, esto se podrá comprobar con AFC() o con gringer

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
threshold = 0.85
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


# %%% Evolución del resultado vs threshold [md]
"""
## Evolución del resultado con el threshold
Dado que este resultado es sensible al Threshold, vamos a ver como cambia el resultado con el Threshold
elegido. Se representa la diferencia de la mediana por batch del tiempo que pasa replegado el KO con respecto a su mutante.
(Este resultado puede ser sensible a la normalización de la señal que queda aún pendiente.)

### Solidity Plot
"""


# %%% Construcción del DF de evolución del resultado con el threshold (NUEVO CODIGO)

Variable_plot = "Solidity"

threshold_result = pd.DataFrame(
    columns=["Threshold", "Batch", "Fenotype", "Mean_diff", "CI"]
)
ref = ko44 = ko179 = np.nan

for thr in np.arange(0.1, 1.01, 0.01):  # iteración sobre el threshold
    # data frame con los valores para ese threshold
    time_over_Thr = (
        df.groupby(["Batch", "Fenotype", "Fish"])[Variable_plot]
        .agg(
            contracted=lambda x: (x > thr).sum(),
            contracted_perc=lambda x: 100 * (x > thr).sum() / len(x),
        )
        .reset_index()
        .dropna()
    )

    # generación de los resultados de ese dataframe y añadidos al dataframe de los resultados. La iteración se hace sobre los batch para calcular los valores para cada fenotipo
    for batch, group in time_over_Thr.groupby("Batch"):
        WT = (
            group.dropna()
            .drop("Batch", axis=1)
            .loc[group.Fenotype == "WT"]
            .contracted_perc
        )
        grouped_batch = group.dropna().drop("Batch", axis=1).groupby("Fenotype")
        for fenotype, group in grouped_batch:  # loop over each possible fenotype
            if fenotype != "WT":
                mut = group.contracted_perc
                new_row = {
                    "Threshold": thr,
                    "Batch": batch,
                    "Fenotype": fenotype,
                    "Mean_diff": np.mean(WT) - np.mean(mut),
                    "CI": np.ptp(
                        stats.ttest_ind(WT, mut).confidence_interval(
                            confidence_level=0.80
                        )
                    )
                    / 2,
                }
                threshold_result.loc[len(threshold_result)] = new_row


# %%% Plot del resultado frente al threshold con Intervalos de confianza

threshold_result["hue"] = threshold_result.Batch + " - " + threshold_result.Fenotype

df_plot = threshold_result[
    threshold_result["Fenotype"].isin(["KO44"])
]  # filtro para lo que se quiere repesentar
df_plot["CI_up"] = df_plot.Mean_diff + df_plot.CI
df_plot["CI_down"] = df_plot.Mean_diff - df_plot.CI

g = sns.lineplot(
    data=df_plot,
    x="Threshold",
    y="Mean_diff",
    hue="hue",
)
g.set_title("Diference of Solidity Batch Mean values with Threshold")

for hue in df_plot.hue.unique():
    g.fill_between(
        x="Threshold",
        y1="CI_up",
        y2="CI_down",
        alpha=0.1,
        data=df_plot[df_plot.hue == hue],
    )
plt.show()


# %%%  Circularity vs Threshold Plot [md]
"""
### Circularity plot
Lo mismo para la circularity
"""

# %%%% Construcción del DF

Variable_plot = "Circ"

threshold_result = pd.DataFrame(
    columns=["Threshold", "Batch", "Fenotype", "Mean_diff", "CI"]
)
ref = ko44 = ko179 = np.nan

for thr in np.arange(0.1, 1.01, 0.01):  # iteración sobre el threshold
    # data frame con los valores para ese threshold
    time_over_Thr = (
        df.groupby(["Batch", "Fenotype", "Fish"])[Variable_plot]
        .agg(
            contracted=lambda x: (x > thr).sum(),
            contracted_perc=lambda x: 100 * (x > thr).sum() / len(x),
        )
        .reset_index()
        .dropna()
    )

    # generación de los resultados de ese dataframe y añadidos al dataframe de los resultados. La iteración se hace sobre los batch para calcular los valores para cada fenotipo
    for batch, group in time_over_Thr.groupby("Batch"):
        WT = (
            group.dropna()
            .drop("Batch", axis=1)
            .loc[group.Fenotype == "WT"]
            .contracted_perc
        )
        grouped_batch = group.dropna().drop("Batch", axis=1).groupby("Fenotype")
        for fenotype, group in grouped_batch:  # loop over each possible fenotype
            if fenotype != "WT":
                mut = group.contracted_perc
                new_row = {
                    "Threshold": thr,
                    "Batch": batch,
                    "Fenotype": fenotype,
                    "Mean_diff": np.mean(WT) - np.mean(mut),
                    "CI": np.ptp(
                        stats.ttest_ind(WT, mut).confidence_interval(
                            confidence_level=0.80
                        )
                    )
                    / 2,
                }
                threshold_result.loc[len(threshold_result)] = new_row

# %%%% Plot

threshold_result["hue"] = threshold_result.Batch + " - " + threshold_result.Fenotype

df_plot = threshold_result[
    threshold_result["Fenotype"].isin(["KO44"])
]  # filtro para lo que se quiere repesentar
df_plot["CI_up"] = df_plot.Mean_diff + df_plot.CI
df_plot["CI_down"] = df_plot.Mean_diff - df_plot.CI

g = sns.lineplot(
    data=df_plot,
    x="Threshold",
    y="Mean_diff",
    hue="hue",
)
g.set_title("Diference of Solidity Batch Mean values with Threshold")

for hue in df_plot.hue.unique():
    g.fill_between(
        x="Threshold",
        y1="CI_up",
        y2="CI_down",
        alpha=0.1,
        data=df_plot[df_plot.hue == hue],
    )
plt.show()


# %%%% [md]
"""
Como es de esperar para la circularity, y debido a la alta correlación entre variables, el efecto es el mismo, pero incluso más claro. Los intervalos de confiza son consistentemente diferentes, lo que indica una alta variabilidad entre batches

"""

# %% Peaks - Numero de coletazos [md]
"""
# Peaks - Numero de coletazos

Contando el número de picos de las señales anteriores se evalua el número de coletazos que ejecuta el pez en el tiempo del video.
Para ello usamos la funcion peak finder sobre la magnitud que parece que muestra un mayor rango o SNR, Roundness

### Ejemplo de Peak Finder sobre un pez
"""

# %%% Peak Finder - Ejemplo 1 pez
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

# %%% Peak finder a todos los Zebra

# df["unique_fish"] = (
#     df.Batch.astype(str) + "_" + df.Fenotype.astype(str) + "_" + df.Fish.astype(str)
# )

# # Filtro para ver por batch
# # dfa = df[df.Batch == "batch 6"]

# for f in sorted(set(df.unique_fish)):
#     fish_temp = df[df.unique_fish == f].Round  # ajustar estos parametros
#     peaks, _ = find_peaks(
#         fish_temp, height=0.4, prominence=0.08, threshold=0.0, distance=2, width=1
#     )

#     plt.plot(fish_temp)
#     plt.plot(peaks, fish_temp[peaks], "2", markersize=24)
#     plt.title(f)
#     plt.show()

# %%% [md]
"""
Cuando visualizo las gráficas de todos los Zebra, todos los picos están bien encontrados.

Aplicar un filtro a la gráfica del movimiento NO es necesario, pues los peaks estan bien encontrados
"""
# %%% [md]
"""
## Número de picos por condición

Represento el número de picos por condición y batch. Dado que todos los videos duran el mismo tiempo, se puede asociar a la frecuencia.

### Roundness

Usando la Roundness
"""

# %%% Por condición

Variable_plot = "Round"
peaks_df = (
    df.groupby(["Batch", "Fenotype", "Fish"])[Variable_plot]
    .apply(
        lambda x: len(
            find_peaks(
                x,
                height=0.4,
                prominence=0.08,
                threshold=0.0,
                distance=2,
                width=1,  # para magnitudes 0-1
                # x,height=0.005,prominence=0.002,threshold=0.0,distance=2,width=1,  # para perimetro_inv
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

El método funciona correctamente, pero hay mucha variabilidad interbatch. Aunque los WT son similares en la mayoria de batches.
Recomiendo repasar las gráficas de los peces comparandolas con las fotos de microscopia y ver que videos contienen defectos. También es necesario aumentar la N
"""

# %% Intro FFT y Periodograma [md]
"""
# Fourier Transform & Periodogram

En lugar de contar el número de picos, voy a usar transformar las señales al dominios de las frecuencias. Con esto busco encontrar frecuencias más fuertes para alguna condición. Usaré la transformada de Fourier para encontrar las frecuencias más fuertes en el dominio de frecuencias y lo traspasaremos a periodos. Repetiré el estudio usando el Periodograma de XXXX. Complementaré el estudio tratando de usar una Wavelet Transformation para convertir los picos en una señal más sinosoidal.

"""
# %% FFT [md]
"""
# Fourier Transform & Periodogram

En lugar de contar el número de picos, voy a usar transformar las señales al dominio de las frecuencias. Con esto busco encontrar frecuencias más fuertes para alguna condición experimental

## Ejemplo FFT 1 pez
"""

# %%% Zebra Ejemplo FFT - Analisis de Frecuencias del movimiento para un Pez Zebra
# la duracion debe estar ajustada para cada sample por haber eliminado los NA. comprobar la duración del video y calular como se eliminan los NA: mejor imputarlos que borrarlos. Los videos tienen 1550 frames

zebra_temp = df[
    (df.Batch == "batch 10") & (df.Fenotype == "KO179") & (df.Fish == "ZebraF_6")
]


signal = zebra_temp.Round.values  # signal as array
# smooth signal

# to avoid the signal at the first fourier coefficient F(0) we should substract offset
# ^f(0) =  mean power of signal
# signal = signal - np.mean( signal)
signal = detrend(signal, axis=0)

sample_rate = 9  # frames / s
time_step = 1 / sample_rate  # Delta t
N_points = len(signal)  # lenght signal
# t = np.arange(0, N_points/sample_rate, time_step)
time_points = zebra_temp.Time


def plot_fft_filter(signal, sample_rate, time_step, N_points, time_points, f="Zebra"):
    # FFT
    zebra_fft = fft.fft(signal, norm="backward")
    # Calculos para frequencias en seg pare representar en el eje X
    zebra_freqs = fft.fftfreq(N_points, time_step)
    # Power spectral density
    zebra_psd = np.abs(zebra_fft) ** 2 / (
        sample_rate * N_points
    )  # Power spectral density

    # Filtrado FFT - Filtro por la psd
    zebra_fft_fil = np.array(zebra_fft)
    # # upper power limit
    # zebra_fft_fil[zebra_psd > 0.030] = 0
    # # lower power limit
    zebra_fft_fil[zebra_psd < 0.05e-6] = 0
    # bandpass filter
    zebra_fft_fil[abs(zebra_freqs > 0.5)] = 0
    zebra_fil_psd = np.abs(zebra_fft_fil) ** 2 / (sample_rate * N_points)

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(21, 13))
    fft_plot_x = zebra_freqs[
        (zebra_freqs < 3) & (zebra_freqs > 0)
    ]  # set max freq to plot

    sns.lineplot(x=time_points, y=signal, ax=axes[0, 0], color="r")
    sns.lineplot(x=fft_plot_x, y=zebra_psd[0 : len(fft_plot_x)], ax=axes[0, 1])
    sns.lineplot(x=fft_plot_x, y=zebra_fil_psd[0 : len(fft_plot_x)], ax=axes[1, 0])
    # same limit for the filtered fft
    axes[1, 0].set_ylim(axes[0, 1].get_ylim())
    sns.lineplot(
        x=time_points, y=signal, linewidth=1, alpha=0.6, ax=axes[1, 1], color="r"
    )
    sns.lineplot(x=time_points, y=fft.ifft(zebra_fft_fil), ax=axes[1, 1])

    axes[0, 0].set_title("Signal")
    axes[0, 1].set_title("FFT PSD")
    axes[1, 0].set_title("Filtered FFT")
    axes[1, 1].set_title("IFFT")
    fig.suptitle(f)
    plt.tight_layout()
    plt.show()


plot_fft_filter(signal, sample_rate, time_step, N_points, time_points)

# %%% FFT filter for every zebra

df["unique_fish"] = (
    df.Batch.astype(str) + "_" + df.Fenotype.astype(str) + "_" + df.Fish.astype(str)
)

# Filtro para ver por batch
dfa = df[df.Batch == "batch 7"]

for f in sorted(set(dfa.unique_fish)):
    signal = dfa[dfa.unique_fish == f].Round.values  # ajustar estos parametros
    signal = detrend(signal, axis=0)
    sample_rate = 9  # frames / s
    time_step = 1 / sample_rate
    N_points = len(signal)
    # t = np.arange(0, N_points/sample_rate, time_step)
    time_points = zebra_temp.Time
    plot_fft_filter(signal, sample_rate, time_step, N_points, time_points, f)

# %% CODIGO GUSANOS


# %%%% FFT DF para todos los gusanos
# coimprobar como funciona con el codigo de los gusanos para ver si interesa replicarlo
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


# %% Lomb-Scargle Periodogram
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
