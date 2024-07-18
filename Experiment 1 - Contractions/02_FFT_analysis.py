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
import os, time
import platform
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats, fft
from scipy.ndimage import gaussian_filter1d, median_filter
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
        csv.insert(
            2, "Fish", "ZebraF_" + re.search("(\d+)(.\wif_results)", f.lower()).group(1)
        )
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
    categories=[
        "batch 6",
        "batch 7",
        "batch 8",
        "batch 9",
        "batch 10",
        "batch 11",
        "batch 12",
        "batch 13",
        "batch 14",
    ],
    ordered=True,
)

(pd.crosstab(index=df.Batch, columns=df.Fenotype) / 1550)


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


def interpolate_group(group):
    return group.interpolate(method="linear", limit_direction="both")


# incluyo las categorías como indice para que no se aplique la interpolación a variables categoricas
df.set_index(["Batch", "Fenotype", "Fish"], inplace=True)

df = (
    df.groupby(["Batch", "Fenotype", "Fish"], as_index=False)
    .apply(interpolate_group)
    .droplevel(0)
)


# %%% [md]
"""
El Zebra 10 WT del batch 7 se ha eliminado por contener > 500 NAs
El Zebra WT 1 tiene un video de la mitad de frames. Como esta totalmente quieto, se duplican sus datos para que se ajuste a la misma longitud de los demas
"""

# %% Filtrado de peces anomalos o muertos [md]
"""
Encuentro que hay peces que apenas dan coletazos, de este modo su gráfica tiene ruido, pero no aplica a una apropiada señal de un coletazo. Para ver cuales son y eliminarlos voy a calcular su std y su rango, ya que estos serán anormalmente bajos
"""
# %%% STD y Rango de Circularity

STD_Rango = (
    df.groupby(["Batch", "Fenotype", "Fish"])
    .Round.agg(
        STD=lambda x: x.std(), Rango=lambda x: x.max() - x.min(), n_obs=lambda x: len(x)
    )
    .dropna()
)


# %%%% Plot Rango

grped_bplot = sns.catplot(
    x="Batch",
    y="Rango",
    hue="Fenotype",
    kind="box",
    showfliers=False,
    legend=True,
    height=6,
    aspect=1.9,
    data=STD_Rango,
    hue_order=["WT", "KO44", "KO179"],
)

# make grouped stripplot
grped_bplot = sns.stripplot(
    x="Batch",
    y="Rango",
    hue="Fenotype",
    jitter=0.18,
    dodge=True,
    legend=False,
    marker="o",
    color="black",
    # palette="Set2",
    data=STD_Rango,
    hue_order=["WT", "KO44", "KO179"],
)
# handles, labels = grped_bplot.get_legend_handles_labels()

grped_bplot.set_title("Rango de la Circularity")
# plt.legend(handles[0:3], labels[0:3])
plt.show()

# %%% Rango [md]

"""
En la grafica del rango se ven claramente los peces anómalos, recomiendo explorar el DF y eliminarlos. La STD no creo que aporte esta información, ya que los peces que se mueven poco la tendrán baja a pesar de ser funcionales
"""


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
    df.groupby(["Batch", "Fenotype", "Fish"])
    .dist.agg(dist=lambda x: x.sum(), n_obs=lambda x: len(x))
    .round()
    .dropna()
)

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
Hay un valor anomalo en el batch 7 WT zebra 7, mirarlo Tambien batch 7 zebra 4 WT y zebra 8 KO44

### Suavizado de las columnas
He observado que hay peces que vibran mucho, por lo que su medición muestra que se mueven mucho cuando apenas se han movido realmente. Para solucionarlo voy a aplicar un filtro de ventana gaussiana y recalcular el ultimo gráfico
"""

# %%%% Ventana Gausiana & calculo distancia suavizada


def apply_gaussian_filter(group_df, column, new_column_name, sigma=4.0):
    group_df[new_column_name] = gaussian_filter1d(group_df[column], sigma)
    return group_df


df = df.groupby(["Batch", "Fenotype", "Fish"], group_keys=False).apply(
    apply_gaussian_filter, column="XM", new_column_name="XM_filt"
)
df = df.groupby(["Batch", "Fenotype", "Fish"], group_keys=False).apply(
    apply_gaussian_filter, column="YM", new_column_name="YM_filt"
)

# nuevo calculo de la distancia recorrida
df["X_filt_diff"] = df.groupby(["Batch", "Fenotype", "Fish"]).XM_filt.diff()
df["Y_filt_diff"] = df.groupby(["Batch", "Fenotype", "Fish"]).YM_filt.diff()
df["dist_filt"] = np.sqrt((df.X_filt_diff**2) + (df.Y_filt_diff**2))

# %%%% Boxplot por batch Filtrado
# dataframe con la distancia recorrida por el  gusano
Dist_filt = (
    df.groupby(["Batch", "Fenotype", "Fish"])
    .dist_filt.agg(dist_filt=lambda x: x.sum(), n_obs=lambda x: len(x))
    .round()
    .dropna()
)

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
La aplicación del filtro no cambia los resultados pero si genera una visualización más adecuada. Mirar los extremos
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
df_temp = df.loc[("batch 11", "WT", "ZebraF_4")].melt(
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
Los picos en las señales correlacionan altamente, se ejecuta el análisis solamente sobre una única variable.
"""

# %%% Filtro Gaussiano o Median [md]
"""
Hay picos que debido al ruido, su parte alta tiene alguna irregularidad. También la linea basal es rugosa. Pienso que un filtro gausiano aplicado va a suavizar las irregularidades y facilitar el análisis posterior. Lo gráfico como ejemplo
"""

# %%% Filtro Gaussian
# lo aplico a las 3 señales que me interesan


def apply_median_filter(group_df, column, new_column_name, sigma=4.0):
    group_df[new_column_name] = median_filter(
        group_df[column], size=sigma, mode="reflect"
    )
    return group_df


df_temp = df.loc[("batch 11", "WT", "ZebraF_4")]
df_temp = apply_gaussian_filter(
    df_temp, column="Feret_inv", new_column_name="Feret_inv_filt", sigma=2
)
df_temp = apply_gaussian_filter(
    df_temp, column="Circ", new_column_name="Circ_filt", sigma=3
)
df_temp = apply_gaussian_filter(
    df_temp, column="Perim_inv", new_column_name="Perim_inv_filt", sigma=2
)

df_temp = df_temp.melt(
    id_vars=["Frame"],
    value_vars=[
        "Feret_inv",
        "Feret_inv_filt",
        "Circ",
        "Circ_filt",
        "Perim_inv",
        "Perim_inv_filt",
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
# sns.set(font_scale=2
plt.show()

# %%%% Aplicación al dataset

df = df.groupby(["Batch", "Fenotype", "Fish"], group_keys=False).apply(
    apply_gaussian_filter, column="Circ", new_column_name="Circ_filt", sigma=2
)

df = df.groupby(["Batch", "Fenotype", "Fish"], group_keys=False).apply(
    apply_gaussian_filter, column="Feret_inv", new_column_name="Feret_inv_filt", sigma=3
)

df = df.groupby(["Batch", "Fenotype", "Fish"], group_keys=False).apply(
    apply_gaussian_filter, column="Perim_inv", new_column_name="Perim_inv_filt", sigma=2
)

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

Previamente use la Solidity = area/convex area, pero es una variable muy ruidosa y sensible a irregularidades, voy a en su lugar, usar la inversa del feret (previo filtrado y normalizado min-max. Despues usare un threshold que indicara que el pez esta replegado si su valor es inferior/superior al Threshold. Se pueden usar otras magnitudes, por lo que también voy a realizar lo mismo usando la Circulariry.
"""

# %%% Plot Ejemplo threshold

df_temp = df.loc[("batch 13", "KO44", "ZebraF_1")]
Variable_plot = "Feret_inv_filt"

g = sns.lineplot(data=df_temp, x="Time", y=Variable_plot)
g.axhline(0.008, color="red")
g.set_title("Threshold on " + Variable_plot, size=25)
plt.show()

# %%% [md]
"""
Para no usar el mismo threshold en cada pez, voy a normalizar cada Feret y Circularity entre 0-1, así sí tener un threshold estandar.
"""


# %%% Normalización 0-1
def normalize_group(group, variable):
    # Sort the signal to find the three minimum and three maximum values
    sorted_signal = group[variable].sort_values()
    # Calculate the average of the three minimum values
    three_min_avg = sorted_signal.head(12).mean()
    # Calculate the average of the three maximum values
    three_max_avg = sorted_signal.tail(8).mean()
    # Normalize the signal using the average of the three min and three max values
    group[variable + "_norm"] = (group[variable] - three_min_avg) / (
        three_max_avg - three_min_avg
    )
    # Ensure all values are between 0 and 1
    group[variable + "_norm"] = group[variable + "_norm"].clip(0, 1)
    return group


# Apply normalization 0-1
df = df.groupby(["Batch", "Fenotype", "Fish"], group_keys=False).apply(
    normalize_group, "Feret_inv_filt"
)
df = df.groupby(["Batch", "Fenotype", "Fish"], group_keys=False).apply(
    normalize_group, "Circ_filt"
)
df = df.groupby(["Batch", "Fenotype", "Fish"], group_keys=False).apply(
    normalize_group, "Perim_inv_filt"
)

# %%%% Plot ejemplo
df_temp = df.loc[("batch 11", "WT", "ZebraF_6")]
Variable_plot = "Feret_inv_filt_norm"
# Variable_plot = "Circ_filt_norm"

g = sns.lineplot(data=df_temp, x="Time", y=Variable_plot)
g.axhline(0.6, color="red")
g.set_title("Threshold on " + Variable_plot, size=25)
plt.show()

# %%% Comparación usando un threshold fijo

Variable_plot = "Feret_inv_filt_norm"

threshold = 0.6
time_over_Thr = (
    df.groupby(["Batch", "Fenotype", "Fish"])[Variable_plot]
    .apply(lambda x: (x > threshold).sum())
    .reset_index()
    .rename(columns={Variable_plot: "contracted"})
)

time_over_Thr["contracted_perc"] = 100 * time_over_Thr.contracted / 1550

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
    "Porcentaje del tiempo que pasa replegado (sobre el Threshold = "
    + str(threshold)
    + ")",
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

### Feret_inv_filt_norm
Como he visto que contiene outliers, los he eliminado basandome en la regla del IQR
"""


# %%% Construcción del DF de evolución del resultado con el threshold

Variable_plot = "Feret_inv_filt_norm"

threshold_result = pd.DataFrame(
    columns=["Threshold", "Batch", "Fenotype", "Mean_diff", "CI"]
)
ref = ko44 = ko179 = np.nan

for thr in np.arange(0.2, 1.01, 0.01):  # iteración sobre el threshold
    # data frame con los valores para ese threshold
    time_over_Thr = (
        df.groupby(["Batch", "Fenotype", "Fish"], group_keys=False)[Variable_plot].agg(
            contracted=lambda x: (x > thr).sum(),
            contracted_perc=lambda x: 100 * (x > thr).sum() / len(x),
        )
        # .reset_index()
        # .dropna()
    )

    # generación de los resultados de ese dataframe y añadidos al dataframe de los resultados. La iteración se hace sobre los batch para calcular los valores para cada fenotipo
    for batch, group in time_over_Thr.groupby("Batch", group_keys=False):
        WT = (
            group
            # .dropna()
            # .drop("Batch", axis=1)
            .loc[
                group.index.get_level_values("Fenotype") == "WT"
            ].contracted_perc.dropna()
        )
        grouped_batch = group.dropna().groupby("Fenotype")
        for fenotype, group in grouped_batch:  # loop over each possible fenotype
            if fenotype != "WT":
                mut = group.contracted_perc
                new_row = {
                    "Threshold": thr,
                    "Batch": batch,
                    "Fenotype": fenotype,
                    "Mean_diff": np.mean(remove_outliers_iqr(WT))
                    - np.mean(remove_outliers_iqr(mut)),
                    "Median_diff": np.median(WT) - np.median(mut),
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
g.set_title("Diference of " + Variable_plot + " Batch Mean values along Threshold")

for hue in df_plot.hue.unique():
    g.fill_between(
        x="Threshold",
        y1="CI_up",
        y2="CI_down",
        alpha=0.1,
        data=df_plot[df_plot.hue == hue],
    )
g.axhline(0, color="black", linestyle="--", linewidth=0.5)
plt.show()


# %%%  Circularity vs Threshold Plot [md]
"""
### Circularity plot
Lo mismo para la circularity
"""

# %%%% Construcción del DF y plot

Variable_plot = "Circ_filt_norm"

threshold_result = pd.DataFrame(
    columns=["Threshold", "Batch", "Fenotype", "Mean_diff", "CI"]
)
ref = ko44 = ko179 = np.nan

for thr in np.arange(0.2, 1.01, 0.01):  # iteración sobre el threshold
    # data frame con los valores para ese threshold
    time_over_Thr = (
        df.groupby(["Batch", "Fenotype", "Fish"], group_keys=False)[Variable_plot].agg(
            contracted=lambda x: (x > thr).sum(),
            contracted_perc=lambda x: 100 * (x > thr).sum() / len(x),
        )
        # .reset_index()
        # .dropna()
    )

    # generación de los resultados de ese dataframe y añadidos al dataframe de los resultados. La iteración se hace sobre los batch para calcular los valores para cada fenotipo
    for batch, group in time_over_Thr.groupby("Batch", group_keys=False):
        WT = (
            group
            # .dropna()
            # .drop("Batch", axis=1)
            .loc[
                group.index.get_level_values("Fenotype") == "WT"
            ].contracted_perc.dropna()
        )
        grouped_batch = group.dropna().groupby("Fenotype")
        for fenotype, group in grouped_batch:  # loop over each possible fenotype
            if fenotype != "WT":
                mut = group.contracted_perc
                new_row = {
                    "Threshold": thr,
                    "Batch": batch,
                    "Fenotype": fenotype,
                    "Mean_diff": np.mean(remove_outliers_iqr(WT))
                    - np.mean(remove_outliers_iqr(mut)),
                    "Median_diff": np.median(WT) - np.median(mut),
                    "CI": np.ptp(
                        stats.ttest_ind(WT, mut).confidence_interval(
                            confidence_level=0.80
                        )
                    )
                    / 2,
                }
                threshold_result.loc[len(threshold_result)] = new_row

# el plot
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
g.set_title("Diference of " + Variable_plot + " Batch Mean values along Threshold")

for hue in df_plot.hue.unique():
    g.fill_between(
        x="Threshold",
        y1="CI_up",
        y2="CI_down",
        alpha=0.1,
        data=df_plot[df_plot.hue == hue],
    )

g.axhline(0, color="black", linestyle="--", linewidth=0.5)
plt.show()


# %%%  Perim_inv vs Threshold Plot [md]
"""
### Perimeter inverse plot
Lo mismo para el perimetro, que tenia unos plots claros
"""

# %%%% Construcción del DF y plot

Variable_plot = "Perim_inv_filt_norm"

threshold_result = pd.DataFrame(
    columns=["Threshold", "Batch", "Fenotype", "Mean_diff", "CI"]
)
ref = ko44 = ko179 = np.nan

for thr in np.arange(0.2, 1.01, 0.01):  # iteración sobre el threshold
    # data frame con los valores para ese threshold
    time_over_Thr = (
        df.groupby(["Batch", "Fenotype", "Fish"], group_keys=False)[Variable_plot].agg(
            contracted=lambda x: (x > thr).sum(),
            contracted_perc=lambda x: 100 * (x > thr).sum() / len(x),
        )
        # .reset_index()
        # .dropna()
    )

    # generación de los resultados de ese dataframe y añadidos al dataframe de los resultados. La iteración se hace sobre los batch para calcular los valores para cada fenotipo
    for batch, group in time_over_Thr.groupby("Batch", group_keys=False):
        WT = (
            group
            # .dropna()
            # .drop("Batch", axis=1)
            .loc[
                group.index.get_level_values("Fenotype") == "WT"
            ].contracted_perc.dropna()
        )
        grouped_batch = group.dropna().groupby("Fenotype")
        for fenotype, group in grouped_batch:  # loop over each possible fenotype
            if fenotype != "WT":
                mut = group.contracted_perc
                new_row = {
                    "Threshold": thr,
                    "Batch": batch,
                    "Fenotype": fenotype,
                    "Mean_diff": np.mean(remove_outliers_iqr(WT))
                    - np.mean(remove_outliers_iqr(mut)),
                    "Median_diff": np.median(WT) - np.median(mut),
                    "CI": np.ptp(
                        stats.ttest_ind(WT, mut).confidence_interval(
                            confidence_level=0.80
                        )
                    )
                    / 2,
                }
                threshold_result.loc[len(threshold_result)] = new_row

# el plot
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
g.set_title("Diference of " + Variable_plot + " Batch Mean values along Threshold")

for hue in df_plot.hue.unique():
    g.fill_between(
        x="Threshold",
        y1="CI_up",
        y2="CI_down",
        alpha=0.1,
        data=df_plot[df_plot.hue == hue],
    )

g.axhline(0, color="black", linestyle="--", linewidth=0.5)
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

Variable = "Circ_filt"
# Variable = "Perim_inv_filt"
# Variable = "Feret_inv_filt"

df_temp = df.loc[("batch 7", "KO44", "ZebraF_1"), Variable]

peaks, _ = find_peaks(
    df_temp, height=0.4, prominence=0.03, threshold=0.0, distance=2, width=1
)  # ajustar estos parametros
# he comprobado que algun gusano realiza contracciones y extensiones en 2 frames, por lo que distance = 1 & width = 1

plt.plot(df_temp.values)
plt.plot(peaks, df_temp.iloc[peaks], "2", markersize=24)
plt.title("Picos encontrados sobre la Roundness", size=20)
plt.show()

# %%% [md]
"""
Es interesante ver como funciona sobre todos los gusanos
"""

# %%% Peak finder a todos los Zebra

df["unique_fish"] = [" ".join(tup) for tup in df.index.values]

# i = 0
# for f in sorted(set(df.loc["batch 7"].unique_fish)):
#     i += 1
#     if i > 30:
#         break
#     fish_temp = df[df.unique_fish == f].Circ_filt  # ajustar estos parametros
#     peaks, _ = find_peaks(
#         fish_temp, height=0.4, prominence=0.02, threshold=0.0, distance=2, width=1
#     )

#     plt.plot(fish_temp.values)
#     plt.plot(peaks, fish_temp.iloc[peaks], "2", markersize=24)
#     plt.title(f)
#     plt.show()
# time.sleep(0.5)

# %%% [md]
"""
Cuando visualizo las gráficas de todos los Zebra, todos los picos están bien encontrados.

Aplicar un filtro a la gráfica del movimiento NO es necesario, pues los peaks estan bien encontrados
"""
# %%% [md]
"""
## Número de picos por condición

Represento el número de picos por condición y batch. Dado que todos los videos duran el mismo tiempo, se puede asociar a la frecuencia.

### Circularity

Usando la Circularity filtrada
"""

# %%% Por condición

Variable_plot = "Circ_filt"
peaks_df = (
    df.groupby(["Batch", "Fenotype", "Fish"])[Variable_plot]
    .apply(
        lambda x: len(
            find_peaks(
                x,
                height=0.4,
                prominence=0.02,  # no cambia mucho si pongo 0.08
                threshold=0.0,
                distance=2,
                width=1,  # para magnitudes 0-1
                # x,height=0.005,prominence=0.002,threshold=0.0,distance=2,width=1,  # para perimetro_inv
            )[0]
        )
    )
    .dropna()
    .reset_index()
    .rename(columns={Variable_plot: "N_peaks"})
)

grped_bplot = sns.catplot(
    x="Batch",
    y="N_peaks",
    hue="Fenotype",
    kind="box",
    legend=True,
    showfliers=False,
    height=6,
    aspect=1.9,
    data=peaks_df,
    hue_order=["WT", "KO44", "KO179"],
)

# handles, labels = grped_bplot.get_legend_handles_labels()
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
no cambia hasta 0.8, el cual es ya un valor extremo para Circ. Lo mismo con `prominence`en el
intervalo 0.01-0,6 y es invariante hasta valores extremos a partir de 0.4 (altura del pico)

### Conclusión

El método funciona correctamente, pero hay mucha variabilidad interbatch. Aunque los WT son similares en la mayoria de batches.
Recomiendo repasar las gráficas de los peces comparandolas con las fotos de microscopia y ver que videos contienen defectos. 
"""

# %% FFT [md]
"""
# Fourier Transform & Periodogram

En lugar de contar el número de picos, voy a usar transformar las señales al dominios de las frecuencias. Con esto busco encontrar frecuencias más fuertes para alguna condición. Usaré la transformada de Fourier para encontrar las frecuencias más fuertes en el dominio de frecuencias y lo traspasaremos a periodos. Repetiré el estudio usando el Periodograma de XXXX. Complementaré el estudio tratando de usar una Wavelet Transformation para convertir los picos en una señal más sinosoidal.

## Ejemplo FFT 1 pez
"""

# %%% Zebra Ejemplo FFT - Analisis de Frecuencias del movimiento para un Pez Zebra


zebra_temp = df.loc[("batch 7", "WT", "ZebraF_3")]

signal = zebra_temp.Circ_filt.values  # signal as array
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


def plot_fft_filter(
    signal, sample_rate, time_step, N_points, time_points, f="Zebra", filtro_psd=0.1
):
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
    # zebra_fft_fil[zebra_psd > 0.2] = 0
    # # lower power limit
    zebra_fft_fil[zebra_psd < filtro_psd * max(zebra_psd)] = 0
    # bandpass filter
    zebra_fft_fil[abs(zebra_freqs > 3)] = 0  # frecuencia máxima
    zebra_fft_fil[abs(zebra_freqs < 0.02)] = 0  # frecuencia máxima
    zebra_fil_psd = np.abs(zebra_fft_fil) ** 2 / (sample_rate * N_points)

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(21, 13))
    fft_plot_x = zebra_freqs[
        (zebra_freqs < 3) & (zebra_freqs > 0)
    ]  # set max freq to plot

    sns.lineplot(x=time_points, y=signal, ax=axes[0, 0], color="r")
    sns.lineplot(x=fft_plot_x, y=zebra_psd[0 : len(fft_plot_x)], ax=axes[0, 1])
    sns.lineplot(
        x=fft_plot_x, y=zebra_fil_psd[0 : len(fft_plot_x)], ax=axes[1, 0]
    ).set_xlim(-0.05, 1.1)
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


plot_fft_filter(signal, sample_rate, time_step, N_points, time_points, filtro_psd=0.1)

# %%%%


# df["unique_fish"] = [" ".join(tup) for tup in df.index.values]

# # Filtro para ver por batch
# dfa = df.loc["batch 11"]


# for f in sorted(set(dfa.unique_fish)):
#     signal = dfa[dfa.unique_fish == f].Circ_filt_norm.values  # ajustar estos parametros
#     signal = detrend(signal, axis=0)
#     sample_rate = 9  # frames / s
#     time_step = 1 / sample_rate
#     N_points = len(signal)
#     # t = np.arange(0, N_points/sample_rate, time_step)
#     time_points = zebra_temp.Time
#     plot_fft_filter(signal, sample_rate, time_step, N_points, time_points, f)

# %%%% [md]
"""
Puede que haya un problema con el filtrado pensando en detectar las frecuencias . para señales que tienen pocos picos, la frecuencia se puede calcular mas o menos, pero las interferencias de las señas tratan de ajustar los picos, dando una PSD algo más alta para frecuencias mayores que posiblemente no sean significativas. tengo que estudiar esto un poco más, pero el efecto es que la FFT modela la forma de la señal y no las frecuencias de los coletazos

para un pez que da 1 coletazo la freq es 1/173=0.006Hz seria la fr min
para la maxima, un pez que da1 coletazo en 3 frames = 3Hz seria fr max, aunque arriba la he filtrado para 

"""
# %%%% [md]
"""
No tengo muy claro que este sistema funcione para de verdad detectar las frecuencias a las que da los coletazos. La FT trata de ajustar la forma de la señal y no sus frecuencias, por lo que hay interferencias. Aún así puedo probar el agrupamiento por batch, pero creo que es mejor probar la autocorrelation y el periodograma antes
"""


# %%%% Calculo de la FFT para cada zebra

sample_rate = 9  # frames / s
time_step = 1 / sample_rate  # Delta t
N_points = 1550  # lenght signal


def fft_filter(group, variable, sample_rate, time_step, N_points, filtro_psd=0.1):
    signal = detrend(group[variable].values, axis=0)
    # signal /= np.std(signal)
    # FFT
    zebra_fft = fft.fft(signal, norm="backward")
    # Calculos para frequencias en seg pare representar en el eje X
    zebra_freqs = fft.fftfreq(N_points, time_step)
    # Power spectral density
    zebra_psd = np.abs(zebra_fft) ** 2 / (sample_rate * N_points)
    # Filtrado FFT - Filtro por la psd
    zebra_fft_fil = np.array(zebra_fft)
    # # upper power limit
    # zebra_fft_fil[zebra_psd > 0.2] = 0
    # # lower power limit
    zebra_fft_fil[zebra_psd < filtro_psd * max(zebra_psd)] = 0
    # bandpass filter
    zebra_fft_fil[abs(zebra_freqs > 3)] = 0  # frecuencia máxima
    zebra_fft_fil[abs(zebra_freqs < 0.01)] = 0  # frecuencia máxima
    zebra_fil_psd = np.abs(zebra_fft_fil) ** 2 / (sample_rate * N_points)
    group["PSD"] = zebra_fil_psd  # / max(zebra_fil_psd)
    group["freqs"] = zebra_freqs
    return group


df_fft = df.groupby(["Batch", "Fenotype", "Fish"], group_keys=False).apply(
    fft_filter,
    Variable_plot,
    sample_rate,
    time_step,
    N_points,
    filtro_psd=0.000001,
)

df_fft = (
    df_fft.loc[(df_fft.freqs >= 0.01) & (df_fft.freqs <= 3)]
    .drop("Circ_filt", axis=1)
    .reset_index()
)

# Variable = "Circ_filt"
# df_temp = df2.loc[("batch 7", "KO44", "ZebraF_1")]
# plt.plot(df_temp.PSD.values)
# plt.title("Señal cuadrada modulada", size=20)
# plt.show()

# %%%% Bining de la FFT para simplificarla


def aggregate_n_rows(group, n=4):
    # Aggregate columns differently
    aggregated = group.groupby(group.index // n).agg(
        {
            "PSD": "sum",  # Sum the 'Value1' column
            "freqs": "mean",  # Average the 'Value2' column
        }
    )
    return aggregated


df_bin = (
    df_fft.groupby(["Batch", "Fenotype", "Fish"], as_index=True)
    .apply(lambda x: aggregate_n_rows(x, n=10))
    .dropna()
    .reset_index()
)

# %%%% plot por batch
# for b in df_bin.Batch.sort_values().unique():
#     temp_plot = df_bin.loc[df_bin.Batch == b]

#     sns.lineplot(
#         data=temp_plot,
#         x="freqs",
#         y="PSD",
#         hue="Fenotype",
#         hue_order=["WT", "KO44", "KO179"],
#         estimator="median",
#         errorbar=("ci", 80),
#     ).set_xlim(-0.05, 1.5)
#     plt.title(b, size=20)
#     plt.show()

# %%% FT sobre señal cuadrada [md]
"""
Con objeto de evitar los problemas dados por que la FT encuentra frecuencias que modelan las señales y enmascara las frecuencias de movimientos, voy a generar una señal cuadrada donde cada valor de señal será un coletazo detectado con los picos
"""

# %%%% Señal cuadrada a 1 Zebra

Variable = "Circ_filt"
# Variable = "Perim_inv_filt"
# Variable = "Feret_inv_filt"

df_temp = df.loc[("batch 7", "WT", "ZebraF_3"), Variable]

peaks, _ = find_peaks(
    df_temp, height=0.3, prominence=0.01, threshold=0.0, distance=2, width=1
)

signal2 = np.zeros(len(df_temp))
signal2[peaks] = 1
signal2[peaks + 1] = 1  # para que tengan una longitud los altos de la señal
signal2 = signal2 - np.mean(signal2)

plt.plot(df_temp.values)
plt.plot(signal2, "r", alpha=0.6)
plt.title("Señal cuadrada modulada", size=20)
plt.show()

# %%%% [md]
"""
Ahora se puede hacer la FT sobre esta señal y buscamos las frecuencias
"""

# %%%% FFT sobre la señal cuadrada

plot_fft_filter(signal2, sample_rate, time_step, N_points, time_points, filtro_psd=0.2)

# %%%% [md]
"""
En este caso tenemos unas frecuencias, vamos a calcularlo para cada zebra y agregar por condición
"""
# %%%% Construcción de DF con señales cuadradas

Variable_plot = "Circ_filt"


def señal_cuadrada(group):
    peaks, _ = find_peaks(
        group, height=0.3, prominence=0.01, threshold=0.0, distance=2, width=1
    )
    signal2 = np.zeros(len(group))
    signal2[peaks] = 1
    signal2[peaks + 1] = 1
    # signal2 = signal2 - np.mean(signal2)
    return signal2


df2 = (
    df.groupby(["Batch", "Fenotype", "Fish"], as_index=True)[Variable_plot]
    .apply(señal_cuadrada)
    .dropna()
    .reset_index()
)
df2 = df2.explode(Variable_plot, ignore_index=True)
df2.set_index(["Batch", "Fenotype", "Fish"], inplace=True)

# %%%% Calculo de la FFT para cada zebra

sample_rate = 9  # frames / s
time_step = 1 / sample_rate  # Delta t
N_points = 1550  # lenght signal


def fft_filter(group, variable, sample_rate, time_step, N_points, filtro_psd=0.1):
    signal = detrend(group[variable].values, axis=0)
    signal /= np.std(signal)
    # FFT
    zebra_fft = fft.fft(signal, norm="backward")
    # Calculos para frequencias en seg pare representar en el eje X
    zebra_freqs = fft.fftfreq(N_points, time_step)
    # Power spectral density
    zebra_psd = np.abs(zebra_fft) ** 2 / (sample_rate * N_points)
    # Filtrado FFT - Filtro por la psd
    zebra_fft_fil = np.array(zebra_fft)
    # # upper power limit
    # zebra_fft_fil[zebra_psd > 0.2] = 0
    # # lower power limit
    zebra_fft_fil[zebra_psd < filtro_psd * max(zebra_psd)] = 0
    # bandpass filter
    zebra_fft_fil[abs(zebra_freqs > 3)] = 0  # frecuencia máxima
    zebra_fft_fil[abs(zebra_freqs < 0.01)] = 0  # frecuencia máxima
    zebra_fil_psd = np.abs(zebra_fft_fil) ** 2 / (sample_rate * N_points)
    group["PSD"] = zebra_fil_psd  # / max(zebra_fil_psd)
    group["freqs"] = zebra_freqs
    return group


df2_fft = df2.groupby(["Batch", "Fenotype", "Fish"], group_keys=False).apply(
    fft_filter,
    Variable_plot,
    sample_rate,
    time_step,
    N_points,
    filtro_psd=0.1,
)

df2_fft = (
    df2_fft.loc[(df2_fft.freqs >= 0.01) & (df2_fft.freqs <= 3)]
    # .drop("Circ_filt", axis=1)
    .reset_index()
)

# Variable = "Circ_filt"
# df_temp = df2.loc[("batch 7", "KO44", "ZebraF_1")]
# plt.plot(df_temp.PSD.values)
# plt.title("Señal cuadrada modulada", size=20)
# plt.show()

# %%%% Bining de la FFT para simplificarla


def aggregate_n_rows(group, n=4):
    # Aggregate columns differently
    aggregated = group.groupby(group.index // n).agg(
        {
            "PSD": "sum",  # Sum the 'Value1' column
            "freqs": "mean",  # Average the 'Value2' column
        }
    )
    return aggregated


df2_bin = (
    df2_fft.groupby(["Batch", "Fenotype", "Fish"], as_index=True)
    .apply(lambda x: aggregate_n_rows(x, n=15))
    .dropna()
    .reset_index()
)


# %%%% plot por batch
# for b in df2_bin.Batch.sort_values().unique():
#     temp_plot = df2_bin.loc[df2_bin.Batch == b]

#     sns.lineplot(
#         data=temp_plot,
#         x="freqs",
#         y="PSD",
#         hue="Fenotype",
#         hue_order=["WT", "KO44", "KO179"],
#         estimator="median",
#         errorbar=("ci", 80),
#     )
#     plt.title(b, size=20)
#     plt.show()

# %%% [md]
""" 
No veo que haya demasiada información en la FT, ni siquiera tengo claro que sea correcto aplicarla, ya que para la señal sin manipular, la FT tiende a 'modelarla' mas que a encontrar sus frecuencias basicas, y la señal de picos tiene angulos, lo que hace que las altas frecuencias se disparen. Como alternativas voy a aplicar el periodograma y el Lomb-Scarle, que pienso que pueden reflejar mejor las frecuencias
"""
# %% Periodograma
"""
Según el libro de FFT analisis este método es más adecuado para señales random, así que voy a aplicarlo
"""

# %%% selección de ventana

Variable = "Circ_filt"

df_temp = df.loc[("batch 7", "WT", "ZebraF_2"), Variable]

ventanas = ["boxcar", "blackman", "hann", "flattop", "nuttall"]
for ventana in ventanas:
    x_period, y_period = periodogram(
        df_temp.values,
        detrend="constant",
        fs=sample_rate,
        return_onesided=True,
        axis=0,
        window=ventana,
        scaling="spectrum",
    )

    plt.plot(x_period, y_period)
    plt.xlabel("Frecuency (Hz)")
    plt.title("Función periodograma (usa FFT) y una ventana " + ventana)
    plt.show()

# %%%[md]
"""
Parece que según la ventana puede funcionar mejor que la FFT, lo voy a probar
"""

# %%% Periodograma a todos los cebra

sample_rate = 9  # frames / s
N_points = 1550  # lenght signal


def periodograma(group, variable, sample_rate, filtro_psd=0.1):
    x_period, y_period = periodogram(
        group[variable],
        detrend="constant",
        fs=sample_rate,
        return_onesided=False,
        axis=0,
        window="flattop",
        scaling="spectrum",
    )
    # power filter
    y_period[y_period < filtro_psd * max(y_period)] = 0
    # bandpass filter
    # y_period[abs(x_period > 3)] = 0  # frecuencia máxima
    # y_period[abs(x_period < 0.01)] = 0  # frecuencia máxima
    group[("PSD_" + variable)] = y_period  # / max(zebra_fil_psd)
    group["freqs"] = x_period
    return group


Variable = "Circ_filt"

df_periodogram = df.groupby(["Batch", "Fenotype", "Fish"], group_keys=False).apply(
    periodograma,
    Variable,
    sample_rate,
    filtro_psd=0.000001,
)

df_periodogram = df_periodogram.loc[
    (df_periodogram.freqs >= 0.01) & (df_periodogram.freqs <= 3),
    ("PSD_" + Variable, "freqs"),
].reset_index()

# %%%% Bining de la FFT para simplificarla


def aggregate_n_rows(group, n=4):
    # Aggregate columns differently
    aggregated = group.groupby(group.index // n).agg(
        {
            "PSD_" + Variable: "sum",  # Sum the 'Value1' column
            "freqs": "mean",  # Average the 'Value2' column
        }
    )
    return aggregated


df_periodogram_bin = (
    df_periodogram.groupby(["Batch", "Fenotype", "Fish"], as_index=True)
    .apply(lambda x: aggregate_n_rows(x, n=8))
    .dropna()
    .reset_index()
)


# %%%% plot por batch
for b in df_periodogram.Batch.sort_values().unique():
    temp_plot = df_periodogram.loc[df_periodogram.Batch == b]

    sns.lineplot(
        data=temp_plot,
        x="freqs",
        y="PSD_" + Variable,
        hue="Fenotype",
        hue_order=["WT", "KO44", "KO179"],
        estimator="median",
        errorbar=("ci", 80),
    )
    plt.title(b, size=20)
    plt.show()

# %%%% [md]
"""
Del mismo modo que la FFT, el resultado es ruidoso y poco concluyente, a pesar de que en este caso si se aprecian más altas frecuencias para los WT en algunos batches. Voy a probar con el Lomb-Scarle, ya que pienso que su ajuste por frecuencias es el más adecuado para encontrar frecuencias en las señales que tenemos, debido a que ajusta por minimos cuadrados señales periodicas y no mediante el paso al dominio de las frecuencias
"""
# %% Lomb-scarle Periodogram

Variable = "Circ_filt"

df_temp = df.loc[("batch 7", "WT", "ZebraF_2"), ("Time", Variable)]

periods = np.linspace(0.01, 3, 100)
T_periodogram = 2 * np.pi / periods  # Periodos en segundos
y_periodogram = lombscargle(
    df_temp.Time.values,
    detrend(df_temp[Variable].values, axis=0),
    periods,
    normalize=True,
)

plt.plot(periods, y_periodogram)
plt.xlabel("Period (Seg)")
plt.title("Función lombscargle -> ajusta por LSE")
plt.show()

# %%% LS to every Zebra


def LS(group, Variable, N_freqs=100, filtro_psd=0.1):
    freqs = np.linspace(0.01, 3, N_freqs)
    y_periodogram = lombscargle(
        group.Time.values,
        detrend(group[Variable].values, axis=0),
        freqs,
        normalize=True,
    )
    # power filter
    y_periodogram[y_periodogram < filtro_psd * max(y_periodogram)] = 0
    result_df = pd.DataFrame({"freqs": freqs, "y_periodogram": y_periodogram})
    return result_df


Variable = "Circ_filt"

df_LS = (
    df.groupby(["Batch", "Fenotype", "Fish"], group_keys=True)
    .apply(
        LS,
        Variable,
        N_freqs=80,
        filtro_psd=0.2,
    )
    .reset_index()
)

# %%% Plot every condition
for b in df_LS.Batch.sort_values().unique():
    temp_plot = df_LS.loc[df_LS.Batch == b]

    sns.lineplot(
        data=temp_plot,
        x="freqs",
        y="y_periodogram",
        hue="Fenotype",
        hue_order=["WT", "KO44", "KO179"],
        estimator="median",
        errorbar=("ci", 80),
    )
    plt.title(b, size=20)
    plt.show()

# %%% [md]
"""
No ha funcionado para nada, posiblemente por la cantidad de ruido que hay y que estamos ajustando una señal sinosoida. Quizas tendria más sentido una señal cuadrada. Creo que el error esta en que los picos no se ajustan bien. Se podría arreglar usando una transformada wavelet o mirando la ACF o mirando la ACF de la señal cuadrada. Pero de momento lo aparco
"""

# %% Autocorrelation

Variable = "Circ_filt"

df_temp = df.loc[("batch 13", "WT", "ZebraF_3"), ("Time", Variable)]

from statsmodels.graphics.tsaplots import plot_acf

plot_acf(df_temp.Circ_filt, lags=200)  # data: your time series
# lags: number of 'periods' you will like to investigate

plt.show()

sns.lineplot(data=df_temp, x="Time", y="Circ_filt")
plt.show()

# %%% [md]
"""
No tengo muy claro que este sistema vaya a funcionar tampoco, ya que esta costando sacar las frecuencias de coletazos del pez. Pero se podría probar la ACF agregada o filtrando los lags altos. tambien la ACF de la señal cuadrada. Pero estoy un poco harto y lo voy a dejar aparcado a un futuro.

"""
# %% Estadistica sobre tiempo promedio entre coletazos [md]
"""
Como ultimo recurso, voy a calcular el tiempo promedio entre coletazos y con el realizar una estadistica con las distribuciones exponenciales, gamma o poisson
"""
