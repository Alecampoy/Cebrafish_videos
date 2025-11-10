# -*- coding: utf-8 -*-

# %% Intro [md]
"""
# Tracking Experiment
# **Análisis del Tracking del Movimiento de Zebrafish en una placa**
### Author: Alejandro Campoy Lopez  
"""

# %% Librerias
from IPython import get_ipython

# limpia variables y librerias antiguas
get_ipython().magic("reset -sf")

import warnings, math
import pandas as pd
import re
import os, platform
import glob
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as st
from scipy import stats
from statannotations.Annotator import Annotator



def get_files_in_folder(folder_path):
    file_list = glob.glob(folder_path + "/**", recursive=True)
    files = [file for file in file_list if not os.path.isdir(file)]
    return files


plt.rcParams["figure.figsize"] = (15, 8)

# %matplotlib inline

warnings.filterwarnings("ignore")

# %% Lectura Archivos [md]
"""
El archivo WT2 del batch 6 esta roto, cogerlo del disco duro y otros mas, recopiarlos del disco duro
# Lectura de Archivos
Lectura de todos los archivos csv con los resultados de los diferentes batches.
Se añade una columna representando el gusano y el batch mediante el uso de regex
"""

# %%% Load batches 

if platform.system() == "Windows":
    folder_path = "P:\CABD\Lab Ozren\datos csv\Experimento Tracking\Batches Laura\\No Strict"
    #folder_path = "P:\CABD\Lab Ozren\datos csv\Experimento Tracking\\test"
else:
    folder_path = "/home/ale/pCloudDrive/CABD/Lab Ozren/datos csv\Experimento Tracking\Batches Laura"

files = get_files_in_folder(folder_path)

df = []
for f in files:
    if ".csv" in f:
        csv = pd.read_csv(f, sep=";")
        csv.insert(0, "Batch", re.search("batch \d+", f.lower()).group(0))
        csv.insert(1, "Fenotype", re.search("(KO\d*|WT)", f.upper()).group(1))
        csv.insert(2, "Fish", "ZebraF_" + re.search("(\d+)(.tif)", f.lower()).group(1))
        df.append(csv)
        del (csv, f)

df = pd.concat(df)
df.Time = (df.Frame - 1) / 6  # 6 fps
# df = df.drop(df.columns[-1], axis=1)

elements = round(
    pd.crosstab(index=df.Batch, columns=df.Fenotype) / 5601
)  # divided by lengh of the video
print(str(elements).replace(".0", "").replace("],", "]\n"))


df["Batch"] = pd.Categorical(
    df["Batch"],
    categories=[
        "batch 1",
        "batch 2",
        "batch 3",
        "batch 4",
        "batch 5",
        "batch 6",
        "batch 7",
        "batch 8",
    ],
    ordered=True,
)

# # %% solo para el test

# df["Batch"] = pd.Categorical(
#     df["Batch"],
#     categories=[
#         "batch 3",
#         "batch 44",
#         "batch 70",
#     ],
#     ordered=True,
# )

# %% Elimino peces muertos del NO STRICT. estos videos pueden estar bien en el otro analisis y substiirse

pez_video_mal = (
    (df['Batch'] == 'batch 1')
    & (df['Fenotype'] == "KO179")
    & (df['Fish'] == 'ZebraF_12'))
df = df[~pez_video_mal]

pez_video_mal = (
    (df['Batch'] == 'batch 1')
    & (df['Fenotype'] == "WT")
    & (df['Fish'] == 'ZebraF_3'))
df = df[~pez_video_mal]

pez_video_mal = (
    (df['Batch'] == 'batch 2')
    & (df['Fenotype'] == "WT")
    & (df['Fish'] == 'ZebraF_9'))
df = df[~pez_video_mal]

pez_video_mal = (
    (df['Batch'] == 'batch 2')
    & (df['Fenotype'] == "WT")
    & (df['Fish'] == 'ZebraF_11'))
df = df[~pez_video_mal]

pez_video_mal = (
    (df['Batch'] == 'batch 3')
    & (df['Fenotype'] == "KO44")
    & (df['Fish'] == 'ZebraF_9'))
df = df[~pez_video_mal]

pez_video_mal = (
    (df['Batch'] == 'batch 3')
    & (df['Fenotype'] == "KO44")
    & (df['Fish'] == 'ZebraF_10'))
df = df[~pez_video_mal]

pez_video_mal = (
    (df['Batch'] == 'batch 3')
    & (df['Fenotype'] == "KO179")
    & (df['Fish'] == 'ZebraF_7'))
df = df[~pez_video_mal]

pez_video_mal = (
    (df['Batch'] == 'batch 3')
    & (df['Fenotype'] == "WT")
    & (df['Fish'] == 'ZebraF_3'))
df = df[~pez_video_mal]

pez_video_mal = (
    (df['Batch'] == 'batch 4')
    & (df['Fenotype'] == "KO179")
    & (df['Fish'] == 'ZebraF_10'))
df = df[~pez_video_mal]

pez_video_mal = (
    (df['Batch'] == 'batch 4')
    & (df['Fenotype'] == "KO185")
    & (df['Fish'] == 'ZebraF_2'))
df = df[~pez_video_mal]

pez_video_mal = (
    (df['Batch'] == 'batch 5')
    & (df['Fenotype'] == "KO185")
    & (df['Fish'] == 'ZebraF_8'))
df = df[~pez_video_mal]

pez_video_mal = (
    (df['Batch'] == 'batch 5')
    & (df['Fenotype'] == "WT")
    & (df['Fish'] == 'ZebraF_12'))
df = df[~pez_video_mal]

pez_video_mal = (
    (df['Batch'] == 'batch 6')
    & (df['Fenotype'] == "KO179")
    & (df['Fish'] == 'ZebraF_5'))
df = df[~pez_video_mal]

pez_video_mal = (
    (df['Batch'] == 'batch 6')
    & (df['Fenotype'] == "KO179")
    & (df['Fish'] == 'ZebraF_11'))
df = df[~pez_video_mal]

pez_video_mal = (
    (df['Batch'] == 'batch 6')
    & (df['Fenotype'] == "KO185")
    & (df['Fish'] == 'ZebraF_5'))
df = df[~pez_video_mal]

pez_video_mal = (
    (df['Batch'] == 'batch 6')
    & (df['Fenotype'] == "KO185")
    & (df['Fish'] == 'ZebraF_10'))
df = df[~pez_video_mal]

pez_video_mal = (
    (df['Batch'] == 'batch 6')
    & (df['Fenotype'] == "KO185")
    & (df['Fish'] == 'ZebraF_11'))
df = df[~pez_video_mal]

pez_video_mal = (
    (df['Batch'] == 'batch 6')
    & (df['Fenotype'] == "KO185")
    & (df['Fish'] == 'ZebraF_12'))
df = df[~pez_video_mal]

pez_video_mal = (
    (df['Batch'] == 'batch 6')
    & (df['Fenotype'] == "WT")
    & (df['Fish'] == 'ZebraF_1'))
df = df[~pez_video_mal]

pez_video_mal = (
    (df['Batch'] == 'batch 6')
    & (df['Fenotype'] == "WT")
    & (df['Fish'] == 'ZebraF_2'))
df = df[~pez_video_mal]

pez_video_mal = (
    (df['Batch'] == 'batch 7')
    & (df['Fenotype'] == "KO179")
    & (df['Fish'] == 'ZebraF_5'))
df = df[~pez_video_mal]

pez_video_mal = (
    (df['Batch'] == 'batch 7')
    & (df['Fenotype'] == "KO179")
    & (df['Fish'] == 'ZebraF_9'))
df = df[~pez_video_mal]

pez_video_mal = (
    (df['Batch'] == 'batch 7')
    & (df['Fenotype'] == "KO185")
    & (df['Fish'] == 'ZebraF_1'))
df = df[~pez_video_mal]

pez_video_mal = (
    (df['Batch'] == 'batch 7')
    & (df['Fenotype'] == "KO185")
    & (df['Fish'] == 'ZebraF_12'))
df = df[~pez_video_mal]

pez_video_mal = (
    (df['Batch'] == 'batch 7')
    & (df['Fenotype'] == "WT")
    & (df['Fish'] == 'ZebraF_1'))
df = df[~pez_video_mal]

pez_video_mal = (
    (df['Batch'] == 'batch 7')
    & (df['Fenotype'] == "WT")
    & (df['Fish'] == 'ZebraF_2'))
df = df[~pez_video_mal]

pez_video_mal = (
    (df['Batch'] == 'batch 8')
    & (df['Fenotype'] == "KO44")
    & (df['Fish'] == 'ZebraF_3'))
df = df[~pez_video_mal]

pez_video_mal = (
    (df['Batch'] == 'batch 8')
    & (df['Fenotype'] == "WT")
    & (df['Fish'] == 'ZebraF_9'))
df = df[~pez_video_mal]

elements = round(
    pd.crosstab(index=df.Batch, columns=df.Fenotype) / 5601
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

# Loop through all axes in the FacetGrid
for ax in NAs_barplot.axes.flatten():
    # Force update of tick labels from the data
    ax.set_xticks(range(len(NAs["Fish"].unique())))
    ax.set_xticklabels(NAs["Fish"].unique(), rotation=45, ha="right")
    ax.tick_params(axis="x", which="major", pad=10)
    ax.set_xlabel("Fish", fontsize=14)

plt.tight_layout()
plt.show()

# %%% NA Impute


# Ayuda codigo chatgpt
# Define a function to interpolate within each group
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
Este análisis se ha realizado usando strict. Hay NAs pero nigun pez tiene demasiados si consideramos que hemos medido miles de frames. Se han imputado
"""

# %% Calculo de la distancia - Variables auxiliares
"""
Con el dataset limpio genero unas variables auxiliares. La más importante es la distancia al borde (o al centro) normalizada a 1
"""

# Para evaluar la distancia al borde de 0 (border) a 1 (center)
df["Dist_border_dmap"] = (
    df["Mean-Distance"] / 172
)  # valor del radio del pocillo medido de las imagenes

df["Dist_center_dmap"] = abs(df["Dist_border_dmap"] - 1)

# Calculo la distancia y la normalizo con respecto al radio del pocillo -calculado con el Feret-
df["Dist_center_feret"] = np.sqrt(
    (df["X"] - df["Pocillo_XM"]) ** 2 + (df["Y"] - df["Pocillo_YM"]) ** 2
) / (df["Pocillo_diam"] / 2)

df["Dist_border_feret"] = abs(df["Dist_center_feret"] - 1)

# Ocurre que debido a los margenes de error, hay algunos pocillos que tienen el centro con un valor algo mayor que 1. Voy a imputar esto.
# Tambien, para evitar problemas con las distribuciones, voy a quitar los valores límite 0 y 1
df.loc[df.Dist_border_dmap >= 1.000, "Dist_border_dmap"] = 0.9999
df.loc[df.Dist_border_dmap <= 0, "Dist_border_dmap"] = 0.0001
df.loc[df.Dist_center_feret >= 1.000, "Dist_center_feret"] = 0.9999
df.loc[df.Dist_center_feret <= 0, "Dist_center_feret"] = 0.0001

# %% Filtrado de datos debido a detección de otras particulas[md]
"""
# Filtrado de los datos
Voy a eliminar los frames en los que se ha detectado un salto de posición demasiado alto ya que posiblemente se debera al haber detectado una posición anomala en ese frame. P.ejem el pez de ejemplo tiene uno de estos eventos
"""
# %%% Histograma de las distancias rocorridas
# Calculo las distancias en cada paso
df["X_diff"] = df.groupby(["Batch", "Fenotype", "Fish"]).X.diff()
df["Y_diff"] = df.groupby(["Batch", "Fenotype", "Fish"]).Y.diff()
df["dist"] = np.sqrt((df.X_diff**2) + (df.Y_diff**2))

# ejemplo de la distribución de los saltos entre frames
df_temp = df.loc[("batch 8", "KO44", "ZebraF_6")].reset_index()

g = sns.histplot(data=df_temp, x="dist", stat="count", binwidth=9)
g.set_title("Distribution of frame movement of a single zebra")
# g.set_yscale('log')
plt.yscale("log")
plt.show()
# %%% [md]
"""
Se aprecia como hay movimientos en los que la distancia recorrida es muy alta. Esto solo puede deberse a un error en la detección del gusano u otra particula. Considero 180px es el límite de salto por frame, aunque bien podrían ser menos.
Se eliminan e imputan.
"""

# %%% Imputación anomalias

"""
Como se recalculan las distancias entre 2 frames, es posible que una posición outlier aparezca en varios frames, pero solo se elimina una vez, por lo que al recalcular las distancias apareceran de nuevo estos saltos. Para corregirlo se realiza el proceso de filtrado varias veces.
"""

# Para ambos DF
for i in range(15):
    # Imputación en las columnas que tienen medidas
    df.loc[(df.dist > 180), ("X", "Y", "Mean-Distance")] = np.nan
    # imputación por interpolación de los cercanos
    df[["X", "Y", "Mean-Distance"]] = df[["X", "Y", "Mean-Distance"]].interpolate(
        method="linear"
    )

    # Recalculo las distancias
    df["X_diff"] = df.groupby(["Batch", "Fenotype", "Fish"]).X.diff()
    df["Y_diff"] = df.groupby(["Batch", "Fenotype", "Fish"]).Y.diff()
    df["dist"] = np.sqrt((df.X_diff**2) + (df.Y_diff**2))


df_temp = df.loc[("batch 8", "KO44", "ZebraF_6")].reset_index()

g = sns.histplot(data=df_temp, x="dist", stat="density", binwidth=9)
g.set_title("Distribution of frame movement of a single zebra")
g.set_yscale("log")
plt.show()

"""
Se aprecia en la gráfica como se han eliminado los saltos más allá del valor 180
"""

# %%% Recalculo las distancias después de la imputación

df["Dist_border_dmap"] = df["Mean-Distance"] / 170

df["Dist_center_dmap"] = abs(df["Dist_border_dmap"] - 1)

df["Dist_center_feret"] = np.sqrt(
    (df["X"] - df["Pocillo_XM"]) ** 2 + (df["Y"] - df["Pocillo_YM"]) ** 2
) / (df["Pocillo_diam"] / 2)

df.loc[df.Dist_border_dmap >= 1.000, "Dist_border_dmap"] = 0.9999
df.loc[df.Dist_border_dmap <= 0, "Dist_border_dmap"] = 0.0001
df.loc[df.Dist_center_feret >= 1.000, "Dist_center_feret"] = 0.9999
df.loc[df.Dist_center_feret <= 0, "Dist_center_feret"] = 0.0001

# %% Distancia Recorrida [md]
"""
## Distancia Recorrida
Se calcula la distancia que recorre el pez a lo largo del video y se gráfica por batch
"""


# %%% Calculo de la distancia recorrida

# dataframe con la distancia recorrida por el  gusano
Distancia_recorrida = (
    df.dropna()
    .groupby(["Batch", "Fenotype", "Fish"], observed=True, as_index=True)
    .dist.sum(min_count=100)
    .round()
).reset_index()


Distancia_recorrida = Distancia_recorrida.rename(columns={'dist': 'Distancia_recorrida'})


# %%% Box-plot

batches_2_plot = ["batch 1","batch 2","batch 3", "batch 4","batch 5","batch 6","batch 7","batch 8"]
#batches_2_plot = ["batch 1","batch 2","batch 3", "batch 4"]
#batches_2_plot = ["batch 5","batch 6","batch 7","batch 8"]
#batches_2_plot = ["batch 4","batch 5","batch 6","batch 7","batch 8"]
df_plot = Distancia_recorrida.loc[Distancia_recorrida.Batches.isin(batches_2_plot)]
df_plot["Batch"] = df_plot["Batch"].cat.remove_unused_categories()

grped_bplot = sns.catplot(
    x="Batch",
    y="Distancia_recorrida",
    hue="Fenotype",
    kind="box",
    showfliers=False,
    height=6,
    aspect=1.9,
    data=df_plot,
    hue_order=["WT", "KO44", "KO179", "KO185"],
)
# make grouped stripplot
grped_bplot = sns.stripplot(
    x="Batch",
    y="Distancia_recorrida",
    hue="Fenotype",
    jitter=0.18,
    dodge=True,
    legend=False,
    marker="o",
    color="black",
    # palette="Set2",
    data=df_plot,
    hue_order=["WT", "KO44", "KO179", "KO185"],
)

handles, labels = grped_bplot.get_legend_handles_labels()
# grped_bplot.axes.legend(handles[0:3], labels[0:3], title="Fenotype", loc="upper right")

# Set title for the plot
grped_bplot.axes.set_title("Distancia Total Recorrida por el Zebrafish (px)")

plt.show()

# %% Posición del Pez en el pocillo a lo largo del video [md]
""" 
# Posición del Zebra en el Pocillo
A lo largo del video, el pez se posiciona en algún lugar de la placa. Se piensa que la posición que mantiene el Zebra es comportamental, por lo que vamos a estudiar que posición mantiene con respecto al borde d=0 hasta el centro d=1. (Dado que existe simetria radial, solo usamos la distancia respecto al borde)
"""

# %%% Histograma 1 Zebra [md]
"""
## Histograma de 1 solo Pez
Vamos a evaluar el histograma de 1 Zebra. Este nos va a indicar donde se posiciona el Zebra a lo largo del tiempo del video.
"""
# %%% Grafico Histograma 1 Zebra
df_temp = df.loc[("batch 6", "WT", "ZebraF_5")].reset_index()

nbins = 12
g = sns.histplot(
    data=df_temp, x="Dist_border_feret", stat="density", binrange=[0, 1], bins=nbins
)
g.set_title("Distribution of radial position relative to edge of a single zebra")

plt.show()

# %%%% [md]
"""
Se observa como el Zebra se posiciona a lo largo del video. Al estar normalizado a densidad, el tiempo total del video es 1

`density: normalize such that the total area of the histogram equals 1`
"""

# %%% Distribución Radial[md]

"""
## Distribución Radial
Dado que el pocillo es circular y estamos viendo la distribución de su posición radial, hay que compensar la diferencia del area sobre el que el gusano puede distribuirse para una misma posición radial. Se puede visualizar como la diferencia de area que ocupan los anillos para la posición radial. El código de ejemplo se encuentra en `weight area.py`

"""

# %%% Histograma con pesos

# metodo pesos 1
weights_r = 1 / (np.pi * (1 - df_temp.Dist_border_feret.values))
# dada la naturaleza radial de los datos y que la distribución va de 0 siendo el anillo más grande al anillo de menos area 1. Este método no funciona apropiadamente, ya que los valores limite (1, en el centro) generan un peso infinito (o que tiende a) y da problemas en la representación del histograma. Aunque esto se ha corregido imputando para que la distribución este contenida en [0.001, 0.999]

# metodo pesos 2
nbins = 12
bins = np.arange(0, 1 + (1 / (nbins)), 1 / (nbins))
weights_a = (np.pi * np.arange(0, 1 + (1 / (nbins)), 1 / (nbins)) ** 2)[::-1]
# Diferencias del area de las circunferencias, este es el peso
weights_a = -np.diff(weights_a)
# bin al que pertenece cada observación, resto 1 ya que el limite es el valor de la derecha del bin, luego todos los valores mayores a 0 caen en el bin 1 y su peso corresponde al [0]. Por este motivo he eliminado los valores 0 y 1 y sustituido por 0.0001 y 0.9999
bin_of_dist = np.searchsorted(bins, df_temp.Dist_border_feret) - 1
# peso que le corresponde a cada observación
weights_a_ind = 1 / weights_a[bin_of_dist]
# este método de pesos le da el mismo peso a cada observación según en que bin del histograma caiga, de este modo no hay problema en las observaciones de valor límite. Creo que es para hacer la inversa del area del anillo

g = sns.histplot(
    data=df_temp,
    x="Dist_border_feret",
    stat="density",
    weights=weights_a_ind,
    binrange=[0, 1],
    bins=nbins,
)
g.set_title(
    "Weighted Distribution of radial position relative to edge of a single zebra"
)

plt.show()

# %%%% [md]
"""
Un par de ejemplos significativos de la diferencia entre el histograma con y sin ponderar, que puede verse bien en las imagenes son los zebra 4,5, y 6 del batch 7

Hay que definir que pesos se usan. el problema con los pesos individuales `weights_r`es que para los valores extremos calcula un peso muy alto, por lo que creo que es mejor generar unos pesos para los bins y usar el mismo peso para todas las observaciones que caigan dentro del bin.
"""

# %%% Hist por Condición [md]
""" 
## Histograma por Condición 
Voy a ver si, en media, un fenotipo cambia su modo de distribuirse en el pocillo acumulando los histogramas. Esto no es un problema ya que los histogramas estan normalizados con density. 

"""

# %%%% Hitograma por batch sin ponderar

nbins=12
i=0
for j in ["KO44", "KO179", "KO185"]:
    g = sns.FacetGrid(
        data=df.reset_index(),
        row="Batch",
        hue="Fenotype",
        hue_order=["WT", j],
        palette= ["blue", ["red", "green", "yellow"][i]],#"pastel",
        sharex="col",
        sharey=False,
        height=5,
        aspect=4,
    )
    
    # g.fig.suptitle("Evolución temporal de todas las variables para un pez de ejemplo",
    # fontsize=24, fontdict={"weight": "bold"})
    
    g.map_dataframe(
        sns.histplot,
        x="Dist_border_feret",
        element="step",
        edgecolor="black",
        alpha=0.4,
        binrange=[0, 1],
        # cumulative=True,
        bins=nbins,
        stat="density",
        common_norm=False,
        kde=False,
        kde_kws={"bw_adjust": 1},
    )
    g.add_legend()
    g.set_axis_labels(fontsize=20)
    g.fig.suptitle("Accumulated distribution of radial position relative to edge")
    plt.subplots_adjust(top=0.95)
    plt.show()
    i=i+1

# %%%% Calculo Pesos para Ponderar por la distribución radial por Batch

# metodo pesos 1
weights_r = 1 / (np.pi * (1 - df.Dist_border_feret.values))
df["weights_r"] = weights_r

# metodo pesos 1
bins = np.arange(0, 1 + (1 / (nbins)), 1 / (nbins))
weights_a = (np.pi * np.arange(0, 1 + (1 / (nbins)), 1 / (nbins)) ** 2)[::-1]
weights_a = -np.diff(weights_a)
bin_of_dist = np.searchsorted(bins, df.Dist_border_feret) - 1
df["weights_a_ind"] = 1 / weights_a[bin_of_dist]


# %%%% Hist. ponderado por distribución radial

nbins=12
i=0
for j in ["KO44", "KO179", "KO185"]:
    g = sns.FacetGrid(
        data=df.reset_index(),
        row="Batch",
        hue="Fenotype",
        hue_order=["WT", j],
        palette= ["blue", ["red", "green", "yellow"][i]],#"pastel",
        sharex="col",
        sharey=False,
        height=5,
        aspect=4,
    )
    
    # g.fig.suptitle("Evolución temporal de todas las variables para un pez de ejemplo",
    # fontsize=24, fontdict={"weight": "bold"})
    
    g.map_dataframe(
        sns.histplot,
        x="Dist_border_feret",
        element="step",
        edgecolor="black",
        alpha=0.1,
        binrange=[0, 1],
        # cumulative=True,
        bins=nbins,
        stat="density",
        weights="weights_a_ind",
        common_norm=False,
        kde=False,
        kde_kws={"bw_adjust": 0.8},
    )
    g.add_legend()
    g.set_axis_labels(fontsize=20)
    g.fig.suptitle("Accumulated distribution of radial position relative to edge")
    plt.subplots_adjust(top=0.95)
    plt.show()
    i=i+1

# %%%% [md]
"""
Este gráfico muestra la densidad de probabilidad para cada condición y batch, normalizada para cada condición. He agregado los Zebra ya que como cada uno dura lo mismo, puede hacerse ya que todos tendrán el mismo peso en el caso en el que no están ponderados. Debería estar normalizado por la función 'sns.histplot' 

El problema de la distribución ponderada por ser radial viene, creo, del hecho de que *se están agregando todos los peces conjuntamente*. Un pez anomalo va a tener un peso muy fuerte en este gráfico.

*Se hace necesario ver los peces independientemente*

En los casos en los que la distribución de una condición es significativamente diferente a las otras habria que estudiar que la contribución individual de cada Zebra sea razonablemente similar, y no que un solo Zebra sea el que produce la desviacion del histograma o PDF. Lo vemos
"""

# %% Hist apilado ponderado por Zebra [md]

"""
Resulta un plot bastante sucio pero se aprecia la diferencia. Creo que una mejor alternativa será representar unicamente la condición y el batch de interes para ver la intra-distribución de los Zebra y ver que son razonablemente homogeneos y no se debe a un Zebra Oulier.   

Haciendolo de manera ponderada. *Recomiendo verlo para cada batch*
"""
# %%% Hist apilado ponderado por Zebra
"""OJO que me he dado cuenta que la funcion plot considera la density total, y creo que es diferente que la density de la de un zebra en concreto, y por eso tengo seguramente variabilidad en los resutlados. Comprobarlo y ver como seria lo mas adecuado. hacerme una tabla en papel. esto tambien puede ser la diferencia con los intervalos significativos y que no sea lo mismo"""
    
nbins=12
i=0
for Fenotype in ["KO44", "KO179", "KO185"]:
    batch = "batch 5"
    df_temp = df.loc[(batch, "WT")].reset_index()
    #df_temp2 = df.loc[(batch, Fenotype, ["ZebraF_8", "ZebraF_9"])].reset_index()
    df_temp2 = df.loc[(batch, Fenotype)].reset_index()
    
    nbins = 12
    
    # Create a new figure and axis
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # # Plot 1: Blue histogram
    sns.histplot(
        data=df_temp,
        x="Dist_border_feret",
        hue="Fish",
        multiple="stack", # puede cambiarse a stack
        common_norm=True,
        element="poly",
        weights="weights_a_ind",
        stat="density",
        binrange=[0, 1],
        bins=nbins,
        palette="Blues",
        alpha=0.4,
        ax=ax,
        legend=False,  # Disable default legend for manual handling
    )
    
    # Get unique categories for Fish in df_temp
    fish_categories_temp = df_temp["Fish"].unique()
    
    # Manually create legend for Plot 1 (Blue)
    blue_palette = sns.color_palette("Blues", len(fish_categories_temp))
    blue_legend_elements = [
        plt.Line2D([0], [0], color=blue_palette[i], lw=4, label=cat)
        for i, cat in enumerate(fish_categories_temp)
    ]
    legend1 = ax.legend(
        handles=blue_legend_elements,
        loc="upper left",
        bbox_to_anchor=(1.05, 1.05),
        title="WT",
    )
    plt.gca().add_artist(legend1)  # Add the first legend manually
    
    # Plot 2: Orange histogram
    sns.histplot(
        data=df_temp2,
        x="Dist_border_feret",
        hue="Fish",
        multiple="stack",
        element="step",
        common_norm=True,
        weights="weights_a_ind",
        stat="density",
        binrange=[0, 1],
        bins=nbins,
        palette=["Reds", "Greens", "Oranges"][i],#"Oranges",
        alpha=0.3,
        ax=ax,
        legend=False,  # Disable default legend for manual handling
    )
    
    # Get unique categories for Fish in df_temp2
    fish_categories_temp2 = df_temp2["Fish"].unique()
    
    # Manually create legend for Plot 2 (Orange)
    orange_palette = sns.color_palette("Oranges", len(fish_categories_temp2))
    orange_legend_elements = [
        plt.Line2D([0], [0], color=orange_palette[i], lw=4, label=cat)
        for i, cat in enumerate(fish_categories_temp2)
    ]
    legend2 = ax.legend(
        handles=orange_legend_elements,
        loc="upper left",
        bbox_to_anchor=(1.05, 0.4),
        title=Fenotype,
    )
    
    # Set the title and layout
    ax.set_title("Ponderated Stacked histogram of radial position for " + batch)
    plt.tight_layout()
    
    # Show the plot
    plt.show()
    i=i+1


# %% Comparación de histogramas [md]
"""
### Comparación histogramas
Otra alternativa es generar manualmente el histograma y representarlo con barras de error (intervalo de confianza o sem).

Nota a posteriori: lo siguiente podría realizarse usando un Density Kernel Estimator, que calcule la distribución y sumar las distribuciones.
"""

# %%%% Generación Data Frame y agregación de histogramas
nbins = int(round(math.log(5601, 2), 1))  # Numero de Bins siguiendo regla

# funcion para calcular el hitograma con pesos
def compute_histogram(group, nbins):
    hist, bin_edges = np.histogram(
        group["Dist_border_feret"].to_numpy(),
        range=[0, 1],
        bins=nbins,
        density=True,
        weights=group["weights_a_ind"].to_numpy(),
    )
    return pd.Series({"hist": hist, "bins": bin_edges})

# Si solo quiero representar algun pez
# keep_fish = ["ZebraF_8", "ZebraF_9"]
# mask = (
#     (df.index.get_level_values("Fenotype") != "KO179") |
#     ((df.index.get_level_values("Fenotype") == "KO179") &
#         (df.index.get_level_values("Fish").isin(keep_fish))))

# DF con la altura del bin para cada zebra
distribution_df = (
    df.groupby(["Batch", "Fenotype", "Fish"])
    .apply(compute_histogram, nbins)
    .reset_index()
    .dropna()
)

distribution_df["bins"] = distribution_df["bins"].apply(
    lambda x: x[:-1]
)  # Para que bins tenga el mismo numero de elementos que 'hist'
distribution_df = distribution_df.explode(["hist", "bins"])

# Prueba con la anova, dibujo valores extremos para tenerlos
# mask=(distribution_df.Fenotype=="WT") & (distribution_df.bins==0)
# distribution_df.loc[mask, "hist"]= np.random.randint(7, 23, size=mask.sum())


# Agregación y calculo de la media de cada bin para cada histograma, aunque no lo uso despues
# Este DF contiene los intervalos de un test estadistico
distribution_df_agg = (
    distribution_df.groupby(["Batch", "Fenotype", "bins"])
    .agg(
        mean_hist=("hist", np.mean),
        sd_hist=("hist", np.std),
        confid_int=(
            "hist",
            lambda x: np.ptp(
                st.t.interval(
                    confidence=0.90,
                    df=len(x) - 1,
                    loc=np.mean(x),
                    scale=st.sem(x)  # scale = 1 significa usar la media en luygar de la sem
                    )
            ),
        ),
    )
    .reset_index()
    .dropna()
)


# DF con el valor de una anova para cada bin.
    
def anova_by_fenotype(x):
    # x es la serie "hist" dentro de cada grupo Batch+bins
    # obtenemos los valores de hist separados por fenotype
    subgroups = [group["hist"].values for _, group in x.groupby("Fenotype")]
    if len(subgroups) > 1:
        return st.f_oneway(*subgroups).pvalue
    else:
        return np.nan  # si solo hay un fenotype no podemos calcular ANOVA

# Aplicamos groupby
distribution_df_anova = (
    distribution_df.groupby(["Batch", "bins"])
    .apply(lambda g: pd.Series({
        "mean_hist": g["hist"].mean(),
        "sd_hist": g["hist"].std(),
        "confid_int": np.ptp(
            st.t.interval(
                confidence=0.95,
                df=len(g["hist"]) - 1,
                loc=g["hist"].mean(),
                scale=st.sem(g["hist"]))),
        "anova_pval": anova_by_fenotype(g)})
        )
    ).reset_index().dropna()

distribution_df_anova["p_val_sig"] = distribution_df_anova["anova_pval"] < 0.05



# %%%% [md]
"""
En el siguiente plot agrego usando la media los histogramas normalizados de cada condición por batch. 
Las barras de error se corresponden a los intervalos al 90% de un t- test
"""

# %%%% Plot de la agregación de histogramas como linea con t.test CI
# Lo represento como una linea para ver los errores


# Si Quiero representar solo  2 peces de unoun mut
# keep_fish = ["ZebraF_8", "ZebraF_9"]
# mask = ((distribution_df["Fenotype"] != "KO179") |      
#     ((distribution_df["Fenotype"] == "KO179") & (distribution_df["Fish"].isin(keep_fish))))
# df_filtered = distribution_df[mask]

g = sns.FacetGrid(
    distribution_df,
    hue="Fenotype",
    hue_order=["KO44", "KO179", "KO185", "WT"],
    row="Batch",
    palette=["red", "green", "yellow", "blue"],
    sharex="col",
    sharey=False,
    height=3,
    aspect=4,
)

# g.fig.suptitle("Evolución temporal de todas las variables para un pez de ejemplo",
# fontsize=24, fontdict={"weight": "bold"})

g.map_dataframe(
    sns.lineplot,
    estimator="mean",
    x="bins",
    y="hist",
    #errorbar="sd"
    errorbar=(lambda x: st.t.interval(
                        confidence=0.90,
                        df=len(x) - 1,
                        loc=np.mean(x),
                        scale=st.sem(x)
                        ))
)

g.set(xticks=bins[:-1], xticklabels=[f"{x:.2f}" for x in bins][:-1])
g.add_legend(title="Fenotype", fontsize=17, markerscale=5)
g.set_axis_labels(fontsize=25)
g.fig.suptitle("Averaged distribution of radial position relative to edge")
plt.subplots_adjust(top=0.95)
plt.show()

# %%%% Mismo Plot pero median
# Lo represento como una linea para ver los errores

g = sns.FacetGrid(
    distribution_df,
    hue="Fenotype",
    hue_order=["WT", "KO44", "KO179", "KO185"],
    row="Batch",
    palette=["blue", "red", "green", "yellow"],
    sharex="col",
    sharey=False,
    height=4,
    aspect=4,
)

# g.fig.suptitle("Evolución temporal de todas las variables para un pez de ejemplo",
# fontsize=24, fontdict={"weight": "bold"})

g.map_dataframe(
    sns.lineplot,
    estimator="median",
    x="bins",
    y="hist",
    errorbar=("pi", 50),
)

g.add_legend()
g.set_axis_labels(fontsize=20)
g.fig.suptitle("Median distribution of radial position relative to edge")
plt.subplots_adjust(top=0.95)
plt.show()


# %%%% [md]
"""
Resumen: XXXX
"""
# %% Porcentaje de tiempo pegado al borde mediante threshold [md]
"""
##  Porcentaje de tiempo pegado al borde
Por último, voy a realizar un análisis del tiempo que pasa cercano al borde. Para ello hay que definir un threshold de cercania, así que estudiaremos la evolución del resultado con respecto a la distancia que consideramos.
"""

# A partir de aqui reusar el codigo de los coletazos

# %%% Plot Ejemplo threshold
batch = "batch 3"

df_temp = df.loc[(batch, "KO179", "ZebraF_2")].reset_index()
# df_temp = df.loc[(batch, "KO44", "ZebraF_4")].reset_index()

g = sns.lineplot(data=df_temp, x="Frame", y="Dist_border_feret")
g.axhline(0.2, color="red")
g.set_title("Tiempo en el borde - Dist = 0", size=25)
plt.show()

# %%% [md]
"""
Contando para cada Zebra el total del tiempo que pasa bajo el Threshold, obtenemos
"""

# %%% Comparación usando un threshold fijo

Variable_plot = "Dist_border_feret"
threshold = 0.05
time_over_Thr_df = (
    df.groupby(["Batch", "Fenotype", "Fish"])[Variable_plot]
    .agg(
        boder_time=lambda x: (x < threshold).sum(),
        border_time_perc=lambda x: 100 * (x < threshold).sum() / len(x),
    )
    .reset_index()
)  # .dropna()

# a = sns.boxplot(x="Fenotype", y="contracted_perc", data=time_over_Thr_df)
# a.set_title("Numero de tiempo replegado con Thr " + str(threshold))
# b = sns.stripplot(
#     x="Fenotype", y="contracted_perc", data=time_over_Thr_df, color="grey", size=8
# )
# plt.show()


# base boxplot - Código de chatgpt
grped_bplot = sns.catplot(
    x="Batch",
    y="border_time_perc",
    data=time_over_Thr_df,
    hue="Fenotype",
    kind="box",
    legend=True,
    showfliers=False,
    height=6,
    aspect=1.9,
    hue_order=["WT", "KO44", "KO179", "KO185"],
)
grped_bplot._legend.remove()
# overlay the stripplot (jittered points)
ax = grped_bplot.ax  # get the matplotlib Axes
sns.stripplot(
    x="Batch",
    y="border_time_perc",
    data=time_over_Thr_df,
    hue="Fenotype",
    dodge=True,
    jitter=True,
    marker="o",
    color="black",
    ax=ax,
    legend=False,
    hue_order=["WT", "KO44", "KO179", "KO185"],
)

# define the pairwise comparisons per Batch
pairs = []
for batch in time_over_Thr_df["Batch"].unique():
    pairs.extend([
        ((batch, "WT"), (batch, "KO44")),
        ((batch, "WT"), (batch, "KO179")),
        ((batch, "WT"), (batch, "KO185")),
    ])

# add the statistical annotations
annotator = Annotator(
    ax,
    pairs,
    data=time_over_Thr_df,
    x="Batch",
    y="border_time_perc",
    hue="Fenotype",
    hue_order=["WT", "KO44", "KO179", "KO185"]
)

annotator.configure(
    test="t-test_ind",       # can also use "Mann-Whitney", "Kruskal"
    text_format="star",      # or "simple" to show ns/p-values
    loc="inside",            # inside or outside bars
    #comparisons_correction="bonferroni"
)

annotator.apply_and_annotate()

ax.set_title(
    f"Porcentaje del tiempo que pasa el pez cerca del borde - Threshold = {threshold}",
    size=20
)

ax.legend(
    title="Fenotype",
    bbox_to_anchor=(0.95, 1),  # place it outside to the right
    loc='upper left'
)

plt.tight_layout()
plt.show()




# %%% Evolución del resultado con el threshold [md]
"""
### Evolución del resultado con el threshold
Dado que este resultado es sensible al Threshold, vamos a ver como evoluciona el resultado con el Threshold
elegido. Se representa la diferencia de la mediana por batch del tiempo que pasa replegado el KO con respecto a su mutante. 
"""

# %%% Construcción del DF de evolución del resultado con el threshold

Variable_plot = "Dist_border_feret"

threshold_result = pd.DataFrame(
    columns=["Threshold", "Batch", "Fenotype", "Mean_diff", "Median_diff", "CI", "IQR_high_dif", "IQR_low_dif"]
)
ref = ko44 = ko179 = np.nan

for thr in np.arange(0.0, 0.2, 0.02):  # iteración sobre el threshold
    # data frame con los valores para ese threshold
    time_over_Thr = (
        df.groupby(["Batch", "Fenotype", "Fish"])[Variable_plot]
        .agg(
            contracted=lambda x: (x < thr).sum(),
            contracted_perc=lambda x: 100 * (x < thr).sum() / len(x),
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
        )  # los valores de los WT

        grouped_batch = group.dropna().drop("Batch", axis=1).groupby("Fenotype")
        for fenotype, group in grouped_batch:  # loop over each possible fenotype
            if fenotype != "WT":
                mut = group.contracted_perc
                new_row = {
                    "Threshold": thr,
                    "Batch": batch,
                    "Fenotype": fenotype,
                    "Mean_diff": np.mean(WT) - np.mean(mut),                    
                    "Median_diff": np.median(WT) - np.median(mut),
                    "CI": np.ptp(
                        stats.ttest_ind(WT, mut).confidence_interval(
                            confidence_level=0.90  # El ancho del intervalo del confianza en el plot subsiguiente
                        )
                    ) / 2,
                    "IQR_high_dif": np.percentile(WT, 75) - np.percentile(mut, 75),
                    "IQR_low_dif": np.percentile(WT, 25) - np.percentile(mut, 25),
                }
                threshold_result.loc[len(threshold_result)] = new_row


# %%% Plot del resultado frente al threshold con Intervalos de confianza

threshold_result["hue"] = threshold_result.Batch + " - " + threshold_result.Fenotype

df_plot = threshold_result
# [
#     threshold_result["Fenotype"].isin(["KO179"])
# ]  # filtro para lo que se quiere repesentar
df_plot["CI_up"] = df_plot.Mean_diff + df_plot.CI
df_plot["CI_down"] = df_plot.Mean_diff - df_plot.CI
df_plot["IQR_dif_up"] = df_plot.Median_diff + df_plot.IQR_high_dif
df_plot["IQR_dif_down"] = df_plot.Median_diff - df_plot.IQR_low_dif

for f in ["KO44", "KO179", "KO185"]:
    df_temp = df_plot.loc[df_plot.Fenotype == f]
    g = sns.lineplot(
        data=df_temp,
        x="Threshold",
        y="Mean_diff",
        hue="hue",
    )
    g.set_title("Evolución del resultado (diferencia de medias con t-test CI) con el Threshold")
    
    for hue in df_plot.hue.unique():
        g.fill_between(
            x="Threshold",
            y1="CI_up",
            y2="CI_down",
            alpha=0.1,
            data=df_temp[df_temp.hue == hue],
        )
    plt.show()
    

# g = sns.lineplot(
#     data=df_plot,
#     x="Threshold",
#     y="Median_diff",
#     hue="hue",
# )
# g.set_title("Evolución del resultado (diferencia de medianas Q75 y Q25) con el Threshold (la sombra significa otra cosa)")

# for hue in df_plot.hue.unique():
#     g.fill_between(
#         x="Threshold",
#         y1="IQR_dif_up",
#         y2="IQR_dif_down",
#         alpha=0.1,
#         data=df_plot[df_plot.hue == hue],
#     )
# plt.show()

# %%%% [md]
"""
En esta gráfica camos a ver si de verdad hay una diferencia en el porcentaje de tiempo que pasan cerca del borde
"""

# %% Distancia como variable extensiva [md]
"""
# Distancia como variable extensiva

Como sugerencia de Tomás, voy a calcular el sumatorio de la distancia radial al centro como si fuera una variable extensiva y representar la distribución mediante boxplots. Tratar la distancia al centro como variable extensiva nos indica la "cantidad de tiempo" que ha pasado lejos del centro/cerca del borde. Posiblemente este comportamiento lo podrémos haber detectado ya con el análisis anterior. 

"""

# %%% Calculo del DF de la variable extensiva

Distancia_acumulada = (
    df.groupby(["Batch", "Fenotype", "Fish"])["Dist_border_feret"]
    .sum()
    .rename("Distancia_acumulada")
    .reset_index()
)
Distancia_acumulada = Distancia_acumulada.loc[
    Distancia_acumulada.Distancia_acumulada != 0
]

grped_bplot = sns.catplot(
    x="Batch",
    y="Distancia_acumulada",
    data=Distancia_acumulada,
    hue="Fenotype",
    kind="box",
    legend=True,
    showfliers=False,
    height=6,
    aspect=1.9,
    hue_order=["WT", "KO44", "KO179", "KO185"],
)

#grped_bplot._legend.remove()
# overlay the stripplot (jittered points)
ax = grped_bplot.ax  # get the matplotlib Axes

# make grouped stripplot
grped_bplot = sns.stripplot(
    x="Batch",
    y="Distancia_acumulada",
    data=Distancia_acumulada,
    hue="Fenotype",
    jitter=True,
    legend=False,
    dodge=True,
    marker="o",
    color="black",
    # palette="Set2",
    hue_order=["WT", "KO44", "KO179", "KO185"],
)
handles, labels = grped_bplot.get_legend_handles_labels()
grped_bplot.set_title(
    "Distancia radial acumulada como variable extensiva",
    size=20,
)

# define the pairwise comparisons per Batch
pairs = []
for batch in Distancia_acumulada["Batch"].unique():
    pairs.extend([
        ((batch, "WT"), (batch, "KO44")),
        ((batch, "WT"), (batch, "KO179")),
        ((batch, "WT"), (batch, "KO185"))
    ])

# add the statistical annotations
annotator = Annotator(
    ax,
    pairs,
    data=Distancia_acumulada,
    x="Batch",
    y="Distancia_acumulada",
    hue="Fenotype",
    hue_order=["WT", "KO44", "KO179", "KO185"]
)

annotator.configure(
    test="t-test_ind",       # can also use "Mann-Whitney", "Kruskal"
    text_format="star",      # or "simple" to show ns/p-values
    loc="inside"            # inside or outside bars

)

annotator.apply_and_annotate()

plt.show()

# %% Correlaciones [md]
"""
Como último recurso voy a graficar las relaciones entre distancia y la posición al borde
"""

# %%% Merge DF
correlaciones_df = Distancia_acumulada.merge(Distancia_recorrida, on=['Batch', 'Fenotype', 'Fish'], how='inner')
correlaciones_df = correlaciones_df.merge(time_over_Thr_df, on=['Batch', 'Fenotype', 'Fish'], how='inner')

# %%% Correlation Plot
batches_2_plot = ["batch 1","batch 2","batch 3", "batch 4","batch 5","batch 6","batch 7","batch 8"]
#batches_2_plot = ["batch 1","batch 2","batch 3", "batch 4"]
batches_2_plot = ["batch 4","batch 5","batch 6","batch 7","batch 8"]

df_plot = correlaciones_df.loc[correlaciones_df.Batch.isin(batches_2_plot)]
df_plot["Batch"] = df_plot["Batch"].cat.remove_unused_categories()


g = sns.lmplot(data=df_plot, x='Distancia_recorrida', y='Distancia_acumulada', hue='Fenotype',
           hue_order=["WT", "KO44", "KO179", "KO185"],
           ci=95)
g.fig.suptitle(
    "Distancia radial acumulada vs Distancia Recorrida",
    size=10)
plt.show()

g = sns.lmplot(data=df_plot, x='Distancia_recorrida', y='border_time_perc', hue='Fenotype',
           hue_order=["WT", "KO44", "KO179", "KO185"],
           ci=95)
g.fig.suptitle(
    "Distancia radial acumulada vs Distancia Recorrida",
    size=10)
plt.show()

g = sns.lmplot(data=df_plot, x='Distancia_acumulada', y='border_time_perc', hue='Fenotype',
           hue_order=["WT", "KO44", "KO179", "KO185"],
           ci=95)
g.fig.suptitle(
    "Distancia radial acumulada vs Distancia Recorrida",
    size=10)
plt.show()




# %% Conclusiones [md]
"""
He comprobado y No hay diferencia entre los dos métodos de análisis de los videos (strict y maximun maxima).
"""
