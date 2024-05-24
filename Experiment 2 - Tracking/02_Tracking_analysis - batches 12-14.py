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
# get_ipython().magic("reset -sf")

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

# %%% Load Files batches 1-11

if platform.system() == "Windows":
    folder_path = "P:\CABD\Lab Ozren\Marta Fernandez\Behavioral Assays Batches 1-11 Results\Experimento Tracking"
else:
    folder_path = "/home/ale/pCloudDrive/CABD/Lab Ozren/Marta Fernandez/Behavioral Assays Batches 1-11 Results/Experimento Tracking/"

files = get_files_in_folder(folder_path)

df1 = []
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
        df1.append(csv)
        del (csv, f)

df1 = pd.concat(df1)
df1.Time = (df1.Frame - 1) / 5  # 5 fps
# df1 = df1.drop(df1.columns[-1], axis=1)


# renombro KO a KO44 para el batch 6 y 7
df1.loc[df1.Fenotype == "KO", "Fenotype"] = "KO44"
df1.loc[df1.Fenotype == "KO 44", "Fenotype"] = "KO44"
df1.loc[df1.Fenotype == "KO 179", "Fenotype"] = "KO179"

df1["Batch"] = pd.Categorical(
    df1["Batch"],
    categories=["batch 6", "batch 7", "batch 8", "batch 11"],
    ordered=True,
)

df1["DF"] = "DF1"  # Para identificar el grupo de batches

elements = round(
    pd.crosstab(index=df1.Batch, columns=df1.Fenotype) / 4200
)  # divided by lengh of the video
print(str(elements).replace(".0", "").replace("],", "]\n"))

# %%%% Eliminar filas que no son una medida correcta por algún motivo de ImageJ
# aparecen como caddena de texto en la columna frame
df1["Frame"] = pd.to_numeric(df1["Frame"], errors="coerce")
df1 = df1.dropna(subset=["Frame"], how="any", axis=0)

# %%%% Variables auxiliares
"""
Con el dataset limpio genero unas variables auxiliares
"""
# Para evaluar la distancia al borde de 0 (border) a 1 (center)
df1["Dist_border"] = (
    df1[["Mean-Distance"]] / 255
)  # Radio del pocillo para estos batches
# Alternativamente la distancia al centro 1=border a 0=center
df1["Dist_center"] = abs(df1["Dist_border"] - 1)
# Variable auxiliar
df1["Feno_Batch"] = df1.Fenotype.astype(str) + "_" + df1.Batch.astype(str)


# %%% Load Files batches 12-14

if platform.system() == "Windows":
    folder_path = "P:\CABD\Lab Ozren\Marta Fernandez\Behavioral Assays Batches 12-14 Results\Experimento Tracking"
else:
    folder_path = "/home/ale/pCloudDrive/CABD/Lab Ozren/Marta Fernandez/Behavioral Assays Batches 12-14 Results/Experimento Tracking/"

files = get_files_in_folder(folder_path)

df2 = []
for f in files:
    if ".csv" in f:
        csv = pd.read_csv(f, sep=";")
        csv.insert(0, "Batch", re.search("batch \d+", f.lower()).group(0))
        csv.insert(1, "Fenotype", re.search("(KO\d*|WT)", f.upper()).group(1))
        csv.insert(2, "Fish", "ZebraF_" + re.search("(\d+)(.tif)", f.lower()).group(1))
        df2.append(csv)
        del (csv, f)

df2 = pd.concat(df2)
df2.Time = (df2.Frame - 1) / 6  # 6 fps
# df2 = df2.drop(df2.columns[-1], axis=1)

df2["DF"] = "DF2"  # Para identificalo luego

df2["Batch"] = pd.Categorical(
    df2["Batch"],
    categories=["batch 12", "batch 13", "batch 14"],
    ordered=True,
)

elements = round(
    pd.crosstab(index=df2.Batch, columns=df2.Fenotype) / 5601
)  # divided by lengh of the video
print(str(elements).replace(".0", "").replace("],", "]\n"))

# %%%% Eliminar filas que no son una medida correcta por algún motivo de ImageJ
# aparecen como cadena de texto en la columna frame
df2["Frame"] = pd.to_numeric(df2["Frame"], errors="coerce")
df2 = df2.dropna(subset=["Frame"], how="any", axis=0)

# %%%% Variables auxiliares
"""
Con el dataset limpio genero unas variables auxiliares. La distancia la normalizo a 1, siendo 0 el borde y 1 el centro del pocillo.
"""

# Para evaluar la distancia al borde de 0 (border) a 1 (center)

df2["Dist_border"] = (
    df2[["Mean-Distance"]] / 170
)  # valor del radio del pocillo medido de las imagenes para batches 12-14
df2["Dist_center"] = abs(df2["Dist_border"] - 1)
# Variable auxiliar
df2["Feno_Batch"] = df2.Fenotype.astype(str) + "_" + df2.Batch.astype(str)

# %%% Union DF de los distintos batches
"""
Ojo, el fps de ambas muestras no es el mismo, por lo que hay que considerarlo a la hora de sacar resultados
"""
df = pd.concat([df1, df2])

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

df[["X", "Y", "Mean-Distance"]] = df[["X", "Y", "Mean-Distance"]].interpolate(
    method="linear"
)

# %%% [md]
"""
Este análisis se ha realizado usando strict. Hay NAs pero nigun pez tiene demasiados si consideramos que hemos medido miles de frames. Se han imputado
"""
# %% Filtrado de datos debido a detección de otras particulas[md]
"""
# Filtrado de los datos
Voy a eliminar los frames en los que se ha detectado un salto de posición demasiado alto ya que posiblemente se debera al haber detectado una posición anomala en ese frame. P.ejem el pez de ejemplo tiene uno de estos eventos
"""
# %%% Histograma de las distancias rocorridas
df["X_diff"] = df.groupby(["Batch", "Fenotype", "Fish"]).X.diff()
df["Y_diff"] = df.groupby(["Batch", "Fenotype", "Fish"]).Y.diff()
df["dist"] = np.sqrt((df.X_diff**2) + (df.Y_diff**2))

df_temp = df[
    (df.Batch == "batch 12") & (df.Fenotype == "KO44") & (df.Fish == "ZebraF_15")
]

g = sns.histplot(
    data=df_temp, x="dist", stat="density", log_scale=(False, True), binwidth=9
)
g.set_title("Distribution of frame movement of a single zebra")
plt.show()
# %%% [md]
"""
Se aprecia como hay movimientos que la distancia recorrida es muy alta. Esto solo puede deberse a un error en la detección del gusano. Considero que 200px es el límite de salto por frame, aunque bien podrían ser menos.
Se eliminan e imputan.
"""

# %%% Imputación anomalias

"""
Como se recalculan las distancias entre 2 frames, es posible que una posición outlier aparezca en varios frames, pero solo se elimina una vez, por lo que al recalcular las distancias apareceran de nuevo estos saltos. Para corregirlo se realiza el proceso de filtrado varias veces.
"""

# Para ambos DF
for i in range(10):
    # Imputación en las columnas que tienen medidas
    df.loc[(df.DF == "DF1") & (df.dist > 220), ("X", "Y", "Mean-Distance")] = np.nan
    df.loc[(df.DF == "DF2") & (df.dist > 180), ("X", "Y", "Mean-Distance")] = (
        np.nan
    )  # ya el tamaño de los pocillso son diferentes
    # imputación por interpolación de los cercanos
    df[["X", "Y", "Mean-Distance"]] = df[["X", "Y", "Mean-Distance"]].interpolate(
        method="linear"
    )

    # Recalculo las distancias
    df["X_diff"] = df.groupby(["Batch", "Fenotype", "Fish"]).X.diff()
    df["Y_diff"] = df.groupby(["Batch", "Fenotype", "Fish"]).Y.diff()
    df["dist"] = np.sqrt((df.X_diff**2) + (df.Y_diff**2))


df_temp = df[
    (df.Batch == "batch 12") & (df.Fenotype == "KO44") & (df.Fish == "ZebraF_15")
]

g = sns.histplot(
    data=df_temp, x="dist", stat="density", log_scale=(False, True), binwidth=9
)
g.set_title("Distribution of frame movement of a single zebra")
plt.show()


# %% Distancia Recorrida [md]
"""
## Distancia Recorrida
Se calcula la distancia que recorre el pez a lo largo del video y se gráfica por batch
"""


# %%% Calculo de la distancia recorrida

# dataframe con la distancia recorrida por el  gusano
Dist = (
    df.dropna()
    .groupby(["Batch", "Fenotype", "Fish"], as_index=False)
    .dist.sum(min_count=0)
    .round()
)  # .reset_index()

Dist = Dist.loc[
    Dist.dist != 0
]  # Importante. Elimina los Zebra que corresponden a categorias de las que no hay datos, ya que el groupby las genera
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

# %% Posición del Pez sobre el video [md]
""" 
# Posición del Zebra en el Pocillo
A lo largo del video, el pez se posiciona en algún lugar de la placa. Se estima que la posición que mantiene el Zebra es comportamental, por lo que vamos a estudiar que posición mantiene con respecto al borde d=0 hasta el centro d=1. (Dado que existe simetria radial, solo usamos la distancia respecto al borde)
"""

# %%% Histograma 1 Zebra [md]
"""
## Histograma de 1 solo Pez
Vamos a evaluar el histograma de 1 Zebra. Este nos va a indicar donde se posiciona el Zebra a lo largo del tiempo del video.
"""
# %%% Grafico Histograma 1 Zebra
df_temp = df[(df.Batch == "batch 14") & (df.Fenotype == "WT") & (df.Fish == "ZebraF_3")]

g = sns.histplot(
    data=df_temp, x="Dist_border", stat="density", binrange=[0, 1], bins=12
)
g.set_title("Distribution of radial position relative to edge of a single zebra")

plt.show()

# %%% [md]
"""
Se observa como el Zebra se posiciona mayormente a lo largo de la franja central, como esperable probabilisticamente, sin posicionarse apenas cerca del borde. 

## El Histograma debe normalizarse y luego lo podremos agrupar

"""
# %%% Histograma por Condición [md]
"""
## Histograma Condición 
Voy a ver si, en media, un fenotipo cambia su modo de distribuirse en el pocillo acumulando los histogramas. Esto no es un problema ya que como todos los videos duran lo mismo, pueden sumarse los counts y sera equivalente a sumar la densidad de probabilidades. 

"""
# %%% Dibujado con Histplot por Batch


g = sns.FacetGrid(
    df,
    row="Batch",
    hue="Fenotype",
    hue_order=["WT", "KO44", "KO179"],
    palette="pastel",
    sharex="col",
    sharey=False,
    height=5,
    aspect=4,
)

# g.fig.suptitle("Evolución temporal de todas las variables para un pez de ejemplo",
# fontsize=24, fontdict={"weight": "bold"})

g.map_dataframe(
    sns.histplot,
    x="Dist_border",
    element="step",
    edgecolor="black",
    binrange=[0, 1],
    # cumulative=True,
    bins=10,
    stat="density",
    common_norm=False,
    kde=True,
    kde_kws={"bw_adjust": 0.8},
)
g.add_legend()
g.set_axis_labels(fontsize=20)
g.fig.suptitle("Accumulated distribution of radial position relative to edge")
plt.subplots_adjust(top=0.95)
plt.show()
# %%% [md]
"""
Este gráfico muestra la densidad de probabilidad para cada condición y batch, normalizada para cada condición. He agregado los Zebra ya que como cada uno dura lo mismo, puede hacerse ya que todos tendrán el mismo peso.

En los casos en los que la distribución de una condición es significativamente diferente a las otras habria que estudiar que la contribución individual de cada Zebra sea razonablemente similar, y no que un solo Zebra sea el que produce la desviacion del histograma o PDF. Lo vemos
"""

# %%% Histograma acumulado por Zebra
# La clave de estos histogramas es que cada sns.hisplot es una capa independiente
batch = "batch 12"
df_temp = df[(df.Batch == batch) & (df.Fenotype == "WT")]
df_temp2 = df[(df.Batch == batch) & (df.Fenotype == "KO44")]

sns.histplot(
    data=df_temp,
    x="Dist_border",
    hue="Fish",
    multiple="stack",
    common_norm=True,
    element="poly",
    stat="density",
    binrange=[0, 1],
    bins=12,
    palette="Blues",
    alpha=0.4,
)

g = sns.histplot(
    data=df_temp2,
    x="Dist_border",
    hue="Fish",
    multiple="stack",
    element="step",
    common_norm=True,
    stat="density",
    binrange=[0, 1],
    bins=12,
    palette="Oranges",
    alpha=0.3,
)
g.set_title("Accumulated histogram of radial position relative to edge")

plt.show()

# %%% [md]
"""
Resulta un plot bastante sucio pero se aprecia la diferencia. Creo que una mejor alternativa será representar unicamente la condición y el batch de interes para ver la intra-distribución de los Zebra y ver que son razonablemente homogeneos y no se debe a un Zebra Oulier.

### Normalización manual histograma
Otra alternativa es generar manualmente el histograma y representarlo con barras de error.

Nota a posteriori: lo siguiente podría realizarse usando un Density Kernel Estimator, que calcule la distribución y sumar las distribuciones.
"""

# %%% Normalización manual histograma


# %%%% Generación Data Frame y agregación de histogramas

nbins = int(round(math.log(5601, 2), 1))  # Numero de Bins siguiendo regla

distribution_df = (
    df.groupby(["Batch", "Fenotype", "Fish"], as_index=False)
    .Dist_border.apply(
        lambda x: np.histogram(x, range=[0, 1], bins=nbins, density=True)
    )
    .dropna()
)

distribution_df["hist"], distribution_df["bins"] = zip(*distribution_df.Dist_border)

distribution_df["bins"] = distribution_df["bins"].apply(
    lambda x: x[:-1]
)  # Para que tenga el mismo numero de elementos que 'hist'
distribution_df = distribution_df.drop("Dist_border", axis=1)
distribution_df = distribution_df.explode(["hist", "bins"])


distribution_df_agg = (
    distribution_df.groupby(["Batch", "Fenotype", "bins"])
    .agg(
        mean_hist=("hist", np.mean),
        sd_hist=("hist", np.std),
        confid_int=(
            "hist",
            lambda x: np.ptp(
                st.t.interval(
                    confidence=0.95,
                    df=len(x) - 1,
                    loc=np.mean(x),
                    scale=st.sem(
                        x
                    ),  # scale = 1 significa usar la media en luygar de la sem
                )
            ),
        ),
    )
    .reset_index()
    .dropna()
)

# %%%% Plot de la agregación de histogramas como linea
# Lo represento como una linea para ver los errores

g = sns.FacetGrid(
    distribution_df,
    hue="Fenotype",
    hue_order=["WT", "KO44", "KO179"],
    row="Batch",
    palette="pastel",
    sharex="col",
    sharey=False,
    height=4,
    aspect=3,
)

# g.fig.suptitle("Evolución temporal de todas las variables para un pez de ejemplo",
# fontsize=24, fontdict={"weight": "bold"})

g.map_dataframe(
    sns.lineplot,
    estimator="mean",
    x="bins",
    y="hist",
    errorbar=("ci", 90),
)

g.add_legend()
g.set_axis_labels(fontsize=20)
g.fig.suptitle("Averaged distribution of radial position relative to edge")
plt.subplots_adjust(top=0.95)
plt.show()

# %%%% [md]
"""
Espero de este gráfico encontrar diferencias consistentes
"""
# %%% Porcentaje de tiempo pegado al borde [md]
"""
###  Porcentaje de tiempo pegado al borde
Por último, voy a realizar un análisis del tiempo que pasa cercano al borde. Para ello hay que definir un threshold de cercania, así que estudiaremos la evolución del resultado con respecto a la distancia que consideramos.
"""

# A partir de aqui reusar el codigo de los coletazos

# %%%% Plot Ejemplo threshold

df_temp = df[
    (df.Batch == "batch 12") & (df.Fenotype == "KO44") & (df.Fish == "ZebraF_14")
]

g = sns.lineplot(data=df_temp, x="Frame", y="Dist_border")
g.axhline(0.2, color="red")
g.set_title("Tiempo en el borde - Dist = 0", size=25)
plt.show()

# %%%% [md]
"""
Contando para cada Zebra el total del tiempo que pasa bajo el Threshold, obtenemos
"""

# %%%% Comparación usando un threshold fijo

Variable_plot = "Dist_border"
threshold = 0.1
time_over_Thr = (
    df.groupby(["Batch", "Fenotype", "Fish"])[Variable_plot]
    .apply(lambda x: (x < threshold).sum())
    .reset_index()
    .rename(columns={Variable_plot: "boder_time"})
).dropna()

time_over_Thr["boder_time_cent"] = (
    100 * time_over_Thr.boder_time / 4202
)  # longitud del video

# a = sns.boxplot(x="Fenotype", y="contracted_perc", data=time_over_Thr)
# a.set_title("Numero de tiempo replegado con Thr " + str(threshold))
# b = sns.stripplot(
#     x="Fenotype", y="contracted_perc", data=time_over_Thr, color="grey", size=8
# )
# plt.show()

grped_bplot = sns.catplot(
    x="Batch",
    y="boder_time_cent",
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
    y="boder_time_cent",
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
    "Porcentaje del tiempo que pasa el gusano cerca del borde - Threshold = "
    + str(threshold),
    size=20,
)
plt.legend(handles[0:3], labels[0:3])
plt.show()


# %%%% Evolución del resultado con el threshold [md]
"""
### Evolución del resultado con el threshold
Dado que este resultado es sensible al Threshold, vamos a ver como evoluciona el resultado con el Threshold
elegido. Se representa la diferencia de la mediana por batch del tiempo que pasa replegado el KO con respecto a su mutante. 
"""

# %%%% Construcción del DF de evolución del resultado con el threshold

Variable_plot = "Dist_border"

threshold_result = pd.DataFrame(
    columns=["Threshold", "Batch", "Fenotype", "Mean_diff", "CI"]
)
ref = ko44 = ko179 = np.nan

for thr in np.arange(0.0, 0.30, 0.02):  # iteración sobre el threshold
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
                    "CI": np.ptp(
                        stats.ttest_ind(WT, mut).confidence_interval(
                            confidence_level=0.80  # El ancho del intervalo del confianza en el plot subsiguiente
                        )
                    )
                    / 2,
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

g = sns.lineplot(
    data=df_plot,
    x="Threshold",
    y="Mean_diff",
    hue="hue",
)
g.set_title("Evolución del resultado (diferencia de medias) con el threshold Threshold")

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
En esta gráfica camos a ver si de verdad hay una diferencia en el porcentaje de tiempo que pasan cerca del borde
"""

# %% Distancia como variable extensiva [md]
"""
# Distancia como variable extensiva

Como sugerencia de Tomás, voy a calcular el sumatorio de la distancia radial al centro como si fuera una variable extensiva y representar la distribución mediante boxplots. Tratar la distancia al centro como variable extensiva nos indica la "cantidad de tiempo" que ha pasado lejos del centro/cerca del borde. Posiblemente este comportamiento lo podrémos haber detectado ya con el análisis anterior. 

"""

# %%% Calculo del DF de la variable extensiva

Distancia_acumulada = (
    df.groupby(["Batch", "Fenotype", "Fish"])["Dist_center"]
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
    legend=False,
    showfliers=False,
    height=6,
    aspect=1.9,
    hue_order=["WT", "KO44", "KO179"],
)
# make grouped stripplot
grped_bplot = sns.stripplot(
    x="Batch",
    y="Distancia_acumulada",
    data=Distancia_acumulada,
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
    "Distancia radial acumulada como variable extensiva",
    size=20,
)
plt.legend(handles[0:3], labels[0:3])
plt.show()

# %%% Conclusiones [md]
"""
He comprobado y No hay diferencia entre los dos métodos de análisis de los videos (strict y maximun maxima).
El batch 11 no se ha posido procesar bien


"""
