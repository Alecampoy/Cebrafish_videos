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

windows = False
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

# Para evaluar la distancia al centro de 0 a 1
df["Dist_center"] = df[["Mean-Distance"]] / 255

# Variable auxiliar
df["Feno_Batch"] = df.Fenotype.astype(str) + "_" + df.Batch.astype(str)


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

df[["X", "Y", "Dist_center"]] = df[["X", "Y", "Dist_center"]].interpolate(
    method="linear"
)

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

# %% Posición del Pez sobre el video [md]
""" 
# Posición del Zebra en el Pocillo
A lo largo del video, el pez se posiciona en algún lugar de la placa. Se estima que la posición que mantiene el Zebra es comportamental, por lo que vamos a estudiar que posición mantiene con respecto al borde d=0 hasta el centro d=255. (Dado que existe simetria radial, solo usamos la distancia respecto al borde)
"""

# %%% Histograma 1 Zebra [md]
"""
## Histograma de 1 solo Pez
Vamos a evaluar el histograma de 1 Zebra. Este nos va a indicar donde se posiciona el Zebra a lo largo del tiempo del video.
"""
# %%% Grafico Histograma 1 Zebra
df_temp = df[(df.Batch == "batch 7") & (df.Fenotype == "WT") & (df.Fish == "ZebraF_1")]

sns.histplot(data=df_temp, x="Dist_center", stat="percent", bins=20)
plt.show()

# %%% [md]
"""
Se observa como el Zebra se posiciona mayormente a lo largo de la franja centra, como esperable probabilisticamente, sin posicionarse apenas cerca del borde. 

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
    x="Dist_center",
    element="step",
    edgecolor="black",
    binrange=[0, 1],
    bins=10,
    stat="density",
    common_norm=False,
    kde=True,
    kde_kws={"bw_adjust": 1},
)
g.add_legend()
g.set_axis_labels(fontsize=20)
# %%% [md]
"""
Este gráfico muestra la densidad de probabilidad para cada condición y batch, normalizada para cada condición. He agregado los Zebra ya que como cada uno dura lo mismo, puede hacerse ya que todos tendrán el mismo peso.

"""
# %%% Normalización manual histograma

perc, bins = np.histogram(df_temp["Dist_center"], range=[0, 1], bins=10, density=True)
sns.barplot(x=bins.round(2)[:-1], y=perc, width=1, edgecolor="black", color="lightblue")
plt.show()

# %%% Generación Data Frame para agregar

a = (
    df.groupby(["Batch", "Fenotype", "Fish"])
    .Dist_center.apply(lambda x: np.histogram(x, range=[0, 1], bins=10))
    .dropna()
    .reset_index()
)

a["d"], a["e"] = zip(*a.Dist_center)
a = a.drop("Dist_center", axis=1)

c = a.explode(list("de"))


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
