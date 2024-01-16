# -*- coding: utf-8 -*-

# %% Intro [md]
"""
# Tracking Experiment
# **Análisis del Tracking del Movimiento de Zebrafish en una placa**
### Author: Alejandro Campoy Lopez  
"""

# %% Librerias
import warnings, math
import pandas as pd
import re
import os, platform
import glob
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as st


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

if platform.system() == "Windows":
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

sns.histplot(data=df_temp, x="Dist_center", stat="density", binrange=[0, 1], bins=12)
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
    # cumulative=True,
    bins=10,
    stat="density",
    common_norm=False,
    kde=True,
    kde_kws={"bw_adjust": 0.8},
)
g.add_legend()
g.set_axis_labels(fontsize=20)
plt.show()
# %%% [md]
"""
Este gráfico muestra la densidad de probabilidad para cada condición y batch, normalizada para cada condición. He agregado los Zebra ya que como cada uno dura lo mismo, puede hacerse ya que todos tendrán el mismo peso.

En los casos en los que la distribución de una condición es significativamente diferente a las otras habria que estudiar que la contribución individual de cada Zebra sea razonablemente similar, y no que un solo Zebra sea el que produce la desviacion del histograma o PDF. Lo vemos
"""

# %%% Histograma acumulado por Zebra
# La clave de estos histogramas es que cada sns.hisplot es una capa independiente

df_temp = df[(df.Batch == "batch 7") & (df.Fenotype == "WT")]
df_temp2 = df[(df.Batch == "batch 7") & (df.Fenotype == "KO44")]

sns.histplot(
    data=df_temp,
    x="Dist_center",
    hue="Fish",
    multiple="stack",
    common_norm=True,
    element="poly",
    stat="density",
    binrange=[0, 1],
    bins=12,
    palette="Reds",
    alpha=0.4,
)

sns.histplot(
    data=df_temp2,
    x="Dist_center",
    hue="Fish",
    multiple="stack",
    element="step",
    common_norm=True,
    stat="density",
    binrange=[0, 1],
    bins=12,
    palette="Greens",
    alpha=0.3,
)
plt.show()

# %%% [md]
"""
Resulta un plot bastante sucio. Creo que una mejor alternativa será representar unicamente la condición y el batch de interes para ver la intra-distribución de los Zebra y ver que no se debe a un Zebra Oulier, sino que son razonablemente homogeneos.

### Normalización manual histograma
Otra alternativa es generar manualmente el histograma y representarlo con barras de error.

Nota a posteriori: lo siguiente podría realizarse usando un Kernel Estimator
"""

# %%% Normalización manual histograma


# %%%% Generación Data Frame y agregación de histogramas

nbins = int(round(math.log(4200, 2), 1))  # Numero de Bins siguiendo regla

distribution_df = (
    df.groupby(["Batch", "Fenotype", "Fish"], as_index=False)
    .Dist_center.apply(
        lambda x: np.histogram(x, range=[0, 1], bins=nbins, density=True)
    )
    .dropna()
)

distribution_df["hist"], distribution_df["bins"] = zip(*distribution_df.Dist_center)

distribution_df["bins"] = distribution_df["bins"].apply(
    lambda x: x[:-1]
)  # Para que tenga el mismo numero de elementos que 'hist'
distribution_df = distribution_df.drop("Dist_center", axis=1)
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
                    confidence=0.95, df=len(x) - 1, loc=np.mean(x), scale=st.sem(x)
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
    (df.Batch == "batch 7") & (df.Fenotype == "KO44") & (df.Fish == "ZebraF_4")
]

g = sns.lineplot(data=df_temp, x="Time", y="Dist_center")
g.axhline(0.2, color="red")
g.set_title("Tiempo en el borde Dist=0", size=25)
plt.show()

# %%%% [md]
"""
Contando para cada Zebra el total del tiempo que pasa bajo el Threshold, obtenemos
"""

# %%%% Comparación usando un threshold fijo

Variable_plot = "Dist_center"
threshold = 0.20
time_over_Thr = (
    df.groupby(["Batch", "Fenotype", "Fish"])[Variable_plot]
    .apply(lambda x: (x < threshold).sum())
    .reset_index()
    .rename(columns={Variable_plot: "boder_time"})
)

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
    "Porcentaje del tiempo que pasa el gusano cerca del borde - bajo el Threshold = "
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

# %%%% plot del threshold
threshold_result = pd.DataFrame(columns=["Threshold", "Batch", "KO44", "KO179"])

i = 0
Variable_plot = "Dist_center"
for thr in np.arange(0, 0.51, 0.01):
    time_over_Thr = (
        df.groupby(["Batch", "Fenotype", "Fish"])[Variable_plot]
        .apply(lambda x: (x < thr).sum())
        .reset_index()
        .rename(columns={Variable_plot: "boder_time"})
    )
    time_over_Thr["boder_time_perc"] = 100 * time_over_Thr.boder_time / 4200

    df_temp_median = time_over_Thr.groupby(["Batch", "Fenotype"])[
        "boder_time_perc"
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
    threshold_result.loc[i + 3] = [
        thr,
        "batch 11",
        df_temp_median["batch 11", "WT"] - df_temp_median["batch 11", "KO44"],
        np.nan,
    ]
    i = i + 4


df_temp = threshold_result.melt(id_vars=["Threshold", "Batch"]).dropna()
df_temp["hue"] = df_temp.Batch + " - " + df_temp.variable
g = sns.lineplot(data=df_temp, x="Threshold", y="value", hue="hue")
g.set_title("Evolution of difference of medians along with threshold")
plt.show()

# https://stackoverflow.com/questions/14529838/apply-multiple-functions-to-multiple-groupby-columns Posiblemente hay una mejor manera de escribir este código basandose en ese post o el siguiente
# https://galea.medium.com/pandas-groupby-with-multiple-columns-c7a94ef0be86

# %%%% [md]
"""
En esta gráfica camos a ver si de verdad hay una diferencia en el porcentaje de tiempo que pasan cerca del borde
"""
