import pandas as pd
import numpy as np
from math import floor, ceil
from matplotlib import pyplot as plt

datos_ = pd.read_csv("sol-burg-vis1DDF-nu-0.50.dat", sep="\t", names=["t", "x", "u"])

archivos = {"0.5": "sol-burg-vis1DDF-nu-0.50.dat",
            "1.6": "sol-burg-vis1DDF-nu-1.60.dat",
            "3.0": "sol-burg-vis1DDF-nu-3.00.dat"}
datos = {}

def asignar_instantes(data: pd.DataFrame) -> list[float]:
    data["Instante"] = data["t"].rank(method='dense').astype(int)
    return data

for archivo in archivos:
    datos[archivo] = pd.read_csv(f"sol-burg-vis1DDF-nu-{archivo}0.dat", sep="\t", names=["t", "x", "u"])
    datos[archivo] = asignar_instantes(datos[archivo])

datos["0"] = asignar_instantes(pd.read_csv("gauss-fija.dat", sep="\t", names=["t", "x", "u"]))

def graficar(dic_datos: dict[pd.DataFrame], instantes_temporales: list, columna="u", eje="x", imprimir=True):
    # Configurar latex
    if imprimir:
        plt.rcParams['font.family'] = 'serif' 
        plt.rcParams['text.usetex'] = True
        plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
    
    max_y = max([frame["u"].max() for frame in dic_datos.values()])
    
    for instante in instantes_temporales:
        plt.clf()
        plt.figure(figsize=(8, 6))
        for nombre in dic_datos:
            data: pd.DataFrame
            data = dic_datos[nombre].dropna()
            data = data[data["Instante"] == (instante)]
            x = data[eje].values[:]
            y = data[columna].values[:]
            label = r"$\varepsilon = "
            label += nombre
            label += r"$"
            plt.plot(x, y, linewidth = 1, label=label)
            

        x_label = "$" + eje + "$"
        y_label = "$" + columna + "$"
        plt.xlabel(x_label, fontsize=23)
        plt.ylabel(y_label, fontsize=23)
        plt.tick_params(axis="both", labelsize=15)
        plt.axhline(0, color='black',linewidth=0.2)
        plt.grid(color = 'gray', linestyle = '--', linewidth = 0.2, alpha=0.2)
        plt.xlim(10,90)
        plt.ylim(0,max_y*1.25)
        plt.legend(fontsize = 20)
        if imprimir:
            plt.savefig(f'graficas/viscosidades-{instante}.pdf', format='pdf')
        else:
            plt.show()


graficar(dic_datos=datos,
         instantes_temporales = [1, 100, 200, 400], 
         columna="u", 
         eje="x", 
         imprimir=False)
