import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter

def generar_subsets(data: pd.DataFrame) -> list[pd.DataFrame]:
    
    vacios = data[pd.isna(data["t"])]
    bounds = np.array([-1] + list(vacios.index))
    desfases = np.zeros(len(bounds), dtype=int)
    desfases[::2] = 1
    bounds = np.add(bounds, desfases)
    bounds = bounds[:-1]

    subsets = [data.iloc[bounds[i]:bounds[i+1], :] for i in range(len(bounds))[::2]]

    return subsets

def animacion(data: pd.DataFrame, subset: pd.DataFrame, cantidad: str, margen: float, ax: plt.Axes, line):
    
    # Obtener datos de x y u para el instante temporal actual
    x = subset['x'].values
    u = subset['u'].values
    
    # Actualizar los datos para la linea a trazar
    line.set_data(x, u)
    
    # Establecer límites del lienzo en donde se graficará
    min_u = data['u'].min()
    max_u = data['u'].max()
    min_y = min_u - (max_u-min_u)*margen
    max_y = max_u + (max_u-min_u)*margen
    min_x = data['x'].min()
    max_x = data['x'].max()
    ax.set_xlim([min_x, max_x])
    ax.set_ylim([min_y, max_y])
    # Obtener el tiempo
    t = subset['t'].values[0]
    # Establecer el título de la gráfica incluyendo el instante temporal
    ax.set_title("{}, t={} s".format(cantidad, t))

# importar datos
folder_data = "transmisiva980"
data_d = pd.read_csv("data/" + folder_data + "/densidad.dat", delimiter='\t', skip_blank_lines=False)
data_d.columns = ["t", "x", "u"]
data_p = pd.read_csv("data/" + folder_data + "/presion.dat", delimiter="\t", skip_blank_lines=False)
data_p.columns = ["t", "x", "u"]
data_u = pd.read_csv("data/" + folder_data + "/velocidad.dat", delimiter="\t", skip_blank_lines=False)
data_u.columns = ["t", "x", "u"]

# Intervalo entre frames en milisegundos
frames_sep = 50

# Crear la figura con tres sub-figuras
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4))
ax1: plt.Axes
ax2: plt.Axes
ax3: plt.Axes
# Primera animación
line1, = ax1.plot([], [], lw=0.8)
subsets_densidad = generar_subsets(data_d)

def init():
    line1.set_data([], [])
    return (line1,)
# Definir animacion de densidad
def animate_1(i):

    animacion(data_d, subsets_densidad[i], r'Densidad $\left(\frac{kg}{m^3}\right)$', 0.25, ax1, line1)

    return line1,
    
num_frames = len(subsets_densidad)
dt = (subsets_densidad[1]['t'].values[0] - subsets_densidad[0]['t'].values[0])*1000

anim1 = animation.FuncAnimation(fig, animate_1, init_func=init, frames=(range(num_frames)), repeat=False, interval=frames_sep)

# Segunda animación
line2, = ax2.plot([], [], lw=1)

subsets_presion = generar_subsets(data_p)

def init2():
    line2.set_data([], [])
    return (line2,)
# Definir animación de presión
def animate_2(i):
    animacion(data_p, subsets_presion[i], r'Presión $(Pa)$', 0.25, ax2, line2)
    return line2,


anim2 = animation.FuncAnimation(fig, animate_2, init_func=init2, frames=(range(num_frames)), repeat=False, interval=frames_sep)

# Tercera animación
line3, = ax3.plot([], [], lw=1)

subsets_velocidad = generar_subsets(data_u)

def init3():
    line3.set_data([], [])
    return (line3,)
# Definir animación de velocidad
def animate_3(i):
    animacion(data_u, subsets_velocidad[i], r'Velocidad $\left(\frac{m}{s}\right)$', 0.25, ax3, line3)
    return line3,

anim3 = animation.FuncAnimation(fig, animate_3, init_func=init3, frames=(range(num_frames)), repeat=False, interval=frames_sep)


plt.show()


