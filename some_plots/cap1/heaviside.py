import numpy as np
import matplotlib.pyplot as plt

def heaviside(x):
    return np.piecewise(x, [x < 50, x >= 50], [40, 10])

# def heaviside

# Crear un rango de valores x
x = np.linspace(0, 100, 1000)

# Calcular los valores de la función paso (Heaviside)
y = heaviside(x)
y_2 = heaviside(x-25)

# Crear el gráfico
plt.plot(x,y_2, color="#7cb8fc")
plt.plot(x, y)

plt.yticks([0,10,40], ["0", "$u_R$", "$u_L$"], )

# Marcar la discontinuidad verticalmente con líneas punteadas
# plt.axvline(0, ymin=0.5, ymax=1)

# Configurar el gráfico
# plt.title(r'$t=0$')
plt.xlim(0,100)
plt.ylim(0,50)
plt.xlabel(r'$x$', fontsize=23)
plt.ylabel(r'$u$', fontsize=23)
plt.tick_params(axis="y", labelsize=15)
plt.grid(True)
plt.show()
