import matplotlib.pyplot as plt
import numpy as np
# plt.rcParams['font.family'] = 'serif' 
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

# Crear datos para las rectas
x = np.linspace(-2.4, 2.4, 100)  # Valores de x de 0 a 10
x_0s = np.linspace(-1.8, 1.8, 8)
m_values = np.exp([-(x_)**2 for x_ in x_0s])  # Pendientes de 0 a 4

# Crear el plot
plt.figure(figsize=(8, 6))

for exp_m, m in zip(m_values, x_0s):
    # x = m * y + m  # Ecuación de la recta: y = mx
    y = (x-m)/exp_m
    plt.plot(x, y, color="purple")
    # plt.arrow(0,0,0,0.5, head_width=0.5, head_length=1, fc='BLUE', ec='blue')

# Configurar el aspecto del plot
# plt.title(r'Características de $u(x,0)=\exp(-x^{2})$')
plt.xlabel(r'$x$', fontsize=20)
plt.ylabel(r'$t$', fontsize=20)
plt.axhline(0, color='black',linewidth=0.2)
plt.axvline(0, color='black',linewidth=0.2)
plt.grid(color = 'gray', linestyle = '--', linewidth = 0.2, alpha=0.2)
plt.ylim(0, 1.4)
# plt.legend()
plt.savefig('caracteristicas_burgers.pdf', format='pdf')
# Mostrar el plot
plt.show()
