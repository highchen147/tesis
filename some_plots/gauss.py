import matplotlib.pyplot as plt
import numpy as np
# plt.rcParams['font.family'] = 'serif' 
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
plt.rcParams['font.family'] = 'serif'

# Crear datos para las rectas
x = np.linspace(-2.4, 2.4, 100)  # Valores de x de 0 a 10

plt.plot(x, np.exp(-np.power(x,2)), color="purple")
plt.title(r"Gráfica de $u(x,0)$")
plt.xlabel(r'$x$', fontsize=20)
plt.ylabel(r'$u$', fontsize=20)
plt.axhline(0, color='black',linewidth=0.2)
plt.axvline(0, color='black',linewidth=0.2)
plt.grid(color = 'gray', linestyle = '--', linewidth = 0.2, alpha=0.7)

plt.xlim(-2.4, 2.4)
plt.savefig('cond_inicial.pdf', format='pdf')
# plt.show()