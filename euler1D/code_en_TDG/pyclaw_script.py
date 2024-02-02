# @title Euler 1D set 1
nombre = "set1"
from clawpack import pyclaw
from clawpack import riemann
from funciones import funcion_paso
# Importar el solucionador
solver = pyclaw.ClawSolver1D(riemann.euler_1D_py.euler_roe_1D)
# Establecer el kernel del algoritmo de integración
solver.kernel_language = "Python"
# Condiciones de frontera.
solver.bc_upper[0] = pyclaw.BC.extrap
solver.bc_lower[0] = pyclaw.BC.extrap
# Dominio
domain = pyclaw.Domain([0], [10], [500])
# Objeto de solución
solution = pyclaw.Solution(solver.num_eqn, domain)
# Objeto de estado de la simulación
state = solution.state
# Centros de celdas de dominio espacial
xc = state.grid.p_centers[0]
from numpy import exp
# Gamma
cda_gamma = 1.4
# Definicion de funciones auxiliares de
# variables independientes y conservadas
densidad = 1 + exp(-((xc-0.01)-5)**2)
velocidad = 1.0
presion = 0.5
energia = presion/(cda_gamma - 1) + 0.5*densidad*velocidad**2
# Asignación de variables conservadas
state.q[0,:] = densidad
state.q[1,:] = densidad*velocidad
state.q[2,:] = energia
# Parámetros del problema
state.problem_data["gamma"] = cda_gamma
state.problem_data["gamma1"] = cda_gamma - 1
state.problem_data["efix"] = False
# Controlador
controller = pyclaw.Controller()
controller.solution = solution
controller.solver = solver
controller.num_output_times = 5
controller.tfinal = 3
controller.output_format = "ascii"
controller.outdir = nombre
# Correr la simulación
status = controller.run()