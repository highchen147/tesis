// Inicialización de arreglos y condiciones iniciales
// Dominio espacial
for (int i = 0; i < Nx; i++)
{
    x[i] = i*dx;
}
// Cantidades físicas
for (int i = 0; i < Nx; i++)
{
    // Densidad
    rho[i] = rho_inicial(x[i], L);
    // Presión
    p[i] = p_inicial(x[i], L);
    // Velocidad
    u[i] = u_inicial(x[i], L);
}

double u_inicial(double x, double L)
{
    return step_neg(x, -1, 1, L/2);
}

/**
 * @brief Asignar valores a las componentes del vector U
 * 
 * @param u1 Componente 1 de U
 * @param u2 Componente 2 de U
 * @param u3 Componente 3 de U
 * @param rho Densidad
 * @param p Presión
 * @param u Velocidad
 * @param Nx Tamaño de los arreglos que almacenan las funciones
 */
void calc_componentes_U(double *u1, 
                        double *u2, 
                        double *u3, 
                        double *rho, 
                        double *p, 
                        double *u, 
                        int Nx)
{
    for (int i = 0; i < Nx; i++)
    {
        u1[i] = rho[i];
        u2[i] = rho[i]*u[i];
        u3[i] = p[i]/(Gamma-1) + 0.5*rho[i]*pow(u[i], 2);
    }
}