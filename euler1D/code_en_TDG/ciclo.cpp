// Comienza a correr el tiempo
tiempo += dt;
// Ciclo principal
for (int k = 0; k < Niter; k++)
{
// Se calculan las componentes del vector U
calc_componentes_U(u1, u2, u3, rho, p, u, Nx);
// Ciclo para integración espacial
for (int i = 1; i < Nx-1; i++)
{
    // Definición del vector U_N que corresponde al vector U en
    // el siguiente instante de tiempo
    vector<double> U_N(3);
    // Definir valores de U
    U = {u1[i], u2[i], u3[i]};
    // Actualizar e integrar U
    U_N = U-((Flujo(flujo_euler(rho[i], p[i], u[i]), 
                    flujo_euler(rho[i+1], p[i+1], u[i+1]),
                    p[i], p[i+1],
                    u[i], u[i+1],
                    rho[i], rho[i+1], correccion_de_entropia)-
              Flujo(flujo_euler(rho[i-1], p[i-1], u[i-1]), 
                    flujo_euler(rho[i], p[i], u[i]), 
                    p[i-1], p[i],
                    u[i-1], u[i],
                    rho[i-1], rho[i], correccion_de_entropia))*
                    (dt/dx));

    // Despejar variables físicas de U
    rho_nueva[i] = U_N[0];
    u_nueva[i] = U_N[1]/rho_nueva[i];
    p_nueva[i] = (U_N[2]-0.5*rho_nueva[i]*pow(u_nueva[i], 2))*
                    (Gamma-1);

    // En caso la velocidad del sonido sea nula, se envía una
    // alerta, indicando el momento y la coordenada.
    if (c_prom(p[i],p[i+1],rho[i],rho[i+1],u[i],u[i+1])==0.0)
    {
        cout << "tiempo: " << k << endl;
        cout << "coordenada : " << i << endl;
    }
}

// Actualizar variables físicas
for (int i = 1; i < Nx-1; i++)
{
    rho[i] = rho_nueva[i];
    u[i] = u_nueva[i];
    p[i] = p_nueva[i];
}

// Condiciones frontera transmisivas
rho[0] = rho[1];
rho[Nx-1] = rho[Nx-2];
u[0] = u[1];
u[Nx-1] = u[Nx-2];
p[0] = p[1];
p[Nx-1] = p[Nx-2];

// Se evalúa si la iteración corresponde a un instante de
// impresión de datos
if (k % out_cada == 0)
{
    salida(file_densidad, rho, x, tiempo, Nx);
    salida(file_presion, p, x, tiempo, Nx);
    salida(file_velocidad, u, x, tiempo, Nx);
    cout << round(100*tiempo/t_total*100)/100 << endl;
}
// Actualizar el tiempo
tiempo += dt;
}