/**
 * @brief Calcula el flujo F de la ecuación de Euler en forma
 * conservativa.
 * 
 * @param rho Densidad
 * @param p Presión
 * @param u Velocidad
 * @return vector<double> 
 */
vector<double> flujo_euler(double rho, double p, double u)
{
    vector<double> f_resultante(3);
    double F1 = rho*u;
    double F2 = p + rho*pow(u, 2);
    double F3 = u*(p/(Gamma-1) + 0.5*rho*pow(u, 2) + p);
    f_resultante = {F1, F2, F3};
    return f_resultante;
}

/**
 * @brief Calcula el flujo entre celdas utilizando el esquema de Roe
 * 
 * @param F_L Flujo exacto en la celda izquierda
 * @param F_R Flujo exacto en la celda derecha
 * @param p_L Presión a la izquierda
 * @param p_R Presión a la derecha
 * @param u_L Velocidad a la izquierda
 * @param u_R Velocidad a la derecha
 * @param rho_L Densidad a la izquierda
 * @param rho_R Densidad a la derecha
 * @param entropy_fix Parámetro para decidir el uso de la
 * correción de entropía
 * @return vector<double> 
 */
vector<double> Flujo(
    vector<double> F_L, 
    vector<double> F_R, 
    double p_L, 
    double p_R, 
    double u_L, 
    double u_R, 
    double rho_L, 
    double rho_R,
    bool entropy_fix)
{
    vector<double> F_prom = (F_L + F_R)*0.5;
    if (entropy_fix)
    {
        return F_prom - (suma_p_fix(p_L,
                                    p_R, 
                                    u_L, 
                                    u_R, 
                                    rho_L, 
                                    rho_R)*0.5);
    }
    else 
    {
        return F_prom - (suma_p(p_L, 
                                p_R, 
                                u_L, 
                                u_R, 
                                rho_L, 
                                rho_R)*0.5);
    }
}