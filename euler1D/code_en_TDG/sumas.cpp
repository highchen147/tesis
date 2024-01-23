/**
 * @brief Suma sobre p de los autovectores de la matriz A del sistema
 * por sus fuerzas y autovalores.
 * 
 * @param p_L Presión a la izquierda
 * @param p_R Presión a la derecha
 * @param u_L Velocidad a la izquierda
 * @param u_R Velocidad a la derecha
 * @param rho_L Densidad a la izquierda
 * @param rho_R Densidad a la derecha
 * @return vector<double> 
 */
vector<double> suma_p(double p_L, double p_R, 
                      double u_L, double u_R, 
                      double rho_L, double rho_R)
{
    // Cálculo de promedios de Roe
    // Velocidad promedio
    double u = u_prom(u_L, u_R, rho_L, rho_R);
    // Densidad promedio
    double rho = rho_prom(rho_L, rho_R);
    // Entalpía promedio
    double h = h_prom(p_L, p_R, u_L, u_R, rho_L, rho_R);
    // Velocidad del sonido promedio
    double c = c_prom(p_L, p_R, rho_L, rho_R, u_L, u_R);
    // Cálculo de diferencias laterales
    double dp = p_R - p_L;
    double du = u_R - u_L;
    double drho = rho_R - rho_L;
    // Coeficientes alfa
    double alfa_1 = 0.5*(dp-rho*c*du)/pow(c, 2);
    double alfa_2 = (pow(c, 2)*drho-dp)/pow(c, 2);
    double alfa_3 = 0.5*(dp+rho*c*du)/pow(c, 2);
    // Vector de coeficientes alfa
    vector<double> alfa = {alfa_1, alfa_2, alfa_3};
    // Vector de autovalores de las ondas centrales
    vector<double> lambda = {u-c, u, u+c};
    // Autovectores
    vector<double> r_1 = {1, u-c, h-u*c};
    vector<double> r_2 = {1, u, 0.5*pow(u,2)};
    vector<double> r_3 = {1, u+c, h+u*c};
    // vector de vectores propios
    vector<vector<double>> r_vec = {r_1, r_2, r_3};

    // Se declara el vector resultante de la suma, 
    // de dimensión 3 y con ceros.
    vector<double> resultado(3, 0);

    // Se realiza la suma 
    for (int i = 0; i < 3; i++)
    {
        // Se define la variable que almacena el valor absoluto de 
        // cada autovalor
        double lambda_i = abs(lambda[i]);
        resultado += r_vec[i]*(alfa[i]*lambda_i);
    }
    return resultado;
}

/**
 * @brief Suma sobre p de los autovectores de la matriz A del 
 * sistema por sus fuerzas y autovalores considerando la 
 * corrección HH de entropía.
 * 
 * @param p_L Presión a la izquierda
 * @param p_R Presión a la derecha
 * @param u_L Velocidad a la izquierda
 * @param u_R Velocidad a la derecha
 * @param rho_L Densidad a la izquierda
 * @param rho_R Densidad a la derecha
 * @return vector<double> 
 */
vector<double> suma_p_fix(double p_L, double p_R, double u_L, double u_R, double rho_L, double rho_R)
{
    // Cálculo de promedios de Roe
    double u = u_prom(u_L, u_R, rho_L, rho_R);
    double rho = rho_prom(rho_L, rho_R);
    double h = h_prom(p_L, p_R, u_L, u_R, rho_L, rho_R);
    double c = c_prom(p_L, p_R, rho_L, rho_R, u_L, u_R);
    // Cálculo de diferencias
    double dp = p_R - p_L;
    double du = u_R - u_L;
    double drho = rho_R - rho_L;
    // Coeficientes alfa
    double alfa_1 = 0.5*(dp-rho*c*du)/pow(c, 2);
    double alfa_2 = (pow(c, 2)*drho-dp)/pow(c, 2);
    double alfa_3 = 0.5*(dp+rho*c*du)/pow(c, 2);
    vector<double> alfa = {alfa_1, alfa_2, alfa_3};
    // Autovalores de las ondas centrales
    vector<double> lambda = {u-c, u, u+c};
    // Construcción de autovalores de las ondas laterales
    double c_L = c_prom(p_L, p_L, rho_L, rho_L, u_L, u_L);
    double c_R = c_prom(p_R, p_R, rho_R, rho_R, u_R, u_R);
    vector<double> lambda_L = {u_L - c_L, u, u_L + c_L};
    vector<double> lambda_R = {u_R - c_R, u, u_R + c_R};
    // Autovectores
    vector<double> r_1 = {1, u-c, h-u*c};
    vector<double> r_2 = {1, u, 0.5*pow(u,2)};
    vector<double> r_3 = {1, u+c, h+u*c};
    // vector de vectores
    vector<vector<double>> r_vec = {r_1, r_2, r_3};

    // Se declara el vector resultante de la suma, 
    // de dimensión 3 y con ceros.
    vector<double> resultado(3, 0);

    // Se realiza la suma 
    for (int i = 0; i < 3; i++)
    {
        // Se define la variable que almacena el valor absoluto del 
        // autovalor iésimo.
        double lambda_i = abs(lambda[i]);
        // Se define delta de la corrección HH.
        double delta_i = max(lambda_i-lambda_L[i], 
                             lambda_R[i]-lambda_i);
        delta_i = max(delta_i, 0.0);
        // Aplicación de la función por partes de la correción HH.
        if (lambda_i < delta_i)
        {
            lambda_i = delta_i;
        }
        else
        {
            lambda_i = abs(lambda[i]);
        }
        resultado += r_vec[i]*(alfa[i]*lambda_i);
    }
    return resultado;
}