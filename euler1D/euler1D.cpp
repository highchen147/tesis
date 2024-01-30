/**
 * @file euler1D.cpp
 * @author Rodrigo Castillo (steverogersavengers@gmail.com)
 * @brief Programa integrador de las ecuaciones de Euler en una dimensión,
 * correspondiente a un gas ideal.
 * @version 0.3
 * @date 2023-05-04
 * 
 * @copyright Copyright (c) 2023
 * 
 * Para compilar: g++ euler1D.cpp funciones.cpp -o solver
 * 
 */
#include <iostream>
#include <string>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <sys/stat.h>
#include <algorithm>
#include "funciones.hpp" // Incluye funciones predefinidas
using namespace std;

int generateRandomNum();
double rho_inicial(double x, double L);
double p_inicial(double x, double L);
double u_inicial(double x, double L);
void calc_componentes_U(double *u1, double *u2, double *u3, double *rho, double *p, double *u, int Nx);
void calc_componentes_F(double *F1, double *F2, double *F3, double *rho, double *p, double *u, int Nx);
double rho_prom(double rho_L, double rho_R);
double u_prom(double u_L, double u_R, double rho_L, double rho_R);
double h_prom(double p_L, double p_R, double u_L, double u_R, double rho_L, double rho_R);
double c_prom(double p_L, double p_R, double rho_L, double rho_R, double u_L, double u_R);
void salida(ofstream &of, double *u, double *x, double tiempo, int N);
vector<double> flujo_euler(double rho, double p, double u);
vector<double> Flujo(vector<double> F_L, vector<double> F_R, double p_L, double p_R, double u_L, double u_R, double rho_L, double rho_R, bool entropy_fix);
vector<double> suma_p(double p_L, double p_R, double u_L, double u_R, double rho_L, double rho_R);
vector<double> operator+(const vector<double>& a, const vector<double>& b);
vector<double> operator-(const vector<double>& a, const vector<double>& b);
vector<double> operator*(const vector<double>& v, double scalar);
vector<double>& operator+=(vector<double>& v1, const vector<double>& v2);

const double Gamma = 1.4;

int main()
{
    // Inicializar el generador de números aleatorios
    srand(static_cast<unsigned>(time(nullptr)));
    int numeroAleatorio = generateRandomNum();

    // Parámetros temporales
    const double t_total = 3; // Tiempo total en segundos
    const double dt = 0.006; // Tamaño de paso temporal en segundos
    int Niter = floor(t_total/dt); // Número total de iteraciones
    const int num_outs = 5; // Número de gráficas de instantes temporales
    int out_cada = floor(Niter / num_outs); // Cada out_cada veces se 
                                            // imprimen los valores
    
    double tiempo = 0.0; // Variable de tiempo en la simulación

    // Parámetros espaciales
    int Nx = 500; // Número de celdas en el eje x
    double L = (10); // Largo del dominio en metros
    double dx = L/(Nx); // Tamaño de paso en el eje x

    // Otros parámetros
    bool correccion_de_entropia = false;

    // Archivos
    ofstream file_densidad;
    ofstream file_presion;
    ofstream file_velocidad;
    string nombreDirectorio;
    cout << "Ingrese el nombre del directorio que desea crear para almacenar los datos: ";
    cin >> nombreDirectorio;
    nombreDirectorio = "data/" + nombreDirectorio + to_string(numeroAleatorio);
    int directorio = mkdir(nombreDirectorio.c_str());
    cout << "Se ha creado " + nombreDirectorio << endl;
    string file_densidad_name = nombreDirectorio + "/densidad.dat";
    string file_presion_name = nombreDirectorio + "/presion.dat";
    string file_velocidad_name = nombreDirectorio + "/velocidad.dat";

    file_densidad.open(file_densidad_name, ios::out );
    file_presion.open(file_presion_name, ios::out);
    file_velocidad.open(file_velocidad_name, ios::out);
     
    // Arreglos
    // Cantidades físicas
    // Densidad
    double *rho = new double[Nx];
    double *rho_nueva = new double[Nx];
    // Velocidad
    double *u = new double[Nx]; 
    double *u_nueva = new double[Nx];
    // Presión
    double *p = new double[Nx]; 
    double *p_nueva = new double[Nx];
    // Componentes del vector U
    double *u1 = new double[Nx];
    double *u2 = new double[Nx];
    double *u3 = new double[Nx];
    // Componentes del vector F
    double *F1 = new double[Nx];
    double *F2 = new double[Nx];
    double *F3 = new double[Nx];
    // Celdas sobre el eje x
    double *x = new double[Nx];

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

    // Se declaran los vectores principales de la integración
    vector<double> U(3);
    vector<double> F(3);
    // Se calculan las componentes del vector U de acuerdo a su definición 
    calc_componentes_U(u1, u2, u3, rho, p, u, Nx);
    // Se calculan las componentes del vector F, que representa el flujo
    calc_componentes_F(F1, F2, F3, rho, p, u, Nx);

    // Se envían los datos iniciales
    salida(file_densidad, rho, x, tiempo, Nx);
    salida(file_presion, p, x, tiempo, Nx);
    salida(file_velocidad, u, x, tiempo, Nx);
    
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
            U_N=U-((Flujo(flujo_euler(rho[i], p[i], u[i]), 
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
        if ((k+1) % out_cada == 0)
        {
            salida(file_densidad, rho, x, tiempo, Nx);
            salida(file_presion, p, x, tiempo, Nx);
            salida(file_velocidad, u, x, tiempo, Nx);
            // cout << round(100*tiempo/t_total*100)/100 << endl;
            cout << tiempo <<endl;
        }
        // Actualizar el tiempo
        tiempo += dt;
    }
    std::cout << "Se ha creado " + nombreDirectorio << endl;
    
}

int generateRandomNum() {
    return rand() % 1000;
}

/**
 * @brief Función inicial de la velocidad
 * 
 * @param x Posición en x
 * @return double 
 */

double u_inicial(double x, double L)
{
    return step_neg(x, 1, 1, L/2);
}

/**
 * @brief Función inicial de la presión
 * 
 * @param x Posición en x
 * @return double 
 */
double p_inicial(double x, double L)
{
    double atm = (1.01325e5);
    // return 100*exp(-0.5*pow((x-L/2), 2));
    // return atm*1/12*(atan(x-L/2)+4.50);
    return step_neg(x, 0.5, 0.5, L/2);
}

/**
 * @brief Función inicial de la densidad
 * 
 * @param x Posición en x
 * @return double 
 */
double rho_inicial(double x, double L)
{
    // return step_neg(x, 1, 1, L/2);
    return 1+exp(-(x-L/2)*(x-L/2));
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
void calc_componentes_U(double *u1, double *u2, double *u3, double *rho, double *p, double *u, int Nx)
{
    for (int i = 0; i < Nx; i++)
    {
        u1[i] = rho[i];
        u2[i] = rho[i]*u[i];
        u3[i] = p[i]/(Gamma-1) + 0.5*rho[i]*pow(u[i], 2);
    }
    
}

/**
 * @brief Asignar valores a las componentes del vector F (flujo)
 * 
 * @param F1 Componente 1 de F
 * @param F2 Componente 2 de F
 * @param F3 Componente 3 de F
 * @param rho Densidad
 * @param p Presión
 * @param u Velocidad
 * @param Nx Tamaño de los arreglos que almacenan las funciones
 */
void calc_componentes_F(double *F1, double *F2, double *F3, double *rho, double *p, double *u, int Nx)
{
    for (int i = 0; i < Nx; i++)
    {
        F1[i] = rho[i]*u[i];
        F2[i] = p[i] + rho[i]*pow(u[i], 2);
        F3[i] = u[i]*(p[i]/(Gamma-1) + 0.5*rho[i]*pow(u[i], 2) + p[i]); 
    }
    
}

/**
 * @brief Calcula media geométrica entre densidades vecinas
 * 
 * @param rho_L Densidad a la izquierda
 * @param rho_R Densidad a la derecha
 * @return double 
 */
double rho_prom(double rho_L, double rho_R)
{
    return sqrt(rho_L*rho_R);
}

/**
 * @brief Calcula media ponderada de velocidades 
 * respecto a la raiz de densidades vecinas
 * 
 * @param u_L velocidad a la izquierda
 * @param u_R velocidad a la derecha
 * @param rho_L densidad a la izquierda
 * @param rho_R densidad a la derecha
 * @return double 
 */
double u_prom(double u_L, double u_R, double rho_L, double rho_R)
{
    return (sqrt(rho_L)*u_L+sqrt(rho_R)*u_R)/
    (sqrt(rho_L) + sqrt(rho_R));
}

/**
 * @brief Calcula media ponderada de la entalpía
 * respecto a la raíz de densidades vecinas
 * 
 * @param p_L presión a la izquierda
 * @param p_R presión a la derecha
 * @param u_L velocidad a la izquierda
 * @param u_R velocidad a la derecha
 * @param rho_L densidad a la izquierda
 * @param rho_R densidad a la derecha
 * @return double 
 */
double h_prom(double p_L, double p_R, double u_L, double u_R, double rho_L, double rho_R)
{
    // Se calcula la entalpía a la izquierda y a la derecha
    double h_L = Gamma/(Gamma-1.0)*(p_L/rho_L) + 0.5*u_L*u_L;
    double h_R = Gamma/(Gamma-1.0)*(p_R/rho_R) + 0.5*u_R*u_R;
    // Se devuelve la misma media ponderada que con u_prom pero ahora con las entalpías
    return (sqrt(rho_L)*h_L+sqrt(rho_R)*h_R)/(sqrt(rho_L) + sqrt(rho_R));
}

/**
 * @brief Calcula la media ponderada de a, correspondiente a velocidad del sonido en casos particulares.
 * También respecto a las densidades vecinas.
 * 
 * @param p_L 
 * @param p_R 
 * @param rho_L 
 * @param rho_R 
 * @return double 
 */
double c_prom(double p_L, double p_R, double rho_L, double rho_R, double u_L, double u_R)
{
    double h = h_prom(p_L, p_R, u_L, u_R, rho_L, rho_R);
    double u = u_prom(u_L,u_R, rho_L, rho_R);
    return sqrt((Gamma - 1)*(h - 0.5*pow(u,2)));
}

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
vector<double> suma_p(double p_L, double p_R, double u_L, double u_R, double rho_L, double rho_R)
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
    // vector de vectores
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
    vector<double> lambda_L = {u_L-c_L, u, u_L+c_L};
    vector<double> lambda_R = {u_R-c_R, u, u_R+c_R};
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

/**
 * @brief Calcula el flujo F de la ecuación de Euler en forma conservativa.
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

/**
 * @brief Función que envía datos a los archivos, con instantes
 * temporales separados por doble enter
 * 
 * @param of Archivo de datos
 * @param u Arreglo que almacena los valores de las funciones
 * @param x Arreglo de dimensión espacial
 * @param tiempo Instante temporal en cuestión
 * @param N Tamaño del arreglo
 */
void salida(ofstream &of, double *u, double *x, double tiempo, int N){
    for (int i = 0; i < N; i++)
    {
        of << tiempo << "\t" << x[i] << "\t" << u[i] << endl;
    }
    of << endl << endl;
}

vector<double> operator+(const vector<double>& a, const vector<double>& b) {
    // Asegurarse de que ambos vectores tengan el mismo tamaño
    if (a.size() != b.size()) {
        throw runtime_error("Los vectores deben tener el mismo tamaño.");
    }
    
    vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); i++) {
        result[i] = a[i] + b[i];
    }
    return result;
}

vector<double> operator-(const vector<double>& a, const vector<double>& b) {
    // Asegurarse de que ambos vectores tengan el mismo tamaño
    if (a.size() != b.size()) {
        throw runtime_error("Los vectores deben tener el mismo tamaño.");
    }
    
    vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); i++) {
        result[i] = a[i] - b[i];
    }
    return result;
}

vector<double> operator*(const vector<double>& v, double scalar) {
    vector<double> result(v.size());
    for (size_t i = 0; i < v.size(); i++) {
        result[i] = v[i] * scalar;
    }
    return result;
}

vector<double>& operator+=(vector<double>& v1, const vector<double>& v2) {
    if (v1.size() != v2.size()) {
        throw invalid_argument("Los vectores deben tener el mismo tamaño.");
    }

    for (size_t i = 0; i < v1.size(); i++) {
        v1[i] += v2[i];
    }

    return v1;
}