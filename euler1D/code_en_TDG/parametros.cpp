const double Gamma = 1.4;
// Parámetros temporales
const double t_total = 4; // Tiempo total en segundos
const double dt = 0.005; // Tamaño de paso temporal en segundos
int Niter = floor(t_total/dt); // Número total de iteraciones
const int num_outs = 400; // Número de instantes temporales producidos
int out_cada = floor(Niter / num_outs); // Cada out_cada veces se 
                                        // imprimen los valores
double tiempo = 0.0; // Variable de tiempo en la simulación

// Parámetros espaciales
int Nx = 500; // Número de celdas en el eje x
double L = (10); // Largo del dominio en metros
double dx = L/(Nx); // Tamaño de paso en el eje x

// Otros parámetros
bool correccion_de_entropia = true;