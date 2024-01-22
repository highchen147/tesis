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

// Se declaran los vectores principales de la integración
vector<double> U(3);
vector<double> F(3);