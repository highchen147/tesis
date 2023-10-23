 #include <iostream>
#include <string>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
using namespace std;

// Sobrecarga de * para multiplicación de escalar por un vector
vector<double> operator*(const vector<double>& v, double scalar) {
    vector<double> result(v.size());
    for (size_t i = 0; i < v.size(); i++) {
        result[i] = v[i] * scalar;
    }
    return result;
}

// Sobrecarga de + para suma de vectores 
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

 int main()
 {
    vector<double> v = {1.0, 0.0, 0.0};
    vector<double> u = {0.0, 1.0, 0.0};
    vector<double> w = {0.0, 0.0, 1.0};

    vector<double> escalares = {1.0, 2.0, 3.0};

    vector<vector<double>> V = {v, u, w};
    
    vector<double> resultado_suma = {0.0, 0.0, 0.0}; // Llamada al operador sobrecargado +
    
    for (int i = 0; i < 3; i++)
    {
        vector<double> vec = V[i]*escalares[i];
        resultado_suma = resultado_suma + vec;
    }
    for (double x : resultado_suma) {
        cout << x << " ";
    }
    cout << endl;

    


 }