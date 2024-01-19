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

double step_neg(double x, double max, double min, double x_0)
{
    if (x <= x_0)
    {
        return max;
    } 
    else 
    {
        return min;
    } 
}