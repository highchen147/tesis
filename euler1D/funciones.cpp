#include "funciones.hpp"

double step_pos(double x, double max, double min, double x_0) {
    if (x <= x_0)
    {
        return min;
    } 
    else 
    {
        return max;
    }
    
}

double step_neg(double x, double max, double min, double x_0) {
    if (x <= x_0)
    {
        return max;
    } 
    else 
    {
        return min;
    }
    
}