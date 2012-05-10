/* 
 * File:   Hermite.cpp
 * Author: zigg
 * 
 * Created on May 10, 2012, 9:04 AM
 */

#include "Hermite.h"

Hermite::Hermite() {
}

Hermite::Hermite(const Hermite& orig) {
}

Hermite::~Hermite() {
}

double Hermite::polynomial(const int degree, const double x) {
    double hermite_polynomial;

    if (degree == 0)
        hermite_polynomial = 1.0;
    else if (degree == 1)
        hermite_polynomial = 2.0 * x;
    else if (degree == 2)
        hermite_polynomial = 4.0 * x * x - 2;
    else if (degree == 3)
        hermite_polynomial = 8.0 * x * x * x - 12;
    else if (degree == 4)
        hermite_polynomial = 16.0 * x * x * x * x - 48.0 * x * x + 12;
    else
        hermite_polynomial = 0;
    
    return hermite_polynomial;
}
