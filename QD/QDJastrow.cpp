/* 
 * File:   QDJastrow.cpp
 * Author: zigg
 * 
 * Created on May 10, 2012, 9:33 AM
 */

#include "QDJastrow.h"

////////////////////////////////////////////////////////////////////////////////
QDJastrow::QDJastrow() {
}

////////////////////////////////////////////////////////////////////////////////
QDJastrow::QDJastrow(const QDJastrow& orig) {
}

////////////////////////////////////////////////////////////////////////////////
QDJastrow::~QDJastrow() {
}

////////////////////////////////////////////////////////////////////////////////
QDJastrow::QDJastrow(int dim, int nParticles, double beta) : Jastrow(dim, nParticles, beta) {

    // Initiating a matrix with all the spin dependant a-values.
    a = zeros(nParticles, nParticles);

    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nParticles; j++) {

            if (i == j)
                a(i, j) = 0;
            else if (i < nParticles / 2 && j >= nParticles / 2 || i >= nParticles / 2 && j < nParticles / 2)
                a(i, j) = 1.0;
            else
                a(i, j) = 1.0 / 3.0;
        }
    }
    
#if 0
    cout
            << "\t dim = " << this->dim 
            << "\t nParticles = " << this->nParticles 
            << "\t beta = " << this->beta
            << "\t a= " << a
            << endl;
#endif
     
}

////////////////////////////////////////////////////////////////////////////////
double QDJastrow::evaluate(const mat &r) {
    double r_norm;
    double value = 0;

    // ------------ Move the loop into the WF ------------
    for (int i = 0; i < nParticles; i++) {
        for (int j = i + 1; j < nParticles; j++) {

            r_norm = 0;
            for (int d = 0; d < dim; d++) {
                r_norm += (r(i, d) - r(j, d))*(r(i, d) - r(j, d));
            }
            r_norm = sqrt(r_norm);
            value += a(i, j) * r_norm / ((1 + beta * r_norm));
        }
    }

    return value;
}


////////////////////////////////////////////////////////////////////////////////
double QDJastrow::getLaplacian(const mat &r, int i) {
    double r_ki;
    double sum = 0;

    // Before i
    for (int k = 0; k < i; k++) {
        r_ki = norm(r.row(k) - r.row(i), 2);
        sum += a(k, i)*(1 - beta * r_ki) / r_ki / pow((1 + beta * r_ki), 3);
    }

    // After i
    for (int k = i + 1; k < nParticles; k++) {
        r_ki = norm(r.row(k) - r.row(i), 2);
        sum += a(k, i)*(1 - beta * r_ki) / r_ki / pow((1 + beta * r_ki), 3);
    }

    double Gj_sq = dot(gradient, gradient);
    return Gj_sq + sum;
}

////////////////////////////////////////////////////////////////////////////////
void QDJastrow::computeGradient(const mat &r, int i) {
    double r_ki;
    gradient = zeros(1, dim);

    // Before i
    for (int k = 0; k < i; k++) {
        r_ki = norm(r.row(k) - r.row(i), 2);
        gradient += a(k, i) / r_ki / pow((1 + beta * r_ki), 2)*(r.row(i) - r.row(k));
    }

    // After i
    for (int k = i + 1; k < nParticles; k++) {
        r_ki = norm(r.row(k) - r.row(i), 2);
        gradient += a(k, i) / r_ki / pow((1 + beta * r_ki), 2)*(r.row(i) - r.row(k));
    }
}

////////////////////////////////////////////////////////////////////////////////
double QDJastrow::getVariationalGradient(const mat &r) {
    double r_sq, r_norm;
    double value = 1;

    for (int i = 0; i < nParticles; i++) {
        for (int j = i + 1; j < nParticles; j++) {

            r_sq = 0;
            for (int d = 0; d < dim; d++)
                r_sq += (r(i, d) - r(j, d))*(r(i, d) - r(j, d));

            r_norm = sqrt(r_sq);
            value += a(i, j) * r_sq / ((1 + beta * r_norm)*(1 + beta * r_norm));
        }
    }
    
    return -value;
}