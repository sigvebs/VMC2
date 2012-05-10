/* 
 * File:   Jastrow.cpp
 * Author: zigg
 * 
 * Created on May 10, 2012, 9:30 AM
 */

#include "Jastrow.h"

////////////////////////////////////////////////////////////////////////////////
Jastrow::Jastrow() {
}

////////////////////////////////////////////////////////////////////////////////
Jastrow::Jastrow(const Jastrow& orig) {
}

////////////////////////////////////////////////////////////////////////////////
Jastrow::~Jastrow() {
}

////////////////////////////////////////////////////////////////////////////////
Jastrow::Jastrow(int dim, int nParticles, double beta) {
    this->dim = dim;
    this->nParticles = nParticles;
    this->beta = beta;
}

////////////////////////////////////////////////////////////////////////////////
double Jastrow::getRatio(const mat &r_new, const mat &r_old ) {
    // This is an inefficent way of solving the problem since 
    // the are evaluating all new distances. Which particle 
    // we are evaluating should be implemented in the WF. 
    double delta_J = evaluate(r_new) - evaluate(r_old);
    return exp(delta_J);
}