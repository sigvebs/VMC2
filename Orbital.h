/* 
 * File:   Orbital.h
 * Author: zigg
 *
 * Created on May 10, 2012, 8:43 AM
 */

#ifndef ORBITAL_H
#define	ORBITAL_H

#include "Hermite.h"
#include <armadillo>
using namespace arma;
        
class Orbital {
public:
    Orbital();
    Orbital(int, double);
    Orbital(const Orbital& orig);
    virtual ~Orbital();
    
    virtual double evaluate(const rowvec &, const int, const int) = 0;
    virtual rowvec getGradient(const rowvec &, const int, const int) = 0;
    virtual double evaluateLaplacian(const rowvec &, const int, const int) = 0;
    virtual double evaluateExp(const rowvec &) = 0;
    virtual double variationalDerivative(const rowvec &, const int, const int) = 0;
    virtual void setNewAlpha(double alpha) = 0;
protected:
    int dim;
    Hermite *H;
    double alpha;
};

#endif	/* ORBITAL_H */

