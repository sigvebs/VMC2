/* 
 * File:   Jastrow.h
 * Author: zigg
 *
 * Created on May 10, 2012, 9:30 AM
 */

#ifndef JASTROW_H
#define	JASTROW_H

#include <armadillo>
using namespace arma;

class Jastrow {
public:
    Jastrow();
    Jastrow(int, int, double);
    Jastrow(const Jastrow& orig);
    virtual ~Jastrow();

    virtual double getVariationalGradient(const mat &) = 0;
    virtual double evaluate(const mat &) = 0;
    virtual void computeGradient(const mat &, int) = 0;
    virtual double getLaplacian(const mat &, int) = 0;
    double getRatio(const mat &, const mat &);
    
    virtual void setNewBeta(double) = 0;

    rowvec getGradient() {
        return gradient;
    };
protected:
    int dim;
    int nParticles;
    double beta;
    rowvec gradient;
};

#endif	/* JASTROW_H */

