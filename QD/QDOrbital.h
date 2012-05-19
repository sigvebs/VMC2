/* 
 * File:   QDOrbital.h
 * Author: zigg
 *
 * Created on May 10, 2012, 8:48 AM
 */

#ifndef QDORBITAL_H
#define	QDORBITAL_H

#include "../Orbital.h"

class QDOrbital: public Orbital {
public:
    QDOrbital();
    QDOrbital(int, double, double);
    QDOrbital(const QDOrbital& orig);
    virtual ~QDOrbital();    
    
    virtual double evaluate(const rowvec &, const int, const int);
    virtual rowvec getGradient(const rowvec &, const int, const int);
    virtual double evaluateLaplacian(const rowvec &, const int, const int);
    virtual double evaluateExp(const rowvec &);
    virtual double variationalDerivative(const rowvec &, const int, const int);
    virtual void setNewAlpha(double);
private:
    double w;
    double sqrtW;
    double wAlpha;
    double sqrtWAlpha;
    double sqrtWDividedAlpha;
    
    rowvec gradient;
};

#endif	/* QDORBITAL_H */

