/* 
 * File:   SGD.h
 * Author: Sigve
 *
 * Created on May 10, 2012, 7:42 PM
 */

#ifndef SGD_H
#define	SGD_H

#include <armadillo>
using namespace arma;

class SGD {
public:
    SGD();
    SGD(const SGD& orig);
    virtual ~SGD();
private:
    double E;
    double Esq;
    double accepted;    
    
    int McSamples;
    int SGDSamples;
    int m;
    
    bool importanceSampling;
    bool usingJastrow;
    int thermalization;
    int nParticles;    
    double alpha, beta, w;
    int dim;
    
    bool writeToFile;
    std::string fileName;
    
    // Minimization parameters.
    double maxStep;
    double A, a, expo;
    double f, fMax, fMin, omega;
};

#endif	/* SGD_H */

