/* 
 * File:   DMC.h
 * Author: Sigve
 *
 * Created on May 12, 2012, 1:33 PM
 */

#ifndef DMC_H
#define	DMC_H

#include "WaveFunction.h"

class DMC {
public:
    DMC();
    DMC(const DMC& orig);
    virtual ~DMC();

    WaveFunction* thermalizedWalker();
private:
    int myRank, nNodes;
    
    double E;
    double ET;
    double Esq;
    double accepted;

    int McSamples;
    int nWalkers;
    int nSteps;

    bool importanceSampling;
    bool usingJastrow;
    int thermalization;
    int nParticles;
    double alpha, beta, w;
    int dim;
    double tau;
};

#endif	/* DMC_H */

