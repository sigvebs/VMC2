/* 
 * File:   DMC.h
 * Author: Sigve
 *
 * Created on May 12, 2012, 1:33 PM
 */

#ifndef DMC_H
#define	DMC_H

#include "WaveFunction.h"

#include <cstdlib>
#include <list>
#include <vector>

class DMC {
public:
    DMC();
    DMC(const DMC& orig);
    virtual ~DMC();

    void writeDistributionToFile(std::string);
    void writeDistributionToFileVMC(WaveFunction*, std::string);
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
    int DMCSamples;

    bool importanceSampling;
    bool usingJastrow;
    int thermalization;
    int nParticles;
    double alpha, beta, w;
    int dim;
    double tau;

    std::vector<WaveFunction*> walkers;
};

#endif	/* DMC_H */

