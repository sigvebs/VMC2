/* 
 * File:   ComputeGrid.h
 * Author: Sigve
 *
 * Created on May 10, 2012, 6:30 PM
 */

#ifndef COMPUTEGRID_H
#define	COMPUTEGRID_H

#include <armadillo>
using namespace arma;

class ComputeGrid {
public:
    ComputeGrid();
    ComputeGrid(const ComputeGrid& orig);
    virtual ~ComputeGrid();
    
    void writeEnergyToFile();
private:
    mat E;
    mat Esq;
    double accepted;    
    std::string fileName;
    
    int McSamples;
    bool importanceSampling;
    int thermalization;
    int nParticles;    
    double alpha, beta, w;
    double alphaStart, betaStart;
    double deltaAlpha, deltaBeta;
    int nAlpha, nBeta;
    bool usingJastrow;
    int dim;
};

#endif	/* COMPUTEGRID_H */

