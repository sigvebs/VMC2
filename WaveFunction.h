/* 
 * File:   WaveFunction.h
 * Author: zigg
 *
 * Created on May 9, 2012, 11:05 PM
 */

#ifndef WAVEFUNCTION_H
#define	WAVEFUNCTION_H

#include <armadillo>
using namespace arma;

#include "Slater.h"
#include "Orbital.h"
#include "Jastrow.h"
#include "Hamiltonian.h"

class WaveFunction {
public:
    WaveFunction();
    WaveFunction(int, int, long, Orbital *, Jastrow*, Hamiltonian*);
    WaveFunction(const WaveFunction& orig);
    virtual ~WaveFunction();

    bool tryNewPosition(int);
    double sampleEnergy();
 
    double WFRatio();
    double evaluate(const mat &);
    mat newQForce();
    void setNewR(const mat &, int);
    void evaluateNew();
    void acceptMove();
    void calculateEnergy();
    void initSlater();
    
    // Brute force
    bool tryNewPositionBF(int);
    void setOptimalStepLength();
    double difference(double);

private:
    bool usingJastrow;
    int activeParticle;
    long idum;

    mat rOld;
    mat rNew;
    
    mat qForce;
    mat qForceOld;
    
    Jastrow *jastrow;
    Slater *slater;
    Hamiltonian *hamiltonian;

    double E;
    int dim;
    int nParticles;
    
    // Importance sampling
    double dt;
    double sqrtDt;
    double D;
    
    // Brute force
    double stepLength;
};

#endif	/* WAVEFUNCTION_H */

