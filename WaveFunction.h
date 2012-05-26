/* 
 * File:   WaveFunction.h
 * Author: Sigve
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
    rowvec getVariationGradient();
    void setNewVariationalParameters(double, double);

    // Brute force
    bool tryNewPositionBF(int);
    void setOptimalStepLength();
    double difference(double);

    // DMC
    bool DMCtryNewPosition(int);

    virtual WaveFunction* clone() const {
        return new WaveFunction(*this);
    }

    void test();    
    mat getROld(){
        return rOld;
    }
private:
    bool usingJastrow;
    int activeParticle;
    long idum;

    mat rOld;
    mat rNew;
    
    mat qForce;
    mat qForceOld;    
    
    rowvec gradientSlater;
    rowvec gradientJastrow;

    Jastrow *jastrow;
    Slater *slater;
    Orbital *orbital;
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

