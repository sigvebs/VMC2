/* 
 * File:   OneBodyDensity.h
 * Author: sigve
 *
 * Created on 11. mai 2012, 17:01
 */

#ifndef ONEBODYDENSITY_H
#define	ONEBODYDENSITY_H

#include "WaveFunction.h"

class OneBodyDensity {
public:
    OneBodyDensity();
    OneBodyDensity(const OneBodyDensity& orig);
    virtual ~OneBodyDensity();
    
    void writeToFile(const mat &, double);
    double McIntegrator( const vec &, int );
    mat normalize(const mat &);
protected:
    int McSamples;
    WaveFunction* wf;
    int dim;
    int nParticles;
    long idum; 
    int myRank;
    int nNodes;
    
    std::string fileName;
};

#endif	/* ONEBODYDENSITY_H */

