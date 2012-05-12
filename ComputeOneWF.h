/* 
 * File:   ComputeOneWF.h
 * Author: zigg
 *
 * Created on May 9, 2012, 10:51 PM
 */

#ifndef COMPUTEONEWF_H
#define	COMPUTEONEWF_H

class ComputeOneWF {
public:
    ComputeOneWF();
    ComputeOneWF(const ComputeOneWF& orig);
    virtual ~ComputeOneWF();
private:
    double E;
    double Esq;
    double accepted;
    
    int McSamples;
    bool importanceSampling;
    int thermalization;
    int nParticles;    
    double alpha, beta, w;
    bool usingJastrow;
    int dim;
    
    bool blocking;
};

#endif	/* COMPUTEONEWF_H */

