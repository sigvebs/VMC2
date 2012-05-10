/* 
 * File:   Hamiltonian.h
 * Author: zigg
 *
 * Created on May 10, 2012, 9:46 AM
 */

#ifndef HAMILTONIAN_H
#define	HAMILTONIAN_H

#include <armadillo>
using namespace arma;

class Hamiltonian {
public:
    Hamiltonian();

    Hamiltonian(int dim, int nParticles, bool usingJastrow) : dim(dim), nParticles(nParticles), usingJastrow(usingJastrow) {
    };
    Hamiltonian(const Hamiltonian& orig);
    virtual ~Hamiltonian();

    double getEnergy(const mat &);
protected:
    virtual double interaction(const mat &) = 0;
    virtual double potential(const mat &) = 0;
    int dim;
    int nParticles;
    bool usingJastrow;
};

#endif	/* HAMILTONIAN_H */

