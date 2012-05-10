/* 
 * File:   QDHamiltonian.h
 * Author: Sigve
 *
 * Created on May 10, 2012, 9:51 AM
 */

#ifndef QDHAMILTONIAN_H
#define	QDHAMILTONIAN_H

#include "../Hamiltonian.h"

class QDHamiltonian : public Hamiltonian {
public:
    QDHamiltonian();
    QDHamiltonian(int, int, double, bool);
    QDHamiltonian(const QDHamiltonian& orig);
    virtual ~QDHamiltonian();

    virtual double interaction(const mat &);
    virtual double potential(const mat &);
private:
    double w;
};

#endif	/* QDHAMILTONIAN_H */

