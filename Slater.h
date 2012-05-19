/* 
 * File:   Slater.h
 * Author: zigg
 *
 * Created on May 10, 2012, 8:41 AM
 */

#ifndef SLATER_H
#define	SLATER_H

#include <armadillo>
using namespace arma;

#include "Orbital.h"

class Slater {
public:
    Slater();
    Slater(int, int, Orbital *);
    Slater(const Slater& orig);
    virtual ~Slater();

    void updateMatrix();
    void setPosition(const mat &, int);
    double getLaplacian(int);
    double getLaplacianNumerical(int);
    void computeGradient(int);
    double getRatio();

    rowvec getGradient() {
        return gradient;
    };

    rowvec getGradientNumerical(int);

    void updateInverse();
    void updateInverseNumeric();
    void acceptPosition();
    void rejectPosition();
    void init();
    double evaluate(const mat &);
    double getVariationalGradient();

    virtual Slater* clone() const {
        return new Slater(*this);
    }

public: // TODO: Change to private
    Orbital *orbital;
    mat rNew;
    mat rOld;

    mat Dp;
    mat Dm;
    mat DpNew;
    mat DmNew;
    mat DpInv;
    mat DmInv;
    mat DpInvNew;
    mat DmInvNew;

    ivec nx;
    ivec ny;

    int dim;
    int N;
    int nParticles;
    int activeParticle;

    rowvec gradient;
};

#endif	/* SLATER_H */

