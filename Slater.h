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
    void computeGradient(int);
    double getRatio();

    rowvec getGradient() {
        return gradient;
    };
    void updateInverse();
    void acceptNewPosition();
    void init();
    double evaluate(const mat &);
    double getVariationalGradient();
private:
    Orbital *orbital;
    mat r_new;
    mat r_old;

    mat Dp;
    mat Dm;
    mat Dp_new;
    mat Dm_new;
    mat Dp_inv;
    mat Dm_inv;
    mat Dp_inv_new;
    mat Dm_inv_new;

    ivec nx;
    ivec ny;

    /*
    mat D; 
    mat DNew;
    mat DInv; 
    mat DInvNew; 
     */
    int dim;
    int N;
    int nParticles;
    int activeParticle;

    rowvec gradient;
};

#endif	/* SLATER_H */

