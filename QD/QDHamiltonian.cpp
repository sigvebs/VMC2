/* 
 * File:   QDHamiltonian.cpp
 * Author: zigg
 * 
 * Created on May 10, 2012, 9:51 AM
 */

#include "QDHamiltonian.h"

////////////////////////////////////////////////////////////////////////////////

QDHamiltonian::QDHamiltonian() {
}

////////////////////////////////////////////////////////////////////////////////

QDHamiltonian::QDHamiltonian(int dim, int nParticles, double w, bool usingJastrow) : w(w), Hamiltonian(dim, nParticles, usingJastrow) {

#if 0
    cout
            << "\tdim = " << this->dim
            << "\tnParticles = " << this->nParticles
            << "\tw = " << this->w
            << "\tusingJastrow = " << this->usingJastrow
            << endl;

#endif
}
////////////////////////////////////////////////////////////////////////////////

QDHamiltonian::QDHamiltonian(const QDHamiltonian& orig) {
}

////////////////////////////////////////////////////////////////////////////////

QDHamiltonian::~QDHamiltonian() {
}

////////////////////////////////////////////////////////////////////////////////

double QDHamiltonian::interaction(const mat &r) {
    double r_ij;
    double ePotential = 0;

    // Contribution from electron-electron potential  
    for (int i = 0; i < nParticles; i++) {
        for (int j = i + 1; j < nParticles; j++) {
            r_ij = norm(r.row(i) - r.row(j), 2);
            ePotential += 1.0 / r_ij;
        }
    }
    return ePotential;
}

////////////////////////////////////////////////////////////////////////////////
// Harmonic Oscillator potential.

double QDHamiltonian::potential(const mat &r) {
    double potential = 0;

    for (int i = 0; i < nParticles; i++) {
        potential += dot(r.row(i), r.row(i));
    }

    return 0.5 * w * w*potential;
}