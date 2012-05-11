/* 
 * File:   Slater.cpp
 * Author: Sigve
 * 
 * Created on May 10, 2012, 8:41 AM
 */

#include "Slater.h"

////////////////////////////////////////////////////////////////////////////////

Slater::Slater() {
}

////////////////////////////////////////////////////////////////////////////////

Slater::Slater(int dim, int nParticles, Orbital *orbital)
: dim(dim), nParticles(nParticles), N(nParticles / 2), orbital(orbital) {
    Dp = zeros(N, N);
    Dm = zeros(N, N);

    // Computing the Quantum Numbers.
    nx.set_size(N);
    ny.set_size(N);
    int l = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            nx(l) = i - j;
            ny(l) = j;
            l++;
            // Breaking the loop if we have enough numbers.
            if (l == N)
                i = j = N;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

Slater::Slater(const Slater& orig) {
}

////////////////////////////////////////////////////////////////////////////////

Slater::~Slater() {
}

////////////////////////////////////////////////////////////////////////////////

void Slater::init() {

    // Updating the whole matrix.
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Dp(i, j) = orbital->evaluate(rNew.row(j), nx(i), ny(i));
            Dm(i, j) = orbital->evaluate(rNew.row(j + N), nx(i), ny(i));
        }
    }
    rOld = rNew;
    DpNew = Dp;
    DmNew = Dm;

    // Calulating the inverse using Armadillo.
    DpInv = inv(Dp);
    DmInv = inv(Dm);

    DpInvNew = DpInv;
    DmInvNew = DmInv;
}

////////////////////////////////////////////////////////////////////////////////

double Slater::getRatio() {
    double R = 0;
    int i = activeParticle;

    if (i < N) { // Spin up
        for (int j = 0; j < N; j++) {
            R += DpNew(j, i) * DpInv(j, i);
        }
    } else {
        for (int j = 0; j < N; j++) {
            R += DmNew(j, i - N) * DmInv(j, i - N);
        }
    }
    return R;
}

////////////////////////////////////////////////////////////////////////////////

void Slater::updateInverse() {

    int i = activeParticle;
    DpInvNew = DpInv;
    DmInv = DmInvNew;

    double S;
    // TMP solution for the inverse.

    double R = getRatio();

    if (i < N) { // Spin up
        for (int j = 0; j < N; j++) {
            if (j != i) {
                S = 0;
                for (int l = 0; l < N; l++)
                    S += DpNew(l, i) * DpInv(l, j);
            }
            for (int l = 0; l < N; l++)
                DpInvNew(l, j) = DpInv(l, j) - (DpInv(l, i) / R) * S;
        }
        for (int l = 0; l < N; l++)
            DpInvNew(l, i) = DpInv(l, i) / R;
    } else { // Spin down
        for (int j = 0; j < N; j++) {
            if (j != i - N) {
                S = 0;
                for (int l = 0; l < N; l++)
                    S += DmNew(l, i - N) * DmInv(l, j);
            }
            for (int l = 0; l < N; l++)
                DmInvNew(l, j) = DmInv(l, j) - (DmInv(l, i - N) / R) * S;
        }
        for (int l = 0; l < N; l++)
            DmInvNew(l, i - N) = DmInv(l, i - N) / R;
    }

}

////////////////////////////////////////////////////////////////////////////////

void Slater::setPosition(const mat & r, int active_particle) {
    this->activeParticle = active_particle;
    rNew = r;
}

////////////////////////////////////////////////////////////////////////////////

void Slater::updateMatrix() {
    if (activeParticle < N) {
        for (int i = 0; i < N; i++) {
            DpNew(i, activeParticle) = orbital->evaluate(rNew.row(activeParticle), nx(i), ny(i));
        }
    } else {
        for (int i = 0; i < N; i++) {
            DmNew(i, activeParticle - N) = orbital->evaluate(rNew.row(activeParticle), nx(i), ny(i));
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void Slater::acceptNewPosition() {

    if (activeParticle < N) {
        for (int i = 0; i < N; i++)
            Dp(i, activeParticle) = DpNew(i, activeParticle);
        DpInv = DpInvNew;
    } else {
        for (int i = 0; i < N; i++)
            Dm(i, activeParticle - N) = DmNew(i, activeParticle - N);
        DmInv = DmInvNew;
    }

    rOld = rNew;
}

////////////////////////////////////////////////////////////////////////////////

double Slater::getLaplacian(int i) {

    double sum = 0;
    if (i < N) {
        for (int j = 0; j < N; j++) { // Spin up.
            sum += orbital->evaluateLaplacian(rNew.row(i), nx(j), ny(j)) * DpInv(j, i);
        }
    } else {
        for (int j = 0; j < N; j++) { // Spin down.
            sum += orbital->evaluateLaplacian(rNew.row(i), nx(j), ny(j)) * DmInv(j, i - N);
        }
    }

    return sum;
}

////////////////////////////////////////////////////////////////////////////////

void Slater::computeGradient(int i) {
    gradient = zeros(1, dim);

    if (i < N) {
        for (int j = 0; j < N; j++) { // Spin up.
            gradient += orbital->getGradient(rNew.row(i), nx(j), ny(j)) * DpInv(j, i);
        }
    } else {
        for (int j = 0; j < N; j++) { // Spin down.
            gradient += orbital->getGradient(rNew.row(i), nx(j), ny(j)) * DmInv(j, i - N);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

double Slater::evaluate(const mat &r) {
    mat D_p = zeros(N, N);
    mat D_m = zeros(N, N);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            D_p(j, i) = orbital->evaluate(r.row(i), nx(j), ny(j));
            D_m(j, i) = orbital->evaluate(r.row(i + N), nx(j), ny(j));
        }
    }

    return det(D_p) * det(D_m);
}

////////////////////////////////////////////////////////////////////////////////

double Slater::getVariationalGradient() {
    double gradient = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
//            gradient += DpInv(j, i) * orbital->variationalDerivative(rNew.row(i), nx(j), ny(j));
//            gradient += DmInv(j, i) * orbital->variationalDerivative(rNew.row(i + N), nx(j), ny(j));
            gradient += DpInv(i,j) * orbital->variationalDerivative(rNew.row(j), nx(i), ny(i));
            gradient += DmInv(i,j) * orbital->variationalDerivative(rNew.row(j + N), nx(i), ny(i));
        }
    }

    return gradient;
}
