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
            Dp(j, i) = orbital->evaluate(r_new.row(i), nx(j), ny(j));
            Dm(j, i) = orbital->evaluate(r_new.row(i + N), nx(j), ny(j));
        }
    }
    Dp_new = Dp;
    Dm_new = Dm;

    // Calulating the inverse using Armadillo.
    Dp_inv = inv(Dp).t();
    Dm_inv = inv(Dm).t();

    Dp_inv_new = Dp_inv;
    Dm_inv_new = Dm_inv;
}

////////////////////////////////////////////////////////////////////////////////
double Slater::getRatio() {
    double R = 0;
    int i = activeParticle;

    if (i < N) { // Spin up
        for (int j = 0; j < N; j++) {
            R += Dp_new(j, i) * Dp_inv(j, i);
        }
    } else {
        for (int j = 0; j < N; j++) {
            R += Dm_new(j, i - N) * Dm_inv(j, i - N);
        }
    }
    return R;
}

////////////////////////////////////////////////////////////////////////////////
void Slater::updateInverse() {

    int i = activeParticle;
    Dp_inv_new = Dp_inv;
    Dm_inv = Dm_inv_new;

    double S;
    // TMP solution for the inverse.

    double R = getRatio();

    if (i < N) { // Spin up
        for (int j = 0; j < N; j++) {
            if (j != i) {
                S = 0;
                for (int l = 0; l < N; l++)
                    S += Dp_new(l, i) * Dp_inv(l, j);
            }
            for (int l = 0; l < N; l++)
                Dp_inv_new(l, j) = Dp_inv(l, j) - (Dp_inv(l, i) / R) * S;
        }
        for (int l = 0; l < N; l++)
            Dp_inv_new(l, i) = Dp_inv(l, i) / R;
    } else { // Spin down
        for (int j = 0; j < N; j++) {
            if (j != i - N) {
                S = 0;
                for (int l = 0; l < N; l++)
                    S += Dm_new(l, i - N) * Dm_inv(l, j);
            }
            for (int l = 0; l < N; l++)
                Dm_inv_new(l, j) = Dm_inv(l, j) - (Dm_inv(l, i - N) / R) * S;
        }
        for (int l = 0; l < N; l++)
            Dm_inv_new(l, i - N) = Dm_inv(l, i - N) / R;
    }

}

////////////////////////////////////////////////////////////////////////////////
void Slater::setPosition(const mat & r, int active_particle) {
    this->activeParticle = active_particle;
    r_new = r;
}

////////////////////////////////////////////////////////////////////////////////
void Slater::updateMatrix() {
    if (activeParticle < N) {
        for (int i = 0; i < N; i++) {
            Dp_new(i, activeParticle) = orbital->evaluate(r_new.row(activeParticle), nx(i), ny(i));
        }
    } else {
        for (int i = 0; i < N; i++) {
            Dm_new(i, activeParticle - N) = orbital->evaluate(r_new.row(activeParticle), nx(i), ny(i));
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
void Slater::acceptNewPosition() {

    if (activeParticle < N) {
        for (int i = 0; i < N; i++)
            Dp(i, activeParticle) = Dp_new(i, activeParticle);
        Dp_inv = Dp_inv_new;
    } else {
        for (int i = 0; i < N; i++)
            Dm(i, activeParticle - N) = Dm_new(i, activeParticle - N);
        Dm_inv = Dm_inv_new;
    }

    // Should be updated more efficiently.
    r_old = r_new;
}

////////////////////////////////////////////////////////////////////////////////
double Slater::getLaplacian(int i) {

    double sum = 0;
    if (i < N) {
        for (int j = 0; j < N; j++) { // Spin up.
            sum += orbital->evaluateLaplacian(r_new.row(i), nx(j), ny(j)) * Dp_inv(j, i);
        }
    } else {
        for (int j = 0; j < N; j++) { // Spin down.
            sum += orbital->evaluateLaplacian(r_new.row(i), nx(j), ny(j)) * Dm_inv(j, i - N);
        }
    }

    return sum;
}

////////////////////////////////////////////////////////////////////////////////
void Slater::computeGradient(int i) {
    gradient = zeros(1, dim);

    if (i < N) {
        for (int j = 0; j < N; j++) { // Spin up.
            gradient += orbital->getGradient(r_new.row(i), nx(j), ny(j)) * Dp_inv(j, i);
        }
    } else {
        for (int j = 0; j < N; j++) { // Spin down.
            gradient += orbital->getGradient(r_new.row(i), nx(j), ny(j)) * Dm_inv(j, i - N);
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
            gradient += Dp_inv(j, i) * orbital->variationalDerivative(r_new.row(i), nx(j), ny(j));
            gradient += Dm_inv(j, i) * orbital->variationalDerivative(r_new.row(i + N), nx(j), ny(j));
        }
    }
    
    return gradient;
}
