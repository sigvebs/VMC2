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
    orbital = orig.orbital;
    rNew = orig.rNew;
    rOld = orig.rOld;

    Dp = orig.Dp;
    Dm = orig.Dm;
    DpNew = orig.DpNew;
    DmNew = orig.DmNew;
    DpInv = orig.DpInv;
    DmInv = orig.DmInv;
    DpInvNew = orig.DpInvNew;
    DmInvNew = orig.DmInvNew;

    nx = orig.nx;
    ny = orig.ny;

    dim = orig.dim;
    N = orig.N;
    nParticles = orig.nParticles;
    activeParticle = orig.activeParticle;

    gradient = orig.gradient;

}

////////////////////////////////////////////////////////////////////////////////

Slater::~Slater() {
    //delete orbital;
}

////////////////////////////////////////////////////////////////////////////////

void Slater::init() {

    // Updating the whole matrix.
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Dp(i, j) = orbital->evaluate(rNew.row(i), nx(j), ny(j));
            Dm(i, j) = orbital->evaluate(rNew.row(i + N), nx(j), ny(j));
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
            R += DpNew(i, j) * DpInv(j, i);
        }
    } else {
        for (int j = 0; j < N; j++) {
            R += DmNew(i - N, j) * DmInv(j, i - N);
        }
    }

    return R;
}

////////////////////////////////////////////////////////////////////////////////

void Slater::updateInverse() {

    double R = getRatio();

    //--------------------------------------------------------------------------
    // Spin up
    //--------------------------------------------------------------------------
    if (activeParticle < N) { 
        int i = activeParticle;
        rowvec S = zeros(1, N);

        for (int j = 0; j < N; j++)
            for (int l = 0; l < N; l++)
                S(j) += DpNew(i, l) * DpInv(l, j);

        for (int k = 0; k < N; k++) {
            for (int j = 0; j < N; j++) {
                if (j != i) {
                    DpInvNew(k, j) = DpInv(k, j) - DpInv(k, i) * S(j) / R;
                } else {
                    DpInvNew(k, j) = DpInv(k, i) / R;
                }
            }
        }
    }
    //--------------------------------------------------------------------------
    // Spin down
    //--------------------------------------------------------------------------
    if (activeParticle >= N) { 
        int i = activeParticle - N;
        rowvec S = zeros(1, N);

        for (int j = 0; j < N; j++)
            for (int l = 0; l < N; l++)
                S(j) += DmNew(i, l) * DmInv(l, j);

        for (int k = 0; k < N; k++) {
            for (int j = 0; j < N; j++) {
                if (j != i) {
                    DmInvNew(k, j) = DmInv(k, j) - DmInv(k, i) * S(j) / R;
                } else {
                    DmInvNew(k, j) = DmInv(k, i) / R;
                }
            }
        }
    }
    //----------------------------------------------------------------------

#define DEBUG_INVERSE_SLATER 0

#if DEBUG_INVERSE_SLATER
    double errorThreshold = 1e-6;

    if (abs(det(DpInvNew * DpNew) - 1) > errorThreshold)
        cout
            << "Inverse Failed: DpInv * Dp = " << endl
            << "Analytical = " << endl
            << DpInvNew * DpNew << endl;

    if (abs(det(DmInvNew * DmNew) - 1) > errorThreshold)
        cout
            << "Inverse Failed: DmInv * Dm = " << endl
            << "Analytical = " << endl
            << DmInvNew * DmNew << endl;
#endif
}

////////////////////////////////////////////////////////////////////////////////

void Slater::updateInverseNumeric() {
    DpInvNew = inv(DpNew);
    DmInvNew = inv(DmNew);
}

////////////////////////////////////////////////////////////////////////////////

void Slater::setPosition(const mat & r, int active_particle) {
    this->activeParticle = active_particle;
    rNew = r;
}

////////////////////////////////////////////////////////////////////////////////

void Slater::updateMatrix() {

    int i = activeParticle;
    if (activeParticle < N) {
        for (int j = 0; j < N; j++) {
            DpNew(i, j) = orbital->evaluate(rNew.row(i), nx(j), ny(j));
        }
    } else {
        for (int j = 0; j < N; j++) {
            DmNew(i - N, j) = orbital->evaluate(rNew.row(i), nx(j), ny(j));
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void Slater::acceptPosition() {

    int i = activeParticle;

    if (activeParticle < N) {
        for (int j = 0; j < N; j++) {
            Dp(i, j) = DpNew(i, j);
        }
        DpInv = DpInvNew;
    } else {
        for (int j = 0; j < N; j++) {
            Dm(i - N, j) = DmNew(i - N, j);
        }
        DmInv = DmInvNew;
    }

    for (int d = 0; d < dim; d++)
        rOld(activeParticle, d) = rNew(activeParticle, d);
}

////////////////////////////////////////////////////////////////////////////////

void Slater::rejectPosition() {

    int i = activeParticle;

    if (activeParticle < N) {
        for (int j = 0; j < N; j++) {
            DpNew(i, j) = Dp(i, j);
        }
        DpInvNew = DpInv;
    } else {
        for (int j = 0; j < N; j++) {
            DmNew(i - N, j) = Dm(i - N, j);
        }
        DmInvNew = DmInv;
    }

    for (int d = 0; d < dim; d++)
        rNew(activeParticle, d) = rOld(activeParticle, d);
}

////////////////////////////////////////////////////////////////////////////////

double Slater::getLaplacian(int i) {

    double sum = 0;
    if (i < N) {
        for (int j = 0; j < N; j++) { // Spin up.
            sum += orbital->evaluateLaplacian(rNew.row(i), nx(j), ny(j)) * DpInvNew(j, i);
        }
    } else {
        for (int j = 0; j < N; j++) { // Spin down.
            sum += orbital->evaluateLaplacian(rNew.row(i), nx(j), ny(j)) * DmInvNew(j, i - N);
        }
    }

    return sum;
}

////////////////////////////////////////////////////////////////////////////////

double Slater::getLaplacianNumerical(int i) {
    double h = 0.001;
    double laplaceNum = 0;

    double wfminus, wfplus, wfold;
    wfold = evaluate(rNew);

    mat rPlus = zeros(nParticles, dim);
    mat rMinus = zeros(nParticles, dim);

    rPlus = rMinus = rNew;

    for (int j = 0; j < dim; j++) {
        rPlus(i, j) = rNew(i, j) + h;
        rMinus(i, j) = rNew(i, j) - h;

        wfminus = evaluate(rMinus);
        wfplus = evaluate(rPlus);

        laplaceNum += (wfplus - 2 * wfold + wfminus) / (h * h * wfold);

        rPlus(i, j) = rNew(i, j);
        rMinus(i, j) = rNew(i, j);
    }

    return laplaceNum;
}
////////////////////////////////////////////////////////////////////////////////

void Slater::computeGradient(int i) {
    gradient = zeros(1, dim);

    if (i < N) {
        for (int j = 0; j < N; j++) { // Spin up.
            gradient += orbital->getGradient(rNew.row(i), nx(j), ny(j)) * DpInvNew(j, i);
        }
    } else {
        for (int j = 0; j < N; j++) { // Spin down.
            gradient += orbital->getGradient(rNew.row(i), nx(j), ny(j)) * DmInvNew(j, i - N);
        }
    }

}

////////////////////////////////////////////////////////////////////////////////

rowvec Slater::getGradientNumerical(int i) {
    double h = 0.001;
    mat gradientSlaterNumerical = zeros(1, dim);
    double wfminus, wfplus;
    double wfold = evaluate(rNew);

    mat rPlus = zeros(nParticles, dim);
    mat rMinus = zeros(nParticles, dim);
    rPlus = rMinus = rNew;

    for (int j = 0; j < dim; j++) {
        rPlus(i, j) = rNew(i, j) + h;
        rMinus(i, j) = rNew(i, j) - h;
        wfminus = evaluate(rMinus);
        wfplus = evaluate(rPlus);

        gradientSlaterNumerical(j) = (wfplus - wfminus) / (2 * h * wfold);
        rPlus(i, j) = rNew(i, j);
        rMinus(i, j) = rNew(i, j);
    }

    return gradientSlaterNumerical;
}

////////////////////////////////////////////////////////////////////////////////

double Slater::evaluate(const mat & r) {
    
    mat D_p = zeros(N, N);
    mat D_m = zeros(N, N);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            D_p(i, j) = orbital->evaluate(r.row(i), nx(j), ny(j));
            D_m(i, j) = orbital->evaluate(r.row(i + N), nx(j), ny(j));
        }
    }

    return det(D_p) * det(D_m);
}

////////////////////////////////////////////////////////////////////////////////

double Slater::getVariationalGradient() {
    double gradient = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            gradient += DpInv(i, j) * orbital->variationalDerivative(rNew.row(j), nx(i), ny(i));
            gradient += DmInv(i, j) * orbital->variationalDerivative(rNew.row(j + N), nx(i), ny(i));
        }
    }

    return gradient;
}
