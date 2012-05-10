/* 
 * File:   WaveFunction.cpp
 * Author: Sigve
 * 
 * Created on May 9, 2012, 11:05 PM
 */

#include "WaveFunction.h"
#include "includes/lib.h"

#include "includes/ziggurat.hpp"
//#include "includes/ziggurat.cpp"
#include "includes/zignor.h"
#include "includes/zignor.c"
#include "includes/zigrandom.h"
#include "includes/zigrandom.c"

////////////////////////////////////////////////////////////////////////////////

WaveFunction::WaveFunction() {
}

////////////////////////////////////////////////////////////////////////////////

WaveFunction::WaveFunction(int dim, int nParticles, long idum, Orbital *orbital, Jastrow *jastrow, Hamiltonian *hamiltonian)
: dim(dim), nParticles(nParticles), idum(idum), jastrow(jastrow), hamiltonian(hamiltonian) {

    int dum = (int) idum;
    RanNormalSetSeedZigVec(&dum, 200);

    slater = new Slater(dim, nParticles, orbital);

    E = 0;
    stepLength = 2.0;

    usingJastrow = true;
    if (jastrow == NULL)
        usingJastrow = false;

    // Initializing positions
    dt = 0.03;
    sqrtDt = sqrt(dt);
    D = 0.5;
    rOld = randn(nParticles, dim) * sqrtDt;
    rNew = rOld;
    slater->setPosition(rNew, 0);
    slater->init();
    qForce = newQForce();
    qForceOld = qForce;
}

////////////////////////////////////////////////////////////////////////////////

WaveFunction::WaveFunction(const WaveFunction& orig) {
}

////////////////////////////////////////////////////////////////////////////////

WaveFunction::~WaveFunction() {
    delete slater;
    delete hamiltonian;
}

////////////////////////////////////////////////////////////////////////////////
// Used for numerical derivatives.

double WaveFunction::evaluate(const mat &r) {
    double psi;

    psi = slater->evaluate(r);

    // Adding the Jastrow part
    if (usingJastrow)
        psi *= exp(jastrow->evaluate(r));

    return psi;
}

////////////////////////////////////////////////////////////////////////////////

bool WaveFunction::tryNewPositionBF(int activeParticle) {


    this->activeParticle = activeParticle;
    bool accepted = false;

    // Calculating new trial position.
    for (int i = 0; i < dim; i++)
        rNew(activeParticle, i) = rOld(activeParticle, i) + stepLength * (ran3(&idum) - 0.5);

    // Updating Slater
    slater->setPosition(rNew, activeParticle);
    slater->updateMatrix();

    //--------------------------------------------------------------------------
    // Metropolis acceptance test.
    double R = WFRatio();
    R = R * R;

    if (ran3(&idum) <= R) {
        accepted = true;
        rOld = rNew;
        slater->updateInverse();
        slater->acceptNewPosition();
        calculateEnergy();
    } else {
        // If the move is not accepted the position and quantum force is reset.
        rNew = rOld;
    }
    //--------------------------------------------------------------------------
    return accepted;
}

////////////////////////////////////////////////////////////////////////////////

bool WaveFunction::tryNewPosition(int active) {
#if 0
    cout
            << "\t dt = " << dt
            << "\t sqrtdt = " << sqrtDt
            << "\t D = " << D
            << endl;
#endif
    activeParticle = active;
    bool accepted = false;

    // Using the quantum force to calculate a new position.
    for (int i = 0; i < dim; i++)
        rNew(activeParticle, i) = rOld(activeParticle, i) + D * qForceOld(activeParticle, i) * dt + DRanNormalZigVec() * sqrtDt;

    // Updating Slater
    slater->setPosition(rNew, activeParticle);
    slater->updateMatrix();

    // Updating the quantum force.
    qForce = newQForce();

    //--------------------------------------------------------------------------
    // Calculating the ratio between the Green's functions.     
    double greens_function = 0;

    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < dim; j++) {
            greens_function +=
                    (qForceOld(i, j) + qForce(i, j)) *
                    (D * dt * 0.5 * (qForceOld(i, j) - qForce(i, j)) - rNew(i, j) + rOld(i, j));
        }
    }

    greens_function = exp(0.5 * greens_function);

    //--------------------------------------------------------------------------
    // Metropolis-Hastings acceptance test.
    double R = WFRatio();
    R = R * R *greens_function;

    if (ran3(&idum) <= R) {
        accepted = true;
        rOld = rNew;
        qForceOld = qForce;
        slater->updateInverse();
        slater->acceptNewPosition();
        calculateEnergy();
    } else {
        // If the move is not accepted the position and quantum force is reset.
        rNew = rOld;
        qForce = qForceOld;
    }
    //--------------------------------------------------------------------------

    return accepted;
}

////////////////////////////////////////////////////////////////////////////////

double WaveFunction::sampleEnergy() {
    return E;
}

////////////////////////////////////////////////////////////////////////////////

double WaveFunction::WFRatio() {
    double R = slater->getRatio();

    if (usingJastrow)
        R *= jastrow->getRatio(rNew, rOld);

    return R;
}

////////////////////////////////////////////////////////////////////////////////
// CHECK

mat WaveFunction::newQForce() {
    mat q_f = zeros(nParticles, dim);
    rowvec gradientSlater;
    rowvec gradientJastrow = zeros(1, dim);
    double R = slater->getRatio();

    for (int i = 0; i < nParticles; i++) {

        // Orbitals' gradient.
        slater->computeGradient(i);
        gradientSlater = slater->getGradient();

        // Jastrow's gradient.
        if (usingJastrow) {
            jastrow->computeGradient(rNew, i);
            gradientJastrow = jastrow->getGradient();
        }

        q_f.row(i) = 2.0 * (gradientJastrow + gradientSlater / R);
    }

    return q_f;
}

////////////////////////////////////////////////////////////////////////////////
// CHECK

void WaveFunction::calculateEnergy() {
    double EKin = 0;
    double EPot = 0;

    // Looping through all particles
    for (int i = 0; i < nParticles; i++) {
        // Finding the Orbitals' Laplacian.
        slater->computeGradient(i);
        EKin += slater->getLaplacian(i);

        if (usingJastrow) {
            // Finding the Jastrow factors Laplacian
            jastrow->computeGradient(rOld, i);
            EKin += jastrow->getLaplacian(rOld, i);

            // Dot product between the gradients of the orbital- and Jastrow function.
            rowvec gradient_orbital = slater->getGradient();
            rowvec gradient_jastrow = jastrow->getGradient();

            EKin += 2 * dot(gradient_orbital, gradient_jastrow);
        }
    }

    EKin = -0.5 * EKin;
    EPot = hamiltonian->getEnergy(rOld);

    // Total energy
    E = EKin + EPot;
}

////////////////////////////////////////////////////////////////////////////////

void WaveFunction::setOptimalStepLength() {
    double min = 0.01;
    double max = 3.5;
    double tolerance = 0.01;

    while ((max - min) > tolerance) {
        if (difference(min) * difference((min + max) / 2) < 0)
            max = (min + max) / 2;
        else
            min = (min + max) / 2;
    }

    stepLength = (min + max) / 2;
}


////////////////////////////////////////////////////////////////////////////////

double WaveFunction::difference(double step_length) {
    double SlSamples = 1e3;
    double accepted;
    double accepted_tmp;
    stepLength = step_length;

    for (int i = 0; i < 2 * SlSamples; i++) {
        for (int j = 0; j < nParticles; j++) {

            accepted_tmp = tryNewPositionBF(j);
            if (i > SlSamples)
                accepted += accepted_tmp;
        }
    }

    return (double) accepted / (SlSamples * nParticles) - 0.5;
}