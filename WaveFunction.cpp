/* 
 * File:   WaveFunction.cpp
 * Author: Sigve
 * 
 * Created on May 9, 2012, 11:05 PM
 */

#include "WaveFunction.h"
#include "includes/lib.h"

#include "includes/ziggurat.hpp"
#include "includes/zignor.h"
#include "includes/zignor.c"
#include "includes/zigrandom.h"
#include "includes/zigrandom.c"

////////////////////////////////////////////////////////////////////////////////

WaveFunction::WaveFunction() {
}

////////////////////////////////////////////////////////////////////////////////

WaveFunction::WaveFunction(int dim, int nParticles, long idum, Orbital *orbital, Jastrow *jastrow, Hamiltonian *hamiltonian)
: dim(dim), nParticles(nParticles), idum(idum), orbital(orbital), jastrow(jastrow), hamiltonian(hamiltonian) {
    int dum = (int) idum;
    RanNormalSetSeedZigVec(&dum, 200);

    slater = new Slater(dim, nParticles, orbital);

    E = 0;
    stepLength = 2.0;

    usingJastrow = true;
    if (jastrow == NULL)
        usingJastrow = false;

    // Initializing positions
    dt = 0.0005;
    sqrtDt = sqrt(dt);
    D = 0.5;

    rOld = randu(nParticles, dim);
    rNew = rOld;
    slater->setPosition(rNew, 0);
    slater->init();
    slater->acceptPosition();
    qForce = newQForce();
    qForceOld = qForce;
}

////////////////////////////////////////////////////////////////////////////////

WaveFunction::WaveFunction(const WaveFunction& orig) {
    usingJastrow = orig.usingJastrow;
    activeParticle = orig.activeParticle;
    idum = orig.idum;

    rOld = orig.rOld;
    rNew = orig.rNew;

    qForce = orig.qForce;
    qForceOld = orig.qForceOld;

    jastrow = orig.jastrow;
    slater = orig.slater->clone();
    orbital = orig.orbital;
    hamiltonian = orig.hamiltonian;

    E = orig.E;
    dim = orig.dim;
    nParticles = orig.nParticles;

    dt = orig.dt;
    sqrtDt = orig.sqrtDt;
    D = orig.D;

    int dum = (int) idum;
    RanNormalSetSeedZigVec(&dum, 200);
}

////////////////////////////////////////////////////////////////////////////////

WaveFunction::~WaveFunction() {
    delete slater;
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
        slater->acceptPosition();
        calculateEnergy();
    } else {
        // If the move is not accepted the position is reset.
        slater->rejectPosition();
        for (int d = 0; d < dim; d++)
            rNew(activeParticle, d) = rOld(activeParticle, d);
    }
    //--------------------------------------------------------------------------
    return accepted;
}

////////////////////////////////////////////////////////////////////////////////

bool WaveFunction::tryNewPosition(int active) {
    this->activeParticle = active;
    bool accepted = false;

    // Calculating new trial position.
    for (int i = 0; i < dim; i++)
        rNew(activeParticle, i) = rOld(activeParticle, i) + D * qForceOld(activeParticle, i) * dt + DRanNormalZigVec() * sqrtDt;

    // Updating Slater
    slater->setPosition(rNew, activeParticle);
    slater->updateMatrix();
    slater->updateInverse();

    // Updating the quantum force.
    qForce = newQForce();

    //--------------------------------------------------------------------------
    // Calculating the ratio between the Green's functions.     
    double greensFunction = 0;

    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < dim; j++) {
            greensFunction +=
                    (qForceOld(i, j) + qForce(i, j)) *
                    (D * dt * 0.5 * (qForceOld(i, j) - qForce(i, j)) - rNew(i, j) + rOld(i, j));
        }
    }

    greensFunction = exp(0.5 * greensFunction);

    //--------------------------------------------------------------------------
    // Metropolis-Hastings acceptance test.
    double R = WFRatio();
    R = R * R * greensFunction;

    if (ran3(&idum) <= R) {
        accepted = true;
        rOld = rNew;
        qForceOld = qForce;
        slater->acceptPosition();
        calculateEnergy();
    } else {
        // If the move is not accepted the position is reset.
        qForce = qForceOld;
        slater->rejectPosition();
        rNew = rOld;
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

mat WaveFunction::newQForce() {
    mat q_f = zeros(nParticles, dim);
    rowvec gradientSlater;
    rowvec gradientJastrow = zeros(1, dim);
    //double R = slater->getRatio();

#define ANALYTICAL_GRADIENT 1
#define NUMERICAL_GRADIENT 0

#if ANALYTICAL_GRADIENT

    //double errorThreshold = 1e-2;
    //rowvec gradientSlaterNum;
    for (int i = 0; i < nParticles; i++) {

        // Orbitals' gradient.
        slater->computeGradient(i);
        gradientSlater = slater->getGradient();
        //gradientSlaterNum = slater->getGradientNumerical(i);

        // Jastrow's gradient.
        if (usingJastrow) {
            jastrow->computeGradient(rNew, i);
            gradientJastrow = jastrow->getGradient();
        }

        /*
        if (max(max(gradientSlater - gradientSlaterNum)) > errorThreshold)
            cout
                << " Active particle = " << activeParticle
                << " i = " << i << endl
                << gradientSlater - gradientSlaterNum << endl;
         */
        q_f.row(i) = 2.0 * (gradientJastrow + gradientSlater);
        //q_f.row(i) = 2.0 * (gradientJastrow + gradientSlaterNum);
    }

#endif

#if NUMERICAL_GRADIENT
#if DEBUG_QFORCE
    mat slaterNum = zeros(nParticles, dim);
    mat jastrowNum = zeros(nParticles, dim);
#endif
    mat slaterNum = zeros(nParticles, dim);
    mat jastrowNum = zeros(nParticles, dim);
    q_f = zeros(nParticles, dim);
    double h = 0.001;
    double wfminus, wfplus, wfold;
    mat rPlus = zeros(nParticles, dim);
    mat rMinus = zeros(nParticles, dim);
    rPlus = rMinus = rNew;
    //--------------------------------------------------------------------------
    wfold = slater->evaluate(rNew);
    // Gradient slater
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < dim; j++) {
            rPlus(i, j) = rNew(i, j) + h;
            rMinus(i, j) = rNew(i, j) - h;
            wfminus = slater->evaluate(rMinus);
            wfplus = slater->evaluate(rPlus);
            slaterNum(i, j) = (wfplus - wfminus) / (2 * h * wfold);
            rPlus(i, j) = rNew(i, j);
            rMinus(i, j) = rNew(i, j);
        }
    }
    //--------------------------------------------------------------------------
    // Gradient slater
    wfold = exp(jastrow->evaluate(rNew));
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < dim; j++) {
            rPlus(i, j) = rNew(i, j) + h;
            rMinus(i, j) = rNew(i, j) - h;
            wfminus = exp(jastrow->evaluate(rMinus));
            wfplus = exp(jastrow->evaluate(rPlus));
            jastrowNum(i, j) = (wfplus - wfminus) / (2 * h * wfold);
            rPlus(i, j) = rNew(i, j);
            rMinus(i, j) = rNew(i, j);
        }
    }
    
    q_f = 2*(slaterNum + jastrowNum);
    //--------------------------------------------------------------------------
#endif

    return q_f;
}
/*
     q_f = zeros(nParticles, dim);
    double h = 0.001;

    double wfminus, wfplus, wfold;

    wfold = evaluate(rNew);

    mat rPlus = zeros(nParticles, dim);
    mat rMinus = zeros(nParticles, dim);

    rPlus = rMinus = rNew;

    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < dim; j++) {
            rPlus(i, j) = rNew(i, j) + h;
            rMinus(i, j) = rNew(i, j) - h;
            wfminus = evaluate(rMinus);
            wfplus = evaluate(rPlus);

            q_f(i, j) = 2 * (wfplus - wfminus) / (2 * h * wfold);
            rPlus(i, j) = rNew(i, j);
            rMinus(i, j) = rNew(i, j);
        }
    }
 */
////////////////////////////////////////////////////////////////////////////////

void WaveFunction::calculateEnergy() {
    double EKin = 0;
    double EPot = 0;
    rowvec gradientOrbital;
    rowvec gradientJastrow;

    double EKinAnalytic;
    double EKinNumeric;

#define ANALYTIC_ENERGY 1
#define NUMERICAL_ENERGY 0
#define PRINT_DIFF_ENERGY 0


#if ANALYTIC_ENERGY
    EKin = 0;

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
            gradientOrbital = slater->getGradient();
            gradientJastrow = jastrow->getGradient();

            EKin += 2 * dot(gradientOrbital, gradientJastrow);
        }
    }

    EKinAnalytic = EKin;

#endif
#if NUMERICAL_ENERGY
    EKin = 0;
    double h = 0.0001;

    double wfminus, wfplus, wfold;

    wfold = evaluate(rNew);

    mat rPlus = zeros(nParticles, dim);
    mat rMinus = zeros(nParticles, dim);

    rPlus = rMinus = rNew;

    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < dim; j++) {
            rPlus(i, j) = rNew(i, j) + h;
            rMinus(i, j) = rNew(i, j) - h;

            wfminus = evaluate(rMinus);
            wfplus = evaluate(rPlus);

            EKin += (wfplus - 2 * wfold + wfminus) / (h * h * wfold);

            rPlus(i, j) = rNew(i, j);
            rMinus(i, j) = rNew(i, j);
        }
    }
#if    PRINT_DIFF_ENERGY
    EKinNumeric = EKin;
    double error = EKinNumeric - EKinAnalytic;
    if (abs(error) > 1e-2) {
        EPot = hamiltonian->getEnergy(rOld);
        cout << "-------------\n";
        cout << error << endl;
        cout << "E anal = " << -0.5 * EKinAnalytic + EPot << endl;
        cout << "E num = " << -0.5 * EKinNumeric + EPot << endl;
    }
#endif

#endif

    EKin = -0.5 * EKin;
    EPot = hamiltonian->getEnergy(rOld);

    // Total energy
    E = EKin + EPot;
}

////////////////////////////////////////////////////////////////////////////////

rowvec WaveFunction::getVariationGradient() {
    rowvec varGradient = zeros(1, 2);

    // Gradient_slater;
    varGradient(0) = slater->getVariationalGradient();

    // Gradient Jastrow.
    if (usingJastrow)
        varGradient(1) = jastrow->getVariationalGradient(rOld);

    return varGradient;
}

////////////////////////////////////////////////////////////////////////////////

void WaveFunction::setNewVariationalParameters(double alpha, double beta) {
    orbital->setNewAlpha(alpha);
    jastrow->setNewBeta(beta);
    slater->setPosition(rOld, 0);
    slater->init();

    // One loop over all particles to adjust the system.
    for (int i = 0; i < nParticles; i++)
        tryNewPosition(i);
}
////////////////////////////////////////////////////////////////////////////////

void WaveFunction::setOptimalStepLength() {
    double min = 0.01;
    double max = 2.8;
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
    double SlSamples = 1e4;
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

////////////////////////////////////////////////////////////////////////////////

bool WaveFunction::DMCtryNewPosition(int active) {
    this->activeParticle = active;
    bool accepted = false;

    // Calculating new trial position.
    for (int i = 0; i < dim; i++)
        rNew(activeParticle, i) = rOld(activeParticle, i) + D * qForceOld(activeParticle, i) * dt + DRanNormalZigVec() * sqrtDt;

    // Updating Slater
    slater->setPosition(rNew, activeParticle);
    slater->updateMatrix();
    slater->updateInverse();

    // Updating the quantum force.
    qForce = newQForce();

    //--------------------------------------------------------------------------
    // Calculating the diffusion ratio between the Green's functions.     
    //--------------------------------------------------------------------------
    double greensFunction = 0;

    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < dim; j++) {
            greensFunction +=
                    (qForceOld(i, j) + qForce(i, j)) *
                    (D * dt * 0.5 * (qForceOld(i, j) - qForce(i, j)) - rNew(i, j) + rOld(i, j));
        }
    }

    greensFunction = exp(0.5 * greensFunction);
    //--------------------------------------------------------------------------
    // Metropolis-Hastings acceptance test.
    //--------------------------------------------------------------------------
    double R = WFRatio();
    R = R * R * greensFunction;

    // Fixed node approximation
    if (R < 0)
        return accepted;

    if (ran3(&idum) <= R) {
        accepted = true;
        rOld = rNew;
        qForceOld = qForce;
        slater->acceptPosition();
        calculateEnergy();
    } else {
        // If the move is not accepted the position is reset.
        qForce = qForceOld;
        slater->rejectPosition();
        rNew = rOld;
    }
    //--------------------------------------------------------------------------
    return accepted;
}

////////////////////////////////////////////////////////////////////////////////