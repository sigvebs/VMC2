/* 
 * File:   WaveFunction.cpp
 * Author: Sigve
 * 
 * Created on May 9, 2012, 11:05 PM
 */

#include <mpi.h>
#include "WaveFunction.h"
#include "includes/lib.h"
#include "includes/ini.h"

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
    ini INIreader("QD.ini");
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    if (myRank == 0) {
        dt = INIreader.GetDouble("main", "dt");
    }
    MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    sqrtDt = sqrt(dt);
    D = 0.5;

    rOld = 2*randn(nParticles, dim);
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
    idum = orig.idum - time(NULL);
    stepLength = orig.stepLength;

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
    for (int i = 0; i < dim; i++) {
        rNew(activeParticle, i) = rOld(activeParticle, i) + stepLength * (ran3(&idum) - 0.5);
    }

    // Updating the Slater matrix
    slater->setPosition(rNew, activeParticle);
    slater->updateMatrix();

    //--------------------------------------------------------------------------
    // Metropolis acceptance test.
    //--------------------------------------------------------------------------
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
        rNew(activeParticle, i) = rOld(activeParticle, i) + D * qForceOld(activeParticle, i) * dt + 2 * D * sqrtDt * DRanNormalZigVec();

    // Updating Slater
    slater->setPosition(rNew, activeParticle);
    slater->updateMatrix();
    double R = WFRatio();
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
                    0.5 * (qForceOld(i, j) + qForce(i, j)) *
                    (D * dt * 0.5 * (qForceOld(i, j) - qForce(i, j)) - (rNew(i, j) - rOld(i, j)));
        }
    }

    greensFunction = exp(greensFunction);
    //--------------------------------------------------------------------------
    // Metropolis-Hastings acceptance test.
    //--------------------------------------------------------------------------
    R = R * R * greensFunction;
    test();
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
    rowvec gradientJastrow;

#define ANALYTICAL 1
#define ANALYTICAL_GRADIENT_DEBUG 0
#define NUMERICAL_GRADIENT 0
#define DEBUG_QFORCE 0

#if ANALYTICAL
    for (int i = 0; i < nParticles; i++) {
        // Slaters' gradient.
        slater->computeGradient(i);
        gradientSlater = slater->getGradient();

        // Jastrow's gradient.
        if (usingJastrow) {
            jastrow->computeGradient(rNew, i);
            gradientJastrow = jastrow->getGradient();
            q_f.row(i) = 2.0 * (gradientJastrow + gradientSlater);
        }
    }
#endif

#if ANALYTICAL_GRADIENT_DEBUG

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
    }
#endif

#if NUMERICAL_GRADIENT
#if DEBUG_QFORCE
    mat qfAnalytic = q_f;
#endif
    mat slaterNum = zeros(nParticles, dim);
    mat jastrowNum = zeros(nParticles, dim);
    q_f = zeros(nParticles, dim);
    double h = 0.00001;
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
    // Gradient Jastrow
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

    q_f = 2 * (slaterNum + jastrowNum);
    //--------------------------------------------------------------------------
#if DEBUG_QFORCE
    if (max(max(abs(qfAnalytic - q_f))) > 1e-4)
        cout
            << "--------------------------" << endl
            << "Max differance: " << max(max(abs(qfAnalytic - q_f))) << endl
            << "analytic = " << endl << qfAnalytic << endl
            << "numeric = " << endl << q_f << endl;
#endif
#endif

    return q_f;
}

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
    if (usingJastrow)
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
        rNew(activeParticle, i) = rOld(activeParticle, i) + D * qForceOld(activeParticle, i) * dt + 2 * D * sqrtDt * DRanNormalZigVec();

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
                    0.5 * (qForceOld(i, j) + qForce(i, j)) *
                    (D * dt * 0.5 * (qForceOld(i, j) - qForce(i, j)) - (rNew(i, j) - rOld(i, j)));
        }
    }

    greensFunction = exp(greensFunction);

    //--------------------------------------------------------------------------
    // Metropolis-Hastings acceptance test.
    //--------------------------------------------------------------------------
    double R = WFRatio();
    R = R * R * greensFunction;

    // Fixed node approximation
    if (R < 0) {
        return accepted;
    }

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

void WaveFunction::test() {
    mat dp = slater->Dp;
    mat dm = slater->Dm;
    mat dpNew = slater->DpNew;
    mat dmNew = slater->DmNew;
    mat dpInv = slater->DpInv;
    mat dmInv = slater->DmInv;
    mat dpInvNew = slater->DpInvNew;
    mat dmInvNew = slater->DmInvNew;
    /*
    cout 
            << "dp*dpInv = " << det(dp*dpInv) << endl
            << "dm*dmInv = " << det(dm*dmInv) << endl
            << "dpNew*dpInvNew = " << det(dpNew*dpInvNew) << endl
            << "dmNew*dmInvNew = " << det(dmNew*dmInvNew) << endl
            << endl;
     */

    // Old
    double errorThreshold = 1e-12;
    if (abs(det(dp * dpInv)) > 1 + errorThreshold || abs(det(dp * dpInv)) < 1 - errorThreshold) {
        cout << "dp*dpInv = " << dp * dpInv << endl;
        cout << "det(dp*dpInv) = " << det(dp * dpInv) << endl;
    }
    if (abs(det(dm * dmInv)) > 1 + errorThreshold || abs(det(dm * dmInv)) < 1 - errorThreshold) {
        cout << "dm*dmInv = " << dm * dmInv << endl;
        cout << "det(dm*dmInv) = " << det(dm * dmInv) << endl;
    }

    // New
    if (abs(det(dpNew * dpInvNew)) > 1 + errorThreshold || abs(det(dpNew * dpInvNew)) < 1 - errorThreshold) {
        cout << "dpNew*dpInvNew = " << dpNew * dpInvNew << endl;
        cout << "det(dpNew*dpInvNew) = " << det(dpNew * dpInvNew) << endl;
    }
    if (abs(det(dmNew * dmInvNew)) > 1 + errorThreshold || abs(det(dmNew * dmInvNew)) < 1 - errorThreshold) {
        cout << "dmNew*dmInvNew = " << dmNew * dmInvNew << endl;
        cout << "det(dmNew*dmInvNew) = " << det(dmNew * dmInvNew) << endl;
    }
}