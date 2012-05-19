/* 
 * File:   newSlaterTest.cpp
 * Author: sigve
 *
 * Created on 14.mai.2012, 11:41:31
 */

#include <stdlib.h>
#include <iostream>
#include <armadillo>

#include "Slater.h"
#include "QD/QDOrbital.h"
#include "QD/QDJastrow.h"
#include "includes/lib.h"

using namespace arma;

/*
 * Simple C++ Test Suite
 */

// Test of the updating algorithm of the inverse Slater determinant.

void testingInverse() {
    std::cout << "newSlaterTest test 1" << std::endl;
    double errorThreshold = 1e-4;

    int dim = 2;
    int nParticles = 12;
    double alpha = 0.987;
    double w = 1.0;
    long idum = -1;
    int cycles = 20;

    mat rInit = randu(nParticles, dim);
    mat rNew = rInit;
    mat rOld = rInit;
    //--------------------------------------------------------------------------
    // Analytic
    Orbital *orbital = new QDOrbital(dim, alpha, w);
    Slater *slater = new Slater(dim, nParticles, orbital);

    slater->setPosition(rInit, 0);
    slater->init();
    slater->acceptPosition();

    //--------------------------------------------------------------------------
    // Numeric
    Orbital *orbitalNumeric = new QDOrbital(dim, alpha, w);
    Slater *slaterNumeric = new Slater(dim, nParticles, orbitalNumeric);

    slaterNumeric->setPosition(rInit, 0);
    slaterNumeric->init();
    slaterNumeric->acceptPosition();

    //--------------------------------------------------------------------------
    // Moving one particle

    for (int n = 0; n < cycles; n++) {
        for (int i = 0; i < nParticles; i++) {

            for (int d = 0; d < dim; d++)
                rNew(i, d) = 2 * (ran3(&idum) - 0.5);

            // Updating analytic
            slater->setPosition(rNew, i);
            slater->updateMatrix();

            // Updating numeric
            slaterNumeric->setPosition(rNew, i);
            slaterNumeric->updateMatrix();

            // 50 percent chance of acceptance

            if (ran3(&idum) >= 0.5) {
                slater->updateInverse();
                slater->acceptPosition();
                slaterNumeric->updateInverseNumeric();
                slaterNumeric->acceptPosition();
                for (int d = 0; d < dim; d++)
                    rOld(i, d) = rNew(i, d);
            } else {
                slater->rejectPosition();
                slaterNumeric->rejectPosition();
                for (int d = 0; d < dim; d++)
                    rNew(i, d) = rOld(i, d);
            }

            mat DpInv = slater->DpInv;
            mat DmInv = slater->DmInv;
            mat Dp = slater->Dp;
            mat Dm = slater->Dm;
            mat DpInvNumeric = slaterNumeric->DpInv;
            mat DmInvNumeric = slaterNumeric->DmInv;

            // Printing difference
            if (abs(det(DpInv * Dp) - 1) > errorThreshold)
                std::cout << "%TEST_FAILED% time=0 testname=test2 (newSlaterTest) message=Analytic solution not matching numeric. " << endl
                    << DpInv * Dp << std::endl;

        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void test2() {
    std::cout << "newSlaterTest test 2" << std::endl;

    double errorThreshold = 1e-5;

    int dim = 2;
    int nParticles = 6;
    double alpha = 0.987;
    double w = 1.0;
    long idum = -1;
    int cycles = 1000;

    mat rInit = randu(nParticles, dim);
    mat rNew = rInit;

    //--------------------------------------------------------------------------
    // Analytic
    Orbital *orbital = new QDOrbital(dim, alpha, w);
    Slater *slater = new Slater(dim, nParticles, orbital);

    slater->setPosition(rInit, 0);
    slater->init();
    slater->acceptPosition();

    //--------------------------------------------------------------------------
    // Moving one particle

    for (int n = 0; n < cycles; n++) {
        for (int i = 0; i < nParticles; i++) {

            for (int d = 0; d < dim; d++)
                rNew(i, d) = 2 * (ran3(&idum) - 0.5);

            slater->setPosition(rNew, i);
            slater->updateMatrix();
            slater->updateInverse();
            slater->acceptPosition();

            // Sum over Laplace elements
            double laplace = 0;
            double laplaceNumeric = 0;

            for (int j = 0; j < nParticles; j++) {
                laplace += slater->getLaplacian(j);
                laplaceNumeric += slater->getLaplacianNumerical(j);
            }

            if (abs(laplace - laplaceNumeric) > errorThreshold)
                std::cout << "%TEST_FAILED% time=0 testname=test2 (newSlaterTest) message=Analytic solution not matching numeric." << endl
                    << "Differeance = " << abs(laplace - laplaceNumeric) << " "
                << "\t laplace = " << laplace << " "
                    << "\t laplaceNumeric = " << laplaceNumeric << " "
                    << std::endl;
        }
    }

}

////////////////////////////////////////////////////////////////////////////////

// Testing the Slater gradient
/*
void test3() {
    std::cout << "newSlaterTest test 3" << std::endl;

    double errorThreshold = 1e-3;

    int dim = 2;
    int nParticles = 2;
    double alpha = 0.987;
    double w = 1.0;
    long idum = -1;
    int cycles = 1000;

    mat rInit = randu(nParticles, dim);
    mat rNew = rInit;
    mat rOld = rInit;

    //--------------------------------------------------------------------------
    // Analytic
    Orbital *orbital = new QDOrbital(dim, alpha, w);
    Slater *slater = new Slater(dim, nParticles, orbital);

    slater->setPosition(rInit, 0);
    slater->init();
    slater->acceptPosition();

    //--------------------------------------------------------------------------
    // Moving one particle
    rowvec gradientSlater;
    rowvec gradientSlaterNumeric;

    for (int n = 0; n < cycles; n++) {
        for (int i = 0; i < nParticles; i++) {

            for (int d = 0; d < dim; d++)
                rNew(i, d) = 2 * (ran3(&idum) - 0.5);

            slater->setPosition(rNew, i);
            slater->updateMatrix();

            double R = slater->getRatio();

            //------------------------------------------------------------------
            // Checking before update

            for (int j = 0; j < nParticles; j++) {
                slater->computeGradient(i);
                gradientSlater = slater->getGradient() / R;
                gradientSlaterNumeric = slater->getGradientNumerical(i);

                rowvec error = gradientSlater - gradientSlaterNumeric;
                if (error.max() > errorThreshold || error.min() < -errorThreshold)
                    std::cout << "%TEST_FAILED% time=0 testname=test3 (newSlaterTest) message= Analytical gradient not matching numerical before update." << endl
                        << "\t Differeance = " << error << " "
                        << "\t anaytical = " << gradientSlater << " "
                        << "\t numerical = " << gradientSlaterNumeric << " "
                        << std::endl;

            }

            slater->updateInverse();
            slater->acceptPosition();

            //------------------------------------------------------------------
            // Checking after update
            for (int j = 0; j < nParticles; j++) {
                slater->computeGradient(i);
                gradientSlater = slater->getGradient();
                gradientSlaterNumeric = slater->getGradientNumerical(i);

                rowvec error = gradientSlater - gradientSlaterNumeric;
                if (error.max() > errorThreshold || error.min() < -errorThreshold)
                    std::cout << "%TEST_FAILED% time=0 testname=test3 (newSlaterTest) message= Analytical gradient not matching numerical after update." << endl
                        << "\t Differeance = " << error << " "
                        << "\t anaytical = " << gradientSlater << " "
                        << "\t numerical = " << gradientSlaterNumeric << " "
                        << std::endl;
            }
            //------------------------------------------------------------------
        }
    }
}
 */
////////////////////////////////////////////////////////////////////////////////

void test4() {
}

int main(int argc, char** argv) {

    std::cout << "%SUITE_STARTING% newSlaterTest" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;
    /*
        std::cout << "%TEST_STARTED% testInverse (newSlaterTest)" << std::endl;
        testingInverse();
        std::cout << "%TEST_FINISHED% time=0 testInverse (newSlaterTest)" << std::endl;
    
            std::cout << "%TEST_STARTED% test2 (newSlaterTest)\n" << std::endl;
            test2();
            std::cout << "%TEST_FINISHED% time=0 test2 (newSlaterTest)" << std::endl;

            std::cout << "%TEST_STARTED% test3 (newSlaterTest)\n" << std::endl;
            test3();
            std::cout << "%TEST_FINISHED% time=0 test3 (newSlaterTest)" << std::endl;
     **/

    std::cout << "%TEST_STARTED% test3 (newSlaterTest)\n" << std::endl;
    test4();
    std::cout << "%TEST_FINISHED% time=0 test3 (newSlaterTest)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

