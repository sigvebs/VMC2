/* 
 * File:   OneBodyDensity.cpp
 * Author: sigve
 * 
 * Created on 11. mai 2012, 17:01
 */

#include "OneBodyDensity.h"

OneBodyDensity::OneBodyDensity(const OneBodyDensity& orig) {
}

OneBodyDensity::~OneBodyDensity() {
}

#include <iostream>
#include <fstream>
#include <string.h>
#include <mpi.h>
#include <time.h>

#include "includes/lib.h"
#include "includes/ini.h"

#include "WaveFunction.h"
#include "QD/QDHamiltonian.h"
#include "QD/QDJastrow.h"
#include "QD/QDOrbital.h"

////////////////////////////////////////////////////////////////////////////////

OneBodyDensity::OneBodyDensity() {
    //--------------------------------------------------------------------------
    // MPI
    MPI_Comm_size(MPI_COMM_WORLD, &nNodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    //--------------------------------------------------------------------------
    double alpha, beta, omega, deltaStep;
    bool usingJastrow;
    int nSteps;

    //--------------------------------------------------------------------------
    // Reading data from QD.ini

    if (myRank == 0) {
        cout << "Reading configuration from file " << endl;
        ini INIreader("QD.ini");

        McSamples = (int) INIreader.GetDouble("OneBodyDensity", "McSamples");
        nSteps = INIreader.GetInt("OneBodyDensity", "nSteps");
        deltaStep = INIreader.GetDouble("OneBodyDensity", "deltaStep");

        alpha = INIreader.GetDouble("OneBodyDensity", "alpha");
        beta = INIreader.GetDouble("OneBodyDensity", "beta");
        omega = INIreader.GetDouble("OneBodyDensity", "omega");

        dim = INIreader.GetInt("OneBodyDensity", "dim");
        nParticles = INIreader.GetInt("OneBodyDensity", "nParticles");
        usingJastrow = INIreader.GetBool("OneBodyDensity", "usingJastrow");

        fileName = INIreader.GetString("OneBodyDensity", "fileName");
    }

    MPI_Bcast(&McSamples, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nSteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nParticles, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&usingJastrow, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(&deltaStep, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&beta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&omega, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    idum = idum - myRank - time(NULL);
    McSamples /= nNodes;
    
    //--------------------------------------------------------------------------
    // Configuring the wave function.    
    if (myRank == 0)
        cout << "Configuring the wave function." << endl;

    Orbital *orbital = new QDOrbital(dim, alpha, omega);
    Jastrow *jastrow = NULL;
    if (usingJastrow)
        jastrow = new QDJastrow(dim, nParticles, beta);

    Hamiltonian *hamiltonian = new QDHamiltonian(dim, nParticles, omega, usingJastrow);
    wf = new WaveFunction(dim, nParticles, idum, orbital, jastrow, hamiltonian);

    //--------------------------------------------------------------------------
    // One body density 
    if (myRank == 0)
        cout << "Calculating one body density." << endl;

    mat rho = zeros(nSteps, nSteps);
    vec r = zeros<vec > (dim);

    int count = 1;
    for (int i = 0; i < nSteps; i++) {
        r(0) = (i - nSteps / 2) * deltaStep;
        for (int j = 0; j < nSteps; j++) {
            r(1) = (j - nSteps / 2) * deltaStep;
            rho(i, j) = McIntegrator(r, count);
            count++;
        }
        
        if (myRank == 0)
            cout << i << endl;
    }

    //--------------------------------------------------------------------------
    if (myRank == 0) {
        rho = normalize(rho);
        writeToFile(rho, deltaStep);
    }
}


////////////////////////////////////////////////////////////////////////////////

mat OneBodyDensity::normalize(const mat& rho) {
    double norm = sum(sum(rho));
    return rho / norm;
}

////////////////////////////////////////////////////////////////////////////////

double OneBodyDensity::McIntegrator(const vec &rActive, int id) {
    double stepLength = 6;
    double rho = 0;
    double waveFunc;
    mat r = zeros(nParticles, dim);
    r.row(0) = rActive.t();

    //--------------------------------------------------------------------------
    // Init of the particle that is at a constant r
    for (int i = 0; i < McSamples; i++) {
        // Generating new coordintes.
        for (int k = 1; k < nParticles; k++)
            for (int j = 0; j < dim; j++)
                r(k, j) = stepLength * (ran2(&idum) - 0.5);

        waveFunc = wf->evaluate(r);
        rho += waveFunc*waveFunc;
    }

    //--------------------------------------------------------------------------

    double tmp = rho;

    if (myRank != 0)
        MPI_Send(&rho, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        // Collecting data from the other processes.
    else {
        for (int i = 1; i < nNodes; i++) {
            MPI_Recv(&tmp, 1, MPI_DOUBLE, i, id, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            rho += tmp;
        }
    }

    rho /= nNodes;
    return rho;
}

////////////////////////////////////////////////////////////////////////////////

void OneBodyDensity::writeToFile(const mat &rho, double deltaStep) {
    double nSteps = rho.n_rows;

    ofstream outStream;
    outStream.open((const char*) &fileName[0]);
    outStream << nSteps << " " << deltaStep << " " << McSamples << endl;

    // Writing coordinates and density to file
    vec r = zeros(dim, 1);
    for (int i = 0; i < nSteps; i++) {
        r(0) = (i - nSteps / 2) * deltaStep;
        for (int j = 0; j < nSteps; j++) {
            r(1) = (j - nSteps / 2) * deltaStep;
            outStream << r(0) << " " << r(1) << " " << rho(i, j) << endl;
        }
    }
    outStream.close();
}
