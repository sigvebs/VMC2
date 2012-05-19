/* 
 * File:   ComputeOneWF.cpp
 * Author: Sigve
 * 
 * Created on May 9, 2012, 10:51 PM
 */

#include "ComputeOneWF.h"

#include "WaveFunction.h"
#include "QD/QDOrbital.h"
#include "QD/QDJastrow.h"
#include "QD/QDHamiltonian.h"
#include "includes/ini.h"

#include <cstdlib>
#include <mpi.h>
#include <fstream>
#include <iostream>

////////////////////////////////////////////////////////////////////////////////

ComputeOneWF::ComputeOneWF(const ComputeOneWF& orig) {
}

////////////////////////////////////////////////////////////////////////////////

ComputeOneWF::~ComputeOneWF() {
}

////////////////////////////////////////////////////////////////////////////////

ComputeOneWF::ComputeOneWF() {
    int accepted_tmp;
    double localE;
    long idum;
    //--------------------------------------------------------------------------
    // MPI
    int myRank, nNodes;
    MPI_Comm_size(MPI_COMM_WORLD, &nNodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    idum = -1 - myRank - time(NULL);

    //--------------------------------------------------------------------------
    // Reading data from QD.ini
    if (myRank == 0) {
        cout << "Reading configuration from file " << endl;
        ini INIreader("QD.ini");

        McSamples = (int) INIreader.GetDouble("ComputeOneWF", "McSamples");
        importanceSampling = INIreader.GetBool("ComputeOneWF", "importanceSampling");
        thermalization = (int) INIreader.GetDouble("ComputeOneWF", "thermalization");

        dim = INIreader.GetInt("ComputeOneWF", "dim");
        alpha = INIreader.GetDouble("ComputeOneWF", "alpha");
        beta = INIreader.GetDouble("ComputeOneWF", "beta");
        w = INIreader.GetDouble("ComputeOneWF", "w");
        nParticles = INIreader.GetInt("ComputeOneWF", "nParticles");

        usingJastrow = INIreader.GetBool("ComputeOneWF", "usingJastrow");

        // Checking whether we have to store every local energy for
        // blocking analysis.
        if ((int) INIreader.GetInt("main", "option") == 4)
            blocking = true;
        else
            blocking = false;
    }
    //--------------------------------------------------------------------------
    // Broadcasting data to all nodes.
    if (myRank == 0)
        cout << "Broadcasting data to nodes." << endl;

    MPI_Bcast(&McSamples, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&importanceSampling, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&thermalization, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&beta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&w, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nParticles, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&usingJastrow, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&blocking, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //--------------------------------------------------------------------------
    // Configuring the wave function.    
    if (myRank == 0)
        cout << "Configuring the wave function." << endl;

    Orbital *orbital = new QDOrbital(dim, alpha, w);
    Jastrow *jastrow = NULL;
    if (usingJastrow)
        jastrow = new QDJastrow(dim, nParticles, beta);

    Hamiltonian *hamiltonian = new QDHamiltonian(dim, nParticles, w, usingJastrow);
    WaveFunction wf(dim, nParticles, idum, orbital, jastrow, hamiltonian);

    // If we are using brute force sampling we have to calculate the optimal step length
    if (!importanceSampling)
        wf.setOptimalStepLength();
    //--------------------------------------------------------------------------
    // If Blocking
    ostringstream filename;
    filename << "DATA/Blocking/blocking_" << myRank << ".dat";
    ofstream blockfile(filename.str().c_str(), ios::out | ios::binary);
    //--------------------------------------------------------------------------
    // Monte Carlo loop.    
    if (myRank == 0)
        cout
            << "Starting Monte Carlo sampling with " << endl
            << "\t McSamples = " << McSamples
            << "\t alpha = " << alpha
            << "\t beta = " << beta
            << "\t w = " << w
            << endl;

    accepted = 0;
    E = 0;
    Esq = 0;
    for (int i = 0; i < McSamples + thermalization; i++) {
        for (int j = 0; j < nParticles; j++) {
            
            if (importanceSampling) {
                accepted_tmp = wf.tryNewPosition(j);
            } else {
                accepted_tmp = wf.tryNewPositionBF(j);
            }
            
            if (i > thermalization) {
                localE = wf.sampleEnergy();
                E += localE;
                Esq += localE*localE;
                accepted += accepted_tmp;

                if (blocking)
                    blockfile.write((char*) &localE, sizeof (double));
            }
        }
    }

    //if (blocking)
    blockfile.close();

    //--------------------------------------------------------------------------
    // Collecting and scaling results from all nodes.
    if (myRank == 0)
        cout << "Collecting and scaling results from all nodes." << endl;

    MPI_Allreduce(MPI_IN_PLACE, &E, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Esq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &accepted, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    E /= (McSamples * nParticles * nNodes);
    Esq /= (McSamples * nParticles * nNodes);
    accepted /= (McSamples * nParticles * nNodes);

    //--------------------------------------------------------------------------
    // Printing results.
    if (myRank == 0)
        cout
            << "\tE = " << E
            << "\tVariance = " << Esq - E * E
            << "\tAccepted = " << accepted
            << endl;
}
////////////////////////////////////////////////////////////////////////////////
