/* 
 * File:   ComputeGrid.cpp
 * Author: Sigve
 * 
 * Created on May 10, 2012, 6:30 PM
 */

#include "ComputeGrid.h"

#include "WaveFunction.h"
#include "QD/QDOrbital.h"
#include "QD/QDJastrow.h"
#include "QD/QDHamiltonian.h"
#include "includes/ini.h"

#include <cstdlib>
#include <mpi.h>

////////////////////////////////////////////////////////////////////////////////

ComputeGrid::ComputeGrid() {
    int accepted_tmp;
    double localE;
    long idum;
    bool writeToFile = false;
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

        McSamples = (int) INIreader.GetDouble("ComputeGrid", "McSamples");
        importanceSampling = INIreader.GetBool("ComputeGrid", "importanceSampling");
        thermalization = (int) INIreader.GetDouble("ComputeGrid", "thermalization");

        dim = INIreader.GetInt("ComputeGrid", "dim");
        alphaStart = INIreader.GetDouble("ComputeGrid", "alphaStart");
        betaStart = INIreader.GetDouble("ComputeGrid", "betaStart");
        deltaAlpha = INIreader.GetDouble("ComputeGrid", "deltaAlpha");
        deltaBeta = INIreader.GetDouble("ComputeGrid", "deltaBeta");
        w = INIreader.GetDouble("ComputeGrid", "w");
        nParticles = INIreader.GetInt("ComputeGrid", "nParticles");
        nAlpha = INIreader.GetInt("ComputeGrid", "nAlpha");
        nBeta = INIreader.GetInt("ComputeGrid", "nAlpha");

        usingJastrow = INIreader.GetBool("ComputeGrid", "usingJastrow");

        writeToFile = INIreader.GetBool("ComputeGrid", "writeToFile");
        fileName = INIreader.GetString("ComputeGrid", "fileName");
    }

    //--------------------------------------------------------------------------
    // Broadcasting data to all nodes.
    if (myRank == 0)
        cout << "Broadcasting data to nodes." << endl;

    MPI_Bcast(&McSamples, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&importanceSampling, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&thermalization, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&alphaStart, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&betaStart, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&deltaAlpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&deltaBeta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&w, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nParticles, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nAlpha, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nBeta, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&usingJastrow, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //--------------------------------------------------------------------------
    // Computing grid of energy values.
    if (myRank == 0)
        cout
            << "Computing grid with " << endl
            << "\t McSamples = " << McSamples
            << "\t alphaStart = " << alphaStart
            << "\t betaStart = " << betaStart
            << "\t deltaAlpha = " << deltaAlpha
            << "\t deltaBeta = " << deltaBeta
            << "\t w = " << w
            << endl;

    E = zeros(nAlpha, nBeta);
    Esq = zeros(nAlpha, nBeta);

    for (int k = 0; k < nAlpha; k++) {
        alpha = alphaStart + k*deltaAlpha;
        beta = betaStart;
        for (int l = 0; l < nBeta; l++) {
            beta += deltaBeta;

            //------------------------------------------------------------------
            // Configuring the wave function.    

            Orbital *orbital = new QDOrbital(dim, alpha, w);
            Jastrow * jastrow = NULL;
            if (usingJastrow)
                jastrow = new QDJastrow(dim, nParticles, beta);

            Hamiltonian *hamiltonian = new QDHamiltonian(dim, nParticles, w, usingJastrow);
            WaveFunction wf(dim, nParticles, idum, orbital, jastrow, hamiltonian);

            // If we are using brute force sampling we have to calculate the optimal step length
            if (!importanceSampling)
                wf.setOptimalStepLength();

            //------------------------------------------------------------------
            // Monte Carlo loop.  

            accepted = 0;

            for (int i = 0; i < McSamples + thermalization; i++) {
                for (int j = 0; j < nParticles; j++) {
                    if (importanceSampling)
                        accepted_tmp = wf.tryNewPosition(j);
                    else
                        accepted_tmp = wf.tryNewPositionBF(j);

                    if (i > thermalization) {
                        localE = wf.sampleEnergy();
                        E(k, l) += localE;
                        Esq(k, l) += localE*localE;
                        accepted += accepted_tmp;
                    }
                }
            }

            //------------------------------------------------------------------
            // Collecting and scaling results from all nodes.
            MPI_Allreduce(MPI_IN_PLACE, &E(k, l), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &Esq(k, l), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &accepted, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            E(k, l) /= (McSamples * nParticles * nNodes);
            Esq(k, l) /= (McSamples * nParticles * nNodes);
            accepted /= (McSamples * nParticles * nNodes);

            //------------------------------------------------------------------
            // Printing results.
            if (myRank == 0)
                cout
                    << "\talpha = " << alpha
                    << "\tbeta = " << beta
                    << "\tE = " << E(k, l)
                << "\tVariance = " << Esq(k, l) - E(k, l) * E(k, l)
                << "\tAccepted = " << accepted
                    << endl;
        }
    }

    //--------------------------------------------------------------------------
    // Writing results to file.
    if (myRank == 0) {
        if (writeToFile) {
            cout << "Writing data to file." << endl;
            writeEnergyToFile();
        }
    }
    //--------------------------------------------------------------------------
}
////////////////////////////////////////////////////////////////////////////////

ComputeGrid::ComputeGrid(const ComputeGrid& orig) {
}

////////////////////////////////////////////////////////////////////////////////

ComputeGrid::~ComputeGrid() {
}

////////////////////////////////////////////////////////////////////////////////

void ComputeGrid::writeEnergyToFile() {
    double a, b;

    ofstream outStream;
    outStream.open((const char*) &fileName[0]);
    outStream << nParticles << " " << nAlpha << " " << nBeta << " " << McSamples << endl;

    for (int i = 0; i < nAlpha; i++) {
        a = alphaStart + i*deltaAlpha;

        for (int j = 0; j < nBeta; j++) {
            b = betaStart + j*deltaBeta;
            outStream << a << " " << b << " " << E(i,j) << " " << Esq(i,j) << "\n";
        }
    }
    outStream.close();
}