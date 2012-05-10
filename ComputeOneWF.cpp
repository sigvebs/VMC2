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

#include <cstdlib>
#include <mpi.h>
#include "includes/ini.h"

////////////////////////////////////////////////////////////////////////////////

ComputeOneWF::ComputeOneWF(const ComputeOneWF& orig) {
}

////////////////////////////////////////////////////////////////////////////////

ComputeOneWF::~ComputeOneWF() {
}

////////////////////////////////////////////////////////////////////////////////

ComputeOneWF::ComputeOneWF() {
    double localE;
    double alpha, beta, w;
    int dim;
    long idum;
    bool usingJastrow;

    //--------------------------------------------------------------------------
    // MPI
    int myRank, nNodes;
    MPI_Comm_size(MPI_COMM_WORLD, &nNodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    //--------------------------------------------------------------------------
    // Parameters.
    int McSamples;
    bool importanceSampling;
    int thermalization;
    int nParticles;

    //--------------------------------------------------------------------------
    McSamples = 1e5;
    importanceSampling = true;
    thermalization = 1e4;

    dim = 2;
    alpha = 0.987;
    alpha = 1.0;
    beta = 0.4;
    w = 1.0;
    idum = -1 - myRank - time(NULL);
    nParticles = 2;

    usingJastrow = false;
    //--------------------------------------------------------------------------
    E = 0;
    Esq = 0;
    accepted = 0;
    int accepted_tmp;

    //--------------------------------------------------------------------------
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
    
    //--------------------------------------------------------------------------
    // Monte Carlo loop.
    if (myRank == 0)
        cout << "Starting Monte Carlo sampling" << endl;

    for (int i = 0; i < McSamples + thermalization; i++) {
        for (int j = 0; j < nParticles; j++) {
            if (importanceSampling)
                accepted_tmp = wf.tryNewPosition(j);
            else
                accepted_tmp = wf.tryNewPositionBF(j);

            if (i > thermalization) {
                localE = wf.sampleEnergy();
                E += localE;
                Esq += localE*localE;
                accepted += accepted_tmp;
            }
        }
    }

    //--------------------------------------------------------------------------
    // Scaling results
    E /= (McSamples * nParticles);
    Esq /= (McSamples * nParticles);
    accepted /= (McSamples * nParticles);

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
