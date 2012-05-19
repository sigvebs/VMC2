/* 
 * File:   SGD.cpp
 * Author: Sigve
 * 
 * Created on May 10, 2012, 7:42 PM
 */

#include "SGD.h"

#include "WaveFunction.h"
#include "QD/QDOrbital.h"
#include "QD/QDJastrow.h"
#include "QD/QDHamiltonian.h"
#include "includes/ini.h"

#include <cstdlib>
#include <mpi.h>

SGD::SGD() {
    double tPrev, t, step;
    rowvec parameter(2), gradient(2), gradientLocal(2), gradientOld(2), nVar(2);
    double localE;
    int correlationLength;
    long idum;
    writeToFile = false;

    //--------------------------------------------------------------------------
    // MPI
    int myRank, nNodes;
    MPI_Comm_size(MPI_COMM_WORLD, &nNodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    //--------------------------------------------------------------------------
    // Reading data from QD.ini
    if (myRank == 0) {
        cout << "Reading configuration from file " << endl;
        ini INIreader("QD.ini");

        McSamples = (int) INIreader.GetDouble("SGD", "McSamples");
        importanceSampling = INIreader.GetBool("SGD", "importanceSampling");
        thermalization = (int) INIreader.GetDouble("SGD", "thermalization");
        SGDSamples = INIreader.GetInt("SGD", "SGDSamples");

        dim = INIreader.GetInt("SGD", "dim");
        alpha = INIreader.GetDouble("SGD", "alpha");
        beta = INIreader.GetDouble("SGD", "beta");
        w = INIreader.GetDouble("SGD", "w");
        nParticles = INIreader.GetInt("SGD", "nParticles");
        m = INIreader.GetInt("SGD", "m");
        correlationLength = INIreader.GetInt("SGD", "correlationLength");

        usingJastrow = INIreader.GetBool("SGD", "usingJastrow");
        writeToFile = INIreader.GetBool("SGD", "writeToFile");
        fileName = INIreader.GetString("SGD", "fileName");

        fMax = INIreader.GetDouble("SGD", "fMax");
        fMin = INIreader.GetDouble("SGD", "fMin");
        omega = INIreader.GetDouble("SGD", "omega");
        expo = INIreader.GetDouble("SGD", "expo");
        A = INIreader.GetDouble("SGD", "A");
        a = INIreader.GetDouble("SGD", "a");
        maxStep = INIreader.GetDouble("SGD", "maxStep");
    }

    //--------------------------------------------------------------------------
    // Broadcasting data to all nodes.
    if (myRank == 0)
        cout << "Broadcasting data to nodes." << endl;

    MPI_Bcast(&McSamples, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&SGDSamples, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&importanceSampling, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&thermalization, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&beta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&w, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nParticles, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&usingJastrow, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&correlationLength, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(&fMax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&fMin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&omega, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&expo, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&A, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&a, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&maxStep, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#if 0
    cout
            << "\t A = " << A
            << "\t a = " << a
            << "\t maxStep = " << maxStep
            << "\t expo = " << expo
            << "\t omega = " << omega
            << "\t fMin = " << fMin
            << "\t fMax = " << fMax
            << endl;
#endif
    //--------------------------------------------------------------------------
    //Initializing and thermalizing walkers.

    if (myRank == 0)
        cout << "Initializing and thermalizing walkers." << endl;

    WaveFunction * walker[m];

    idum = -1 - myRank - time(NULL);

    Orbital *orbital = new QDOrbital(dim, alpha, w);
    Jastrow * jastrow = NULL;
    if (usingJastrow)
        jastrow = new QDJastrow(dim, nParticles, beta);

    Hamiltonian *hamiltonian = new QDHamiltonian(dim, nParticles, w, usingJastrow);
    WaveFunction *wf = new WaveFunction(dim, nParticles, idum, orbital, jastrow, hamiltonian);

    // If we are using brute force sampling we have to calculate 
    // the optimal step length.
    if (!importanceSampling)
        wf->setOptimalStepLength();


    // Thermalizing walker.
    int k = 0;
    int thermCorr = m*correlationLength;

    for (int i = 0; i <= thermalization + thermCorr; i++) {
        for (int j = 0; j < nParticles; j++) {
            if (importanceSampling)
                wf->tryNewPosition(j);
            else
                wf->tryNewPositionBF(j);
        }

        // Storing walker.
        if (i > thermalization) {
            if ((i - thermalization) % (correlationLength) == 0) {
                if (myRank == 0)
                    cout << "k = " << k << endl;
                walker[k] = wf->clone();
                k++;
            }
        }
    }

    /*
     for (int k = 0; k < m; k++) {
            // Giving each walker a unique seed.
            idum = -1 - myRank - time(NULL) - k;
        
            Orbital *orbital = new QDOrbital(dim, alpha, w);
            Jastrow * jastrow = NULL;
            if (usingJastrow)
                jastrow = new QDJastrow(dim, nParticles, beta);

            Hamiltonian *hamiltonian = new QDHamiltonian(dim, nParticles, w, usingJastrow);
            WaveFunction *wf = new WaveFunction(dim, nParticles, idum, orbital, jastrow, hamiltonian);

            // If we are using brute force sampling we have to calculate 
            // the optimal step length.
            if (!importanceSampling)
                wf->setOptimalStepLength();

            // Thermalizing walker.
            for (int i = 0; i < thermalization; i++) {
                for (int j = 0; j < nParticles; j++) {
                    if (importanceSampling)
                        wf->tryNewPosition(j);
                    else
                        wf->tryNewPositionBF(j);
                }
            }

            // Storing walker.
            walker[k] = wf;
        }
     **/
    //--------------------------------------------------------------------------
    ofstream outStream;
    if (myRank == 0) {
        cout << "Starting minimizing process " << endl;
        cout << "SGA samples = " << SGDSamples;
        cout << "\t omega = " << w;
        cout << "\t alpha_0 = " << alpha;
        cout << "\t beta_0 = " << beta;

        outStream.open((const char*) &fileName[0]);
    }
    //--------------------------------------------------------------------------
    tPrev = 0;
    parameter(0) = alpha;
    parameter(1) = beta;
    nVar(0) = 5;
    nVar(1) = 5;
    gradientOld = zeros(1, 2);
    //-------------------------------------------------------------------------- 
    for (int sample = 1; sample <= SGDSamples; sample++) {
 
        // Moving walkers        
        E = 0;
        accepted = 0;
        gradient = zeros(1, 2);
        gradientLocal = zeros(1, 2);

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < McSamples; j++) {
                for (int j = 0; j < nParticles; j++) {

                    if (importanceSampling)
                        accepted += walker[i]->tryNewPosition(j);
                    else
                        accepted += walker[i]->tryNewPositionBF(j);

                    localE = walker[i]->sampleEnergy();
                    E += localE / nParticles;
                    gradient += walker[i]->getVariationGradient();
                    gradientLocal += walker[i]->getVariationGradient() * localE;
                }
            }
        }

        //----------------------------------------------------------------------
        // Collecting results.
        MPI_Allreduce(MPI_IN_PLACE, &E, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &gradient(0), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &gradient(1), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &gradientLocal(0), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &gradientLocal(1), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        E /= (nNodes * m * McSamples);
        gradient /= (nNodes * m * McSamples);
        gradientLocal /= (nNodes * m * McSamples);

        // The total variational gradient
        gradient = 2 * (gradientLocal - gradient * E);
        //----------------------------------------------------------------------
        double x = -dot(gradient, gradientOld);
        f = fMin + (fMax - fMin) / (1 - (fMax / fMin) * exp(-x / omega));
        t = tPrev + f;

        if (t < 0)
            t = 0;
        //----------------------------------------------------------------------
        // Printing progress

        if (myRank == 0) {
            cout << "\r"
                    << "\t SGA cycle = " << sample
                    << "\t alpha = " << parameter(0)
                    << "\t beta = " << parameter(1)
                    << "\t E = " << E
                    << "\t t = " << t
                    << "\t step = " << step
                    //<< "\xd";
                    << endl;
            //fflush(stdout);
        }

        //----------------------------------------------------------------------
        for (int i = 0; i < 2; i++) {
            if (gradient(i) * gradientOld(i) < 0)
                nVar(i)++;

            // Calulating new step lengths.
#if 0
            step = gradient(i) * a / pow(nVar(i), expo);
            if (fabs(step) > maxStep)
                step *= maxStep / fabs(step);

#else
            // Calculating new step lengths.
            double gamma = a / pow(t + A, expo);
            step = gamma * gradient(i);
            if (fabs(step) > maxStep)
                step *= maxStep / fabs(step);
#endif
            parameter(i) -= step;
        }

        //----------------------------------------------------------------------
        // Setting the new parameters
        for (int i = 0; i < m; i++)
            walker[i]->setNewVariationalParameters(parameter(0), parameter(1));

        // Writing to file
        if (myRank == 0)
            outStream << parameter(0) << " " << parameter(1) << endl;
        //----------------------------------------------------------------------
        // Storing previous values.
        tPrev = t;
        gradientOld = gradient;
    }

    if (myRank == 0)
        outStream.close();
    //--------------------------------------------------------------------------
    // Cleaning up
    //for (int i = 0; i < m; i++)
        //delete walker[i];
    
}

////////////////////////////////////////////////////////////////////////////////

SGD::SGD(const SGD& orig) {
}

////////////////////////////////////////////////////////////////////////////////

SGD::~SGD() {
}

////////////////////////////////////////////////////////////////////////////////

