/* 
 * File:   DMC.cpp
 * Author: Sigve
 * 
 * Created on May 12, 2012, 1:33 PM
 */

#include "DMC.h"

#include "WaveFunction.h"
#include "QD/QDOrbital.h"
#include "QD/QDJastrow.h"
#include "QD/QDHamiltonian.h"
#include "includes/ini.h"
#include "includes/lib.h"

#include <cstdlib>
#include <mpi.h>

#include <list>
#include <vector>
using namespace std;

////////////////////////////////////////////////////////////////////////////////

DMC::DMC() {
    long idum;
    string fileName;
    bool writeToFile = false;
    int DMCThermal, blockSize, correlationLength;
    //--------------------------------------------------------------------------
    // MPI
    int myRank, nNodes;
    MPI_Comm_size(MPI_COMM_WORLD, &nNodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    if (myRank == 0)
        cout << "Starting DMC." << endl;

    //--------------------------------------------------------------------------
    // Reading data from QD.ini
    if (myRank == 0) {
        cout << "Reading configuration from file." << endl;
        ini INIreader("QD.ini");

        tau = INIreader.GetDouble("DMC", "tau");
        McSamples = (int) INIreader.GetDouble("DMC", "McSamples");
        importanceSampling = INIreader.GetBool("DMC", "importanceSampling");
        thermalization = (int) INIreader.GetDouble("DMC", "thermalization");
        nWalkers = (int) INIreader.GetDouble("DMC", "nWalkers");
        nSteps = (int) INIreader.GetDouble("DMC", "nSteps");
        DMCThermal = (int) INIreader.GetDouble("DMC", "DMCThermal");
        blockSize = (int) INIreader.GetDouble("DMC", "blockSize");
        correlationLength = (int) INIreader.GetDouble("DMC", "correlationLength");

        dim = INIreader.GetInt("DMC", "dim");
        alpha = INIreader.GetDouble("DMC", "alpha");
        beta = INIreader.GetDouble("DMC", "beta");
        w = INIreader.GetDouble("DMC", "w");
        nParticles = INIreader.GetInt("DMC", "nParticles");

        usingJastrow = INIreader.GetBool("DMC", "usingJastrow");

        writeToFile = INIreader.GetBool("DMC", "writeToFile");
        fileName = INIreader.GetString("DMC", "fileName");
    }

    //--------------------------------------------------------------------------
    // Debug printing
    if (myRank == 66) {
        cout
                << "tau = " << tau << endl
                << "McSamples = " << McSamples << endl
                << "importanceSampling = " << importanceSampling << endl
                << "thermalization = " << thermalization << endl
                << "nWalkers = " << nWalkers << endl
                << "nSteps = " << nSteps << endl
                << "DMCThermal = " << DMCThermal << endl
                << "blockSize = " << blockSize << endl
                << "correlationLength = " << correlationLength << endl
                << "dim = " << dim << endl
                << "alpha = " << alpha << endl
                << "beta = " << beta << endl
                << "w = " << w << endl
                << "nParticles = " << nParticles << endl
                << "usingJastrow = " << usingJastrow << endl
                << "writeToFile = " << writeToFile << endl
                << "fileName = " << fileName << endl
                << endl;
    }
    int maxWalkers = nWalkers * 1.5;
    //--------------------------------------------------------------------------
    // Configuring the wave function.    
    //--------------------------------------------------------------------------
    if (myRank == 0) {
        cout << "Configuring the wave function." << endl;
    }

    idum = -1 - myRank - time(NULL);

    Orbital *orbital = new QDOrbital(dim, alpha, w);
    Jastrow *jastrow = new QDJastrow(dim, nParticles, beta);
    Hamiltonian *hamiltonian = new QDHamiltonian(dim, nParticles, w, usingJastrow);
    WaveFunction wf(dim, nParticles, -1 - myRank - time(NULL), orbital, jastrow, hamiltonian);

    //--------------------------------------------------------------------------
    // Creating walkers from one VMC run. At every correlation length the 
    // WF is copied into a new walker. This way we get a good sampling of
    // the WF distribution.
    //--------------------------------------------------------------------------
    if (myRank == 0) {
        cout
                << "Initializing walkers using VMC,"
                << "\t thermalization = " << thermalization
                << "\t correlationLength = " << correlationLength
                << "\t nWalkers = " << nWalkers
                << endl;
    }

    vector<WaveFunction*> walkers;
    McSamples = correlationLength*nWalkers;
    double ETrial = 0;
    int cWalker = 0;

    for (int i = 0; i < thermalization + McSamples; i++) {

        // Moving walker.
        for (int j = 0; j < nParticles; j++) {
            wf.tryNewPosition(j);
        }

        if (i >= thermalization) {
            ETrial += wf.sampleEnergy();

            // Cloning a walker
            if (i % (correlationLength) == 0) {
                walkers.push_back(wf.clone());
                cWalker++;
                if (myRank == 66) {
                    cout << "Cloning walker " << cWalker << endl;
                }

            }

        }
    }

    ETrial /= McSamples;
    if (myRank == 0) {
        cout
                << walkers.size() << " walkers created" << endl
                << "Trial energy after VMC: E = " << ETrial << endl;
    }
    //--------------------------------------------------------------------------
    // Thermalizing walkers using DMC. The trial energy is updated each time 
    // the walkers are moved. The bad walkers are killed. The good ones get
    // cloned.
    //--------------------------------------------------------------------------
    if (myRank == 0) {
        cout << "Starting DMC thermalization" << endl;
    }

    WaveFunction *toBeKilled;
    int walkersLength, NSamples, copies;
    double EOld, ENew, ETmp, EAvg, bFactor;

    for (int n = 1; n <= DMCThermal; n++) {
        ETmp = 0;
        NSamples = 0;

        // Looping over every walker in walkers once.
        walkersLength = walkers.size();
        for (int i = 0; i < walkersLength; i++) {

            EOld = walkers[i]->sampleEnergy();

            // Moving walker.
            for (int j = 0; j < nParticles; j++) {
                walkers[i]->DMCtryNewPosition(j);

                ENew = walkers[i]->sampleEnergy();

                // Calculating the branching factor.
                bFactor = exp(-tau * (0.5 * (ENew + EOld) - ETrial));

                //------------------------------------------------------------------
                copies = int(bFactor + ran3(&idum));

                // The walker is killed.
                if (copies == 0) {
                    toBeKilled = walkers[i];
                    walkers.erase(walkers.begin() + i);
                    delete toBeKilled;
                    walkersLength--;
                    i--;
                } else {
                    // Cloning walker.
                    if (walkers.size() < maxWalkers) {
                        for (int j = 1; j < copies; j++) {
                            walkers.push_back(walkers[i]->clone());
                        }
                    }
                    // Sampling energy.
                    ETmp += bFactor*ENew;
                    NSamples++;
                }
            }
            //------------------------------------------------------------------
        }
        //----------------------------------------------------------------------
        // Updating the new trial energy and printin progress.
        ETmp /= NSamples;
        EAvg += ETmp;
        ETrial = EAvg / n;

        if (myRank == 0) {
            cout
                    << "TERMALIZATION: "
                    << "\tn = " << n
                    << "\tnWalkers = " << walkers.size()
                    << "\t ETmp = " << ETmp
                    << "\t EAvg = " << EAvg / n
                    << "\xd";
        }
    }

    //--------------------------------------------------------------------------
    // DMC sampling. The trial energy is updated after every walker has moved
    // a block. 
    //--------------------------------------------------------------------------    
    ofstream outStream;
    if (myRank == 0) {
        cout << "Starting DMC sampling" << endl;
        outStream.open((const char*) &fileName[0]);
    }

    int NTotSamples;
    EAvg = 0;

    double EBlock;
    for (int n = 1; n <= DMCThermal; n++) {

        // Looping over every walker in walkers once.
        walkersLength = walkers.size();
        int WL = walkersLength;
        NSamples = 0;
        NTotSamples = 0;
        EBlock = 0;
        int count = 0;
        for (int i = 0; i < walkersLength; i++) {

            // Moving a walker a block
            for (int b = 0; b < blockSize; b++) {
                EOld = walkers[i]->sampleEnergy();

                // Moving walker.
                for (int j = 0; j < nParticles; j++) {
                    walkers[i]->DMCtryNewPosition(j);


                    ENew = walkers[i]->sampleEnergy();

                    // Calculating the branching factor
                    bFactor = exp(-tau * (0.5 * (ENew + EOld) - ETrial));

                    //--------------------------------------------------------------
                    copies = int(bFactor + ran3(&idum));

                    // The walker is killed.
                    if (copies == 0) {
                        toBeKilled = walkers[i];
                        walkers.erase(walkers.begin() + i);
                        delete toBeKilled;
                        walkersLength--;
                        i--;
                        b = blockSize;
                        break;
                    } else {
                        // Cloning walker.
                        if (walkers.size() < maxWalkers) {
                            for (int j = 1; j < copies; j++) {
                                walkers.push_back(walkers[i]->clone());
                            }
                        }
                        // Sampling energy.
                        EBlock += bFactor*ENew;
                        NSamples++;
                    }
                }
                //--------------------------------------------------------------
            }
        }

        //----------------------------------------------------------------------
        // Updating the new trial energy and printing progress.
        EBlock /= NSamples;
        EAvg += EBlock;
        ETrial = EAvg / n;

        if (myRank == 0) {
            cout
                    << "DMC:  = " << n
                    << "\tnWalkers = " << walkers.size()
                    << "\t ETmp = " << EBlock
                    << "\t EAvg = " << EAvg / n
                    << "\t Max samples = " << WL * blockSize
                    << "\t Samples = " << NSamples
                    << "\t Diff = " << WL * blockSize - NSamples
                    << endl;

            // Writing to file
            if (myRank == 0) {
                outStream << EBlock << endl;
            }
        }
    }

    //--------------------------------------------------------------------------
    if (myRank == 0)
        outStream.close();

}

////////////////////////////////////////////////////////////////////////////////

DMC::DMC(const DMC & orig) {
}

////////////////////////////////////////////////////////////////////////////////

DMC::~DMC() {
}

// Renormalizing the number of walkers.
//----------------------------------------------------------------------
/*
        // Renormalizing the number of walkers.
        if (walker.size() > nWalkers) {
            // Kill random walkers until we have nWalkers.
            int nDeathrow = walker.size() - nWalkers;

            for (int i = 0; i < nDeathrow; i++)
                walker.erase(walker.begin() + int(ran3(&idum) * walker.size()));

        } else if (walker.size() < nWalkers) {
            // Spawning new walkers until we have nWalkers.
            int spawns = nWalkers - walker.size();

            for (int i = 0; i < spawns; i++)
                walker.push_back(thermalizedWalker());
        }
 */

