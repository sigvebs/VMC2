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
    double D = 0.5;

    //--------------------------------------------------------------------------
    // Reading data from QD.ini
    if (myRank == 0) {
        cout << "Reading configuration from file." << endl;
        ini INIreader("QD.ini");

        tau = INIreader.GetDouble("main", "dt");
        McSamples = (int) INIreader.GetDouble("DMC", "McSamples");
        DMCSamples = (int) INIreader.GetDouble("DMC", "DMCSamples");
        importanceSampling = INIreader.GetBool("DMC", "importanceSampling");
        thermalization = (int) INIreader.GetDouble("DMC", "thermalization");
        nWalkers = (int) INIreader.GetDouble("DMC", "nWalkers");
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
    int maxWalkers = nWalkers * 5;
    //--------------------------------------------------------------------------
    // Configuring the wave function.    
    //--------------------------------------------------------------------------
    if (myRank == 0) {
        cout << "Configuring the wave function." << endl;
    }

    idum = -1 - myRank - time(NULL);

    Orbital *orbital = new QDOrbital(dim, alpha, w);
    Jastrow *jastrow = NULL;
    if (usingJastrow)
        jastrow = new QDJastrow(dim, nParticles, beta);
    Hamiltonian *hamiltonian = new QDHamiltonian(dim, nParticles, w, usingJastrow);
    WaveFunction wf(dim, nParticles, -1 - myRank - time(NULL), orbital, jastrow, hamiltonian);

    // If we are using brute force sampling we have to calculate the optimal step length
    if (!importanceSampling)
        wf.setOptimalStepLength();
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
#if 1
    McSamples = correlationLength*nWalkers;
    double ETrial = 0;
    int cWalker = 0;

    for (int i = 0; i < McSamples + thermalization; i++) {
        for (int j = 0; j < nParticles; j++) {
            if (importanceSampling) {
                wf.tryNewPosition(j);
            } else {
                wf.tryNewPositionBF(j);
            }
        }

        if (i > thermalization) {
            ETrial += wf.sampleEnergy();

            // Cloning a walker
            if (i % (correlationLength) == 0) {
                walkers.push_back(wf.clone());
                cWalker++;

                if (myRank == 0) {
                    cout << "Cloning walker " << cWalker << "\xd";
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

#else //  Initiating each walker with its own VMC run.
    double ETrial = 0;
    int samples = 0;
    for (int k = 0; k < nWalkers; k++) {
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
        for (int i = 0; i < thermalization + McSamples; i++) {
            for (int j = 0; j < nParticles; j++) {
                if (importanceSampling)
                    wf->tryNewPosition(j);
                else
                    wf->tryNewPositionBF(j);
            }
            if (i >= thermalization) {
                ETrial += wf->sampleEnergy();
                samples++;
            }
        }

        // Storing walker.
        walkers.push_back(wf);
        cout << "walker " << k << " initialized" << "\xd";
    }

    ETrial /= samples;
    if (myRank == 0) {
        cout
                << walkers.size() << " walkers created" << endl
                << "Trial energy after VMC: E = " << ETrial << endl;
    }
#endif

    //--------------------------------------------------------------------------
    // Thermalizing walkers using DMC. The trial energy is updated each time 
    // the walkers are moved. The bad walkers are killed. The good ones get
    // cloned.
    //--------------------------------------------------------------------------
    if (myRank == 0) {
        cout << "Starting DMC thermalization" << endl;
        writeDistributionToFile("DATA/DMC/distributionVMC.dat");
    }
    //exit(0);

    WaveFunction *toBeKilled;
    int walkersLength, NSamples, copies;
    double EOld, ENew, ETmp, EAvg, bFactor;
    int nKilled, nCloned;
    EAvg = 0;
    for (int n = 1; n <= DMCThermal; n++) {
        ETmp = 0;
        NSamples = 0;
        nKilled = 0;
        nCloned = 0;
        // Looping over every walker in walkers once.
        walkersLength = walkers.size();

        for (int i = 0; i < walkersLength; i++) {
            EOld = walkers[i]->sampleEnergy();

            // Moving walker.
            for (int j = 0; j < nParticles; j++) {
                if (importanceSampling) {
                    walkers[i]->DMCtryNewPosition(j);
                } else {
                    walkers[i]->tryNewPositionBF(j);
                }
            }
            ENew = walkers[i]->sampleEnergy();

            // Calculating the branching factor.
            bFactor = exp(-tau * (0.5 * (ENew + EOld) - ETrial));

            //--------------------------------------------------------------
            copies = int(bFactor + ran3(&idum));

            // The walker is killed.
            if (copies == 0 || ENew < ETrial - 1 / sqrt(tau) || ENew > ETrial + 1 / sqrt(tau)) {
                toBeKilled = walkers[i];
                walkers.erase(walkers.begin() + i);
                delete toBeKilled;
                walkersLength--;
                i--;
                nKilled++;
            } else {
                // Cloning walker.
                if (walkers.size() < maxWalkers) {
                    for (int j = 1; j < copies; j++) {
                        walkers.push_back(walkers[i]->clone());
                        nCloned++;
                    }
                }
                // Sampling energy.
                ETmp += bFactor*ENew;
                NSamples++;
            }
            //------------------------------------------------------------------
        }
        //----------------------------------------------------------------------
        // Updating the new trial energy and printing progress.
        ETmp /= NSamples;
        EAvg += ETmp;
        ETrial = EAvg / n;
        //ETrial = EAvg / n - log(walkers.size() / double(nWalkers));

        if (myRank == 0) {
            cout
                    << "TERMALIZATION: "
                    << "\tn = " << n
                    << "\tnWalkers = " << walkers.size()
                    << "\t ETmp = " << ETmp
                    << "\t EAvg = " << EAvg / n
                    << "\t nKilled = " << nKilled
                    << "\t nCloned = " << nCloned
                    << "\xd";
            //<< endl;
        }
    }

    //--------------------------------------------------------------------------
    // DMC sampling. The trial energy is updated after every walker has moved
    // a block. 
    //--------------------------------------------------------------------------    
    ofstream outStream;
    if (myRank == 0) {
        cout
                << "Starting DMC sampling. "
                << "\twith Et = " << ETrial
                << endl;
        outStream.open((const char*) &fileName[0]);
        // Writing to file
        outStream << nParticles << " " << w << " " << DMCSamples << " " << blockSize << endl;

        writeDistributionToFile("DATA/DMC/distributionThermal.dat");
    }

    EAvg = 0;
    double EBlock;

    for (int n = 1; n <= DMCSamples; n++) {

        // Looping over every walker in walkers once.
        walkersLength = walkers.size();
        NSamples = 0;
        EBlock = 0;
        for (int i = 0; i < walkersLength; i++) {

            // Moving one walker a block
            for (int b = 0; b < blockSize; b++) {
                EOld = walkers[i]->sampleEnergy();

                // Moving walker.
                for (int j = 0; j < nParticles; j++) {
                    if (importanceSampling) {
                        walkers[i]->DMCtryNewPosition(j);
                    } else {
                        walkers[i]->tryNewPositionBF(j);
                    }
                }
                ENew = walkers[i]->sampleEnergy();

                // Calculating the branching factor
                bFactor = exp(-tau * (0.5 * (ENew + EOld) - ETrial));

                //--------------------------------------------------------------
                copies = int(bFactor + ran3(&idum));

                // The walker is killed.
                if (copies == 0 || ENew < ETrial - 1 / sqrt(tau) || ENew > ETrial + 1 / sqrt(tau)) {
                    toBeKilled = walkers[i];
                    walkers.erase(walkers.begin() + i);
                    delete toBeKilled;
                    walkersLength--;
                    i--;
                    b = blockSize;
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
                    << "\xd";

            // Writing to file
            if (myRank == 0) {
                outStream << EBlock << endl;
            }
        }
    }

    //--------------------------------------------------------------------------
    if (myRank == 0)
        outStream.close();

    writeDistributionToFile("DATA/DMC/distributionDMC.dat");
}

////////////////////////////////////////////////////////////////////////////////

DMC::DMC(const DMC & orig) {
}

////////////////////////////////////////////////////////////////////////////////

DMC::~DMC() {
}

void DMC::writeDistributionToFile(string fName) {
    // Writing the current distribution of particles to file
    mat r;
    ofstream out;
    out.open((const char*) &fName[0]);
    for (int i = 0; i < walkers.size(); i++) {
        r = walkers[i]->getROld();
        for (int j = 0; j < nParticles; j++) {
            for (int d = 0; d < dim; d++) {
                out << r(j, d) << " ";
            }
            out << endl;
        }
    }
    out.close();
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

