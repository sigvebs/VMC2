/* 
 * File:   Blocking.cpp
 * Author: sigve
 * 
 * Created on 11. mai 2012, 16:06
 */

#include "Blocking.h"

#include "includes/lib.h"
#include <fstream>
#include <string.h>

//#include <armadillo>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
////////////////////////////////////////////////////////////////////////////////

Blocking::Blocking() {

}

////////////////////////////////////////////////////////////////////////////////

Blocking::Blocking(int nNodes) {
    int localN;
    int N;

    cout << "Starting blocking analysis." << endl;

    //--------------------------------------------------------------------------
    // Collecting file size information.
    struct stat result;
    if (stat("DATA/blocking_0.dat", &result) == 0) {
        localN = result.st_size / sizeof (double);
        N = localN*nNodes;
    } else {
        cout << "Could not find file.";
        exit(1);
    }
    
    int deltaBlockSize = 10;
    int maxBlockSize = (int) N/20;
    
    if(maxBlockSize > 1e4){
        maxBlockSize = 1e4;
    }
        
    //--------------------------------------------------------------------------
    // Collecting data
    cout << "Reading data from files." << endl;

    vec McData = zeros(N, 1);
    double* McResults = new double[localN];
    int k = 0;
    
    for (int i = 0; i < nNodes; i++) {
        ostringstream ost;
        ost << "DATA/blocking_" << i << ".dat";
        ifstream infile;
        infile.open(ost.str().c_str(), ios::in | ios::binary);
        infile.read((char*) McResults, result.st_size);
        infile.close();

        for (int j = 0; j < localN; j++)
            McData[k++] = McResults[j];
        }

    delete McResults;
    //--------------------------------------------------------------------------
    // Looping over block sizes and storing results.
    cout << "Looping over block sizes and storing results." << endl;

    int blockSize;
    vec results;
    double mean, sigma;

    ofstream outStream;
    outStream.open("DATA/blockingResults.dat");

    for (int i = 1; i*deltaBlockSize <= maxBlockSize; i++) {
        blockSize = i*deltaBlockSize;
        results = block(McData, blockSize);
        mean = results(0);
        sigma = results(1);
        outStream << blockSize << "\t" << mean << "\t" << sigma << "\n";
    }
    outStream.close();
}

////////////////////////////////////////////////////////////////////////////////

vec Blocking::block(vec McData, int blockSize) {
    int NSamples = McData.size();
    int NBlocks = int(NSamples / blockSize);
    vec results = zeros(2, 1);
    vec EBlock = zeros(NBlocks, 1);

    // Finding the mean values of the blocks
    double deltaE;
    for (int i = 0; i < NBlocks; i++) {
        // Finding the average over one block
        deltaE = 0;
        for (int j = i * blockSize; j < i * blockSize + blockSize; j++)
            deltaE += McData(j);

        EBlock(i) = deltaE / blockSize;
    }

    // Calculating the mean and variance of all the blocks
    double E = 0;
    double E2 = 0;
    for (int i = 0; i < NBlocks; i++) {
        E += EBlock(i);
        E2 += EBlock(i) * EBlock(i);
    }
    E /= NBlocks;
    E2 /= NBlocks;

    double sigma = E2 - E*E;
    sigma = sqrt(sigma/NBlocks);
    results(0) = E;
    results(1) = sigma;

    return results;
}

////////////////////////////////////////////////////////////////////////////////

Blocking::Blocking(const Blocking& orig) {

}

////////////////////////////////////////////////////////////////////////////////

Blocking::~Blocking() {

}
////////////////////////////////////////////////////////////////////////////////