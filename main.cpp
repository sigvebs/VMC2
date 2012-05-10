/* 
 * File:   main.cpp
 * Author: Sigve
 *
 * Created on May 9, 2012, 10:22 PM
 */

#include <cstdlib>
#include <mpi.h>
#include "includes/ini.h"
#include "ComputeOneWF.h"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
    int option;
    
    // Collecting MPI information.
    int myRank;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    // Reading data from QD.ini
    if (myRank == 0) {
        ini INIreader("QD.ini");
        option = (int) INIreader.GetInt("main", "option");
        cout << "Starting option: " << option << endl;
    }

    // Broadcasting option to all nodes.
    MPI_Bcast(&option, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Running option.
    switch (option) {
        case 0:
            new ComputeOneWF();
            break;
    }

    MPI_Finalize();
    return 0;
}

