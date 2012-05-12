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
#include "ComputeGrid.h"
#include "OneBodyDensity.h"
#include "SGD.h"
#include "Blocking.h"

using namespace std;

int main(int argc, char** argv) {
    int option;

    // Collecting MPI information.
    int myRank, nNodes;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &nNodes);

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
        case 1:
            new ComputeGrid();
            break;
        case 2:
            new OneBodyDensity();
            break;
        case 3:
            new SGD();
            break;
        case 4:
            new ComputeOneWF();
            if (myRank == 0)
                new Blocking(nNodes);
            break;
    }

    MPI_Finalize();
    return 0;
}

