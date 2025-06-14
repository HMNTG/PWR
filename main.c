#include "src/cmn.h"

// exercises
// ex 1
#include "src/gauss_elimination.h"
// exercises

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

int main(int argc, char *argv[]) {
    //int pid = atoi(argv[2]);
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int exercise=0; // note, all cores require to evaluate this
    print_args(argc,argv);
    if (argc>1)
        exercise = strtol(argv[1],NULL,10);
#if 0
    printf("Hello from process %d of %d\n", rank, size);
    fflush(stdout);
#endif
    switch (exercise)
    {
    case 0:
        dbg_sec(1); // print helper function
        gauss_elimination_par();
        dbg_sec(2);
        break;
    case 1:
        if (rank==0) {
            dbg_sec(1); // print helper function
            gauss_sequential(NULL,NULL,10);
            dbg_sec(2);
        }
        break;
    default:
        printf("id invalid\n");
        break;
    }

    MPI_Finalize();
    return 0;
}
