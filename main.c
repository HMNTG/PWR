#include "src/cmn.h"

// exercises
// ex 1
#include "src/gauss_elimination.h"
// ex 2
#include "PWR25_Uebung09_Template.h"
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

    int exercise=2; // note, all cores require to evaluate this
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
        dbg_sec_mpi("one"); // print helper function
        gauss_elimination_par();
        dbg_sec_mpi("two");
        break;
    case 1:
        if (rank==0) {
            dbg_sec_mpi("one"); // print helper function
            int n=3;
            double** a = malloc(n * sizeof(double*));
            double* b = malloc(n * sizeof(double));

            for (int i=0;i<n;++i)
                a[i] = malloc(n*sizeof(double));
            // simple example
            a[0][0] =  2;  a[0][1] =  1;  a[0][2] = -1;  b[0] =  8;
            a[1][0] = -3;  a[1][1] = -1;  a[1][2] =  2;  b[1] = -11;
            a[2][0] = -2;  a[2][1] =  1;  a[2][2] =  2;  b[2] = -3;

            //for (int i=0;i<n;++i) { //auto fill
            //    a[i] = malloc(n * sizeof(double));
            //    for (int j=0;j<n;++j)
            //        a[i][j] = (double) i*n + j;
            //    b[i] = (double) i;
            //}
            gauss_sequential(a,b,n);

            for (int i=0;i<n;++i)
                free(a[i]);
            free(a);
            free(b);
            dbg_sec_mpi("two");
        }
        break;
    case 2:
        PWR_ex9_main(argc,argv);
    default:
        printf("id invalid\n");
        break;
    }

    MPI_Finalize();
    return 0;
}
