#include "cmn.h"

#include <mpi.h>

void dbg_sec(int sec) {
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("sec: %d proc: %d\n",sec,rank);fflush(stdout);
}

void print_args(int argc, char *argv[]) {
    printf("argc: %d\n",argc);
    for (int i=0;i<argc;++i) {
        printf("%s\n",argv[i]);
        printf("\n");
    }
    fflush(stdout);
}
