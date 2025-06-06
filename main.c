#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

#define N 1000 // Größe des Gleichungssystems

void fill_matrix(double A[N][N+1]) {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N+1; j++)
            A[i][j] = ((double)rand() / RAND_MAX) * 100.0 - 50.0; // Werte zwischen -50 und 50
}

int main(int argc, char *argv[]) {
    int rank, size;
    static double A[N][N+1];
    static double x[N];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        srand((unsigned int)time(NULL));
        fill_matrix(A);
    }
    MPI_Bcast(A, N*(N+1), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Gauß-Elimination
    for (int k = 0; k < N; k++) {
        if (rank == 0) {
            double pivot = A[k][k];
            for (int j = k; j < N+1; j++)
                A[k][j] /= pivot;
        }
        MPI_Bcast(&A[k], N+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        for (int i = k+1+rank; i < N; i += size) {
            double factor = A[i][k];
            for (int j = k; j < N+1; j++)
                A[i][j] -= factor * A[k][j];
        }
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, A, (N+1)*N/size, MPI_DOUBLE, MPI_COMM_WORLD);
    }

    // Rückwärtseinsetzen (nur von Prozess 0)
    if (rank == 0) {
        for (int i = N-1; i >= 0; i--) {
            x[i] = A[i][N];
            for (int j = i+1; j < N; j++)
                x[i] -= A[i][j] * x[j];
        }
        printf("Fertig! (Loesung nicht ausgegeben, da N=1000)\n");
    }

    MPI_Finalize();
    return 0;
}