#include "gauss_elimination.h"

#include "cmn.h"

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

#include <math.h>

void fill_matrix(double A[N][N+1]) {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N+1; j++)
            A[i][j] = ((double)rand() / RAND_MAX) * 100.0 - 50.0; // Werte zwischen -50 und 50
}

int gauss_elimination_par() {
    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    static double A[N][N+1];
    static double x[N];

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

    return 0;
}

int max_col(double** a, int k, int n) {
    int cm=k; // current max
    for (int i=k;i<n;++i) {
        if (fabs(a[i][k]) > fabs(a[cm][k])) {
            cm = i;
        }
    }
    return cm;
}

void exchange_row(double** a, double* b, int r, int k, int n) {
    double* temp_line = malloc(n*sizeof(double));
    for (int i=0;i<n;++i)
        temp_line[i] = a[r][i];
    double temp_b = b[k];

    for (int i=0;i<n;++i)
        a[r][i] = a[k][i];
    for (int i=0;i<n;++i)
        a[k][i] = temp_line[i];
    b[k] = b[r];
    b[r] = temp_b;
}

double* gauss_sequential(double** a, double* b, int n) {
    if (a || b) {
        printf("we dont use a/b for now\n"); //TODO(skade)
        return (double*) NULL;
    }
    n=3;
    a = malloc(n * sizeof(double*));
    b = malloc(n * sizeof(double));

    for (int i=0;i<n;++i)
        a[i] = malloc(n*sizeof(double));
    // simple example
    a[0][0] =  2;  a[0][1] =  1;  a[0][2] = -1;  b[0] =  8;
    a[1][0] = -3;  a[1][1] = -1;  a[1][2] =  2;  b[1] = -11;
    a[2][0] = -2;  a[2][1] =  1;  a[2][2] =  2;  b[2] = -3;

    // auto fill
    //for (int i=0;i<n;++i) {
    //    a[i] = malloc(n * sizeof(double));
    //    for (int j=0;j<n;++j) {
    //        a[i][j] = (double) i*n + j;
    //    }
    //    b[i] = (double) i;
    //}
    printf("a=\n");
    print_mat_pp(n,n,a);
    printf("b=\n");
    print_vec(n,b);
    printf("\n");

    double* x, sum;
    //double l[MAX_SIZE];
    double* l = malloc(n * sizeof(double));

    int i,j,k,r;
    x = (double*) malloc(n * sizeof(double));
    for (k = 0; k < n-1; k++) {
        r = max_col(a,k,n);
        printf("iter %d:\n",k);
        if (k != r) exchange_row(a,b,r,k,n);
        printf("exchange row:\n");
        printf("a=\n");
        print_mat_pp(n,n,a);
        printf("b=\n");
        print_vec(n,b);
        for (i=k+1; i < n; i++) {
            l[i] = a[i][k]/a[k][k];
            for (j=k; j < n; j++)
                a[i][j] = a[i][j] - l[i] * a[k][j];
            b[i] = b[i] - l[i] * b[k];
        }
        printf("compute:\n");
        printf("a=\n");
        print_mat_pp(n,n,a);
        printf("b=\n");
        print_vec(n,b);
        printf("\n");
    }
    for (k = n-1; k >= 0; k--) {
        sum = 0.0;
        for (j=k+1; j < n; j++)
            sum = sum + a[k][j] * x[j];
        x[k] = 1/a[k][k] * (b[k] - sum);
    }
    
    //return x;
    print_vec(n,x);
    free(x);

    //TODO do not free param
    for (int i=0;i<n;++i)
        free(a[i]);
    free(a);
    free(b);
    return NULL;
}
