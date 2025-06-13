#include "gauss_elimination.h"

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

// move to cmn
void print_mat(double** a, int n) {
    if (a == NULL || n <= 0) {
        printf("Invalid matrix or dimensions.\n");
        return;
    }

    // Calculate max_val_len once before the loop
    int max_val_len = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            int current_len = snprintf(NULL, 0, "%.1f", a[i][j]);
            if (current_len > max_val_len) {
                max_val_len = current_len;
            }
        }
    }

    // Calculate the effective content width for a row
    // This includes numbers and separators (", ")
    int total_row_content_width = 0;
    for (int j = 0; j < n; ++j) {
        total_row_content_width += max_val_len; // Width for the number itself
        if (j < n - 1) {
            total_row_content_width += 2; // Width for ", " separator
        }
    }

    for (int i = 0; i < n; ++i) {
        // --- CORRECTED PART 1: Print top border line before first data row ---
        if (i == 0) {
            printf("+");
            for (int k = 0; k < total_row_content_width; ++k) {
                printf("-");
            }
            printf("+\n"); // Newline after the top border
        }

        // --- Original logic for data rows (now always starts and ends with '|') ---
        printf("|");

        for (int j = 0; j < n; ++j) {
            printf("%*.*f", max_val_len, 1, a[i][j]);

            if (j < n - 1) { // Only print comma and space for non-last elements
                printf(", ");
            }
        }
        printf("|\n"); // Newline after each data row

        // --- CORRECTED PART 2: Print bottom border line after last data row ---
        if (i == n - 1) {
            printf("+");
            for (int k = 0; k < total_row_content_width; ++k) {
                printf("-");
            }
            printf("+\n"); // Newline after the bottom border
        }
    }
}

void print_vec(double* b, int n) {
    printf("(");
    for (int j=0;j<n;++j) {
        if (j==n-1)
            printf("%.1f,",b[j]);
        else
            printf("%.1f,",b[j]);
    }
    printf(")\n");
}

int max_col(double** a, int k, int n) {
    int cm=0; // current max
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
    double temp_b = b[r];

    for (int i=0;i<n;++i)
        a[r][i] = a[k][i];
    for (int i=0;i<n;++i)
        a[k][i] = temp_line[i];
}

double* gauss_sequential(double** a, double* b, int n) {
    if (a || b) {
        printf("we dont use a/b for now\n");
        return (double*) NULL;
    }
    n=10;
    a = malloc(n * sizeof(double*));
    b = malloc(n * sizeof(double));
    for (int i=0;i<n;++i) {
        a[i] = malloc(n * sizeof(double));
        for (int j=0;j<n;++j) {
            a[i][j] = (double) i*n + j;
        }
        b[i] = (double) i;
    }
    printf("a=\n");
    print_mat(a,n);
    printf("b=\n");
    print_vec(b,n);

    double* x, sum;
    //double l[MAX_SIZE];
    double* l = malloc(n * sizeof(double));

    int i,j,k,r;
    x = (double*) malloc(n * sizeof(double));
    for (k = 0; k < n-1; k++) {
        r = max_col(a,k,n);
        if (k != r) exchange_row(a,b,r,k,n);
        for (i=k+1; i < n; i++) {
            l[i] = a[i][k]/a[k][k];
            for (j=k; j < n; j++)
                a[i][j] = a[i][j] - l[i] * a[k][j];
            b[i] = b[i] - l[i] * b[k];
        }
        printf("a=\n");
        print_mat(a,n);
        printf("b=\n");
        print_vec(b,n);
    }
    for (k = n-1; k >= 0; k--) {
        sum = 0.0;
        for (j=k+1; j < n; j++)
            sum = sum + a[k][j] * x[j];
        x[k] = 1/a[k][k] * (b[k] - sum);
    }
    
    //return x;
    print_vec(x,n);
    free(x);

    //TODO do not free param
    for (int i=0;i<n;++i)
        free(a[i]);
    free(a);
    free(b);
    return NULL;
}
