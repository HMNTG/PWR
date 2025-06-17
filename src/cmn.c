#include "cmn.h"

#include <mpi.h>

void print_args(int argc, char *argv[]) {
    printf("argc: %d\n",argc);
    for (int i=0;i<argc;++i) {
        printf("%s\n",argv[i]);
        printf("\n");
    }
    fflush(stdout);
}

void dbg_sec_mpi(char* name) {
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("sec: %s proc: %d\n",name,rank);fflush(stdout);
}

void print_mat_p(int n, int m, double* a) {
    double** b = malloc(n*sizeof(double));
    for (int i=0;i<n;++i)
        b[i] = malloc(n*sizeof(double));

    for (int i=0;i<n;++i)
        for (int j=0;j<m;++j) {
            b[i][j] = a[i*m + j];
        }
    print_mat_pp(n,m,b);
}

void print_mat_pp(int n, int m, double** a) {
    if (a == NULL || n <= 0 || m <= 0) {
        printf("Invalid matrix or dimensions.\n");
        return;
    }

    // Calculate max_val_len once before the loop
    int max_val_len = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            int current_len = snprintf(NULL, 0, "%.1f", a[i][j]);
            if (current_len > max_val_len) {
                max_val_len = current_len;
            }
        }
    }

    // Calculate the effective content width for a row
    int total_row_content_width = 0;
    for (int j = 0; j < m; ++j) {
        total_row_content_width += max_val_len; // Width for the number itself
        if (j < m - 1) {
            total_row_content_width += 2; // Width for ", " separator
        }
    }

    for (int i = 0; i < n; ++i) {
        // Print top border before the first row
        if (i == 0) {
            printf("+");
            for (int k = 0; k < total_row_content_width; ++k) {
                printf("-");
            }
            printf("+\n");
        }

        printf("|");
        for (int j = 0; j < m; ++j) {
            printf("%*.*f", max_val_len, 1, a[i][j]);
            if (j < m - 1) {
                printf(", ");
            }
        }
        printf("|\n");

        // Print bottom border after the last row
        if (i == n - 1) {
            printf("+");
            for (int k = 0; k < total_row_content_width; ++k) {
                printf("-");
            }
            printf("+\n");
        }
    }
}

void print_vec(int n, double* b) {
    printf("(");
    for (int j=0;j<n;++j) {
        if (j==n-1)
            printf("%.1f",b[j]);
        else
            printf("%.1f,",b[j]);
    }
    printf(")\n");
}