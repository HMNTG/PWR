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

void print_mat(double** a, int n, int m) { //TODO(skade) impl m
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
            printf("%.1f",b[j]);
        else
            printf("%.1f,",b[j]);
    }
    printf(")\n");
}