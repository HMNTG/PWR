
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "cmn.h"

#define N_MIN   10
#define N_MAX   10 // matrix size
#define N_MULT  2

#define RANDOM          0
#define DIAGDOM_SCALE   0.5
#define MAX_ITERATIONS  100

#define TOLERANCE     0.01

#define ROOT                0

#define PRINT_GLOBALS       1
#define PRINT_LOCALS        1
#define PRINT_ITERATIONS    1

#define JACOBI 0
#define GAUSS_SEIDEL 1

void print_matrix(int n, int m, double *a) {
    int i, j;

    for (i = 0; i < n; ++i) {
        for (j = 0; j < m; ++j) printf("  %f", a[i * m + j]);
        printf("\n");
    }
}


void random_matrix(int n, int m, double *a, double max_val) {
    int i, j;

    for (i = 0; i < n; ++i) {
        for (j = 0; j < m; ++j) {
#if RANDOM
            a[i * m + j] = random() / (double) RAND_MAX * max_val;
#else
            a[i * m + j] = 1 * max_val;
#endif
        }
    }
}


void random_matrix_diagdom(int n, double *a, double max_diag) {
    const double scale = DIAGDOM_SCALE * max_diag / (n - 1);

    int i, j;

    for (i = 0; i < n; ++i) {
        double x = 0;

        for (j = 0; j < n; ++j) {
#if RANDOM
            a[i * n + j] = random() / (double) RAND_MAX * scale;
#else
            a[i * n + j] = 1 * scale;
#endif

            x += a[i * n + j];
        }

        x -= a[i * n + i];

        a[i * n + i] = x / DIAGDOM_SCALE;
    }
}


void
distribute_rows(int n, int m, double *a, int *local_n, double **local_a, MPI_Comm comm, int comm_size, int comm_rank) {
    *local_n = n / comm_size;

    *local_a = malloc(*local_n * m * sizeof(double));

    MPI_Scatter(a, *local_n * m, MPI_DOUBLE, *local_a, *local_n * m, MPI_DOUBLE, ROOT, comm);
}


void
distribute_cols(int n, int m, double *a, int *local_m, double **local_a, MPI_Comm comm, int comm_size, int comm_rank) {
    *local_m = m / comm_size;

    *local_a = malloc(n * *local_m * sizeof(double));

    MPI_Datatype scols_type = MPI_DATATYPE_NULL;
    MPI_Datatype scols_type_resized = MPI_DATATYPE_NULL;

    if (comm_rank == ROOT) {
        MPI_Type_vector(n, *local_m, m, MPI_DOUBLE, &scols_type);
        MPI_Type_commit(&scols_type);

        MPI_Type_create_resized(scols_type, 0, *local_m * sizeof(double), &scols_type_resized);
        MPI_Type_commit(&scols_type_resized);
    }

    MPI_Scatter(a, 1, scols_type_resized, *local_a, n * *local_m, MPI_DOUBLE, ROOT, comm);

    if (comm_rank == ROOT) {
        MPI_Type_free(&scols_type);
        MPI_Type_free(&scols_type_resized);
    }
}


double distance(double *x0, double *x1, int n, MPI_Comm comm) {
    double l[2] = {0, 0}, g[2];

    // TODO: lokale Abstandsberechnung
    for (int i = 0; i < n; ++i) {
        // S247
        l[0] = (x0[0] - x1[0]); //TODO
        l[1] = 0.;
    }

    MPI_Allreduce(l, g, 2, MPI_DOUBLE, MPI_SUM, comm);

    return sqrt(g[0]) / sqrt(g[1]);
}


int parallel_jacobi(int n, int local_n, double *local_a, double *local_b, double *local_x, double tol, MPI_Comm comm,
                    int comm_size, int comm_rank) {
    double *t, *x_tmp = malloc(2 * n * sizeof(double));
    double *x_old = x_tmp + 0 * n;
    double *x_new = x_old + 1 * n;

    MPI_Allgather(local_b, local_n, MPI_DOUBLE, x_new, local_n, MPI_DOUBLE, comm);

    int it = 0, cont = 0;

    // iterate until error is smaller then tolerance or max iteration count is reached
    do {

        ++it;

        // switch pointer to always work with the latest x vector
        t = x_new;
        x_new = x_old;
        x_old = t;

        // TODO: Schleife über i_local bis local_n
        int i_local;
        for (i_local = 0;i_local<local_n;++i_local) {
            // TODO: Berechne aktuellen globalen Index i_global
            int i_global = i_local + comm_rank*local_n;

            // TODO: nehme aktuellen eintrag von local_b als Startwert für local_x
            // nutze local_x
            local_x[i_local] = local_b[i_local];

            // TODO: Subtrahiere die Werte der zu i_local gehörigen Zeile multipliziert mit dem entsprechendem x_old
            //       Das i_global Element darf NICHT betrachtet werden!
            for (int j=0;j<n;++j) {
                if (i_local==i_global)
                    local_x[i_local] = local_x[i_local] - 0.;
                else
                    local_x[i_local] = local_x[i_local] - local_a[i_local*n + j]*(*x_old);
            }

            // TODO: Dividiere das local_x mit dem aktuellem Wert der Matrix A an Stelle von i_global
            local_x[i_local] /= local_x[i_local*n + i_global];
        }

        MPI_Allgather(local_x, local_n, MPI_DOUBLE, x_new, local_n, MPI_DOUBLE, comm);

        double err = distance(x_old + comm_rank * local_n, local_x, local_n, comm);

#if PRINT_ITERATIONS
# if PRINT_LOCALS
        printf("%d: local_x:\n", comm_rank);
        //print_matrix(local_n, 1, local_x);
        print_mat_p(local_n, 1, local_x);
        printf("\n");
# endif

        if (comm_rank == ROOT) printf("%d: iteration: %d, error: %f\n", comm_rank, it, err);
#endif

        cont = (err > tol);

    } while (it < MAX_ITERATIONS && cont);

    free(x_tmp);

    return it;
}


int parallel_gauss_seidel(int n, int local_m, double *local_a, double *local_b, double *local_x, double tol, MPI_Comm comm,
                          int comm_size, int comm_rank) {
    double *local_x_old = malloc(local_m * sizeof(double));

    int j_local;
    for (j_local = 0; j_local < local_m; ++j_local) local_x[j_local] = 0;

    int it = 0, cont = 0;

    // iterate until error is smaller then tolerance or max iteration count is reached
    do {

        ++it;

        for (j_local = 0; j_local < local_m; ++j_local) local_x_old[j_local] = local_x[j_local];

        int i_global;
        for (i_global = 0; i_global < n; i_global++) {
            double s_k = 0.0;

            for (j_local = 0; j_local < local_m; ++j_local)
                if (j_local + comm_rank * local_m != i_global)
                    s_k += local_a[i_global * local_m + j_local] * local_x[j_local];

            int root = i_global / local_m;

            int i_local = i_global % local_m;

            MPI_Reduce(&s_k, &local_x[i_local], 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);

            if (comm_rank == root)
                local_x[i_local] = (local_b[i_local] - local_x[i_local]) / local_a[i_global * local_m + i_local];
        }

        double err = distance(local_x_old, local_x, local_m, comm);

        cont = (err > tol);

#if PRINT_ITERATIONS
# if PRINT_LOCALS
        printf("%d: local_x:\n", comm_rank);
        //print_matrix(local_m, 1, local_x);
        print_mat_p(local_m, 1, local_x);
        printf("\n");
# endif

        if (comm_rank == ROOT) printf("%d: iteration: %d, error: %f\n", comm_rank, it, err);
#endif

    } while (it < MAX_ITERATIONS && cont);

    free(local_x_old);

    return it;
}


int PWR_ex9_main(int argc, char *argv[]) {
    //MPI_Init(&argc, &argv);

    MPI_Comm comm = MPI_COMM_WORLD;
    int comm_size, comm_rank;
    MPI_Comm_size(comm, &comm_size);
    MPI_Comm_rank(comm, &comm_rank);

    double *a = NULL, *b = NULL, *x = NULL;

    if (comm_rank == ROOT) {
        a = malloc(N_MAX * N_MAX * sizeof(double));
        b = malloc(N_MAX * sizeof(double));
        x = malloc(N_MAX * sizeof(double));
    }

    int nn;
    for (nn = N_MIN; nn <= N_MAX; nn *= N_MULT) {
        int n = nn / comm_size;
        n *= comm_size;

        if (n < comm_size) n = comm_size;

        if (comm_rank == ROOT) {
            random_matrix_diagdom(n, a, 2);
            random_matrix(n, 1, b, 1);

#if PRINT_GLOBALS
            printf("%d: A:\n", comm_rank);
            //print_matrix(n, n, a);
            print_mat_p(n, n, a);
            printf("\n");

            printf("%d: b:\n", comm_rank);
            //print_matrix(n, 1, b);
            print_mat_p(n, 1, b);
            printf("\n");
#endif
        }

        int local_n;
        double *local_b;

        distribute_rows(n, 1, b, &local_n, &local_b, comm, comm_size, comm_rank);

#if PRINT_LOCALS
        printf("%d: local_b:\n", comm_rank);
        //print_matrix(local_n, 1, local_b);
        print_mat_p(local_n, 1, local_b);
        printf("\n");
#endif

        double *local_x = malloc(local_n * sizeof(double));

#if JACOBI
#if PRINT_GLOBALS
        if (comm_rank == ROOT) printf("%d: Jacobi\n", comm_rank);
#endif

        double *local_a_rows;

        distribute_rows(n, n, a, &local_n, &local_a_rows, comm, comm_size, comm_rank);

#if PRINT_LOCALS
        printf("%d: local_a:\n", comm_rank);
        //print_matrix(local_n, n, local_a_rows);
        print_mat_p(local_n, n, local_a_rows);
        printf("\n");
#endif

        MPI_Barrier(comm);
        double t_jac = MPI_Wtime();

        int it_jac = parallel_jacobi(n, local_n, local_a_rows, local_b, local_x, TOLERANCE, comm, comm_size, comm_rank);

        MPI_Barrier(comm);
        t_jac = MPI_Wtime() - t_jac;

#if PRINT_LOCALS
        printf("%d: local_x:\n", comm_rank);
        //print_matrix(local_n, 1, local_x);
        print_mat_p(local_n, 1, local_x);
        printf("\n");
#endif

        MPI_Gather(local_x, local_n, MPI_DOUBLE, x, local_n, MPI_DOUBLE, ROOT, comm);

#if PRINT_GLOBALS
        if (comm_rank == ROOT) {
            printf("%d: x:\n", comm_rank);
            //print_matrix(n, 1, x);
            print_mat_p(n, 1, x);
            printf("\n");
        }
#endif

        free(local_a_rows);
#endif

#if GAUSS_SEIDEL
#if PRINT_GLOBALS
        if (comm_rank == ROOT) printf("%d: Gauss-Seidel\n", comm_rank);
#endif

        int local_m;
        double *local_a_cols;

        distribute_cols(n, n, a, &local_m, &local_a_cols, comm, comm_size, comm_rank);

#if PRINT_LOCALS
        printf("%d: local_a:\n", comm_rank);
        //print_matrix(n, local_m, local_a_cols);
        print_mat_p(n, local_m, local_a_cols);
        printf("\n");
#endif

        MPI_Barrier(comm);
        double t_gs = MPI_Wtime();

        int it_gs = parallel_gauss_seidel(n, local_m, local_a_cols, local_b, local_x, TOLERANCE, comm, comm_size,
                                          comm_rank);

        MPI_Barrier(comm);
        t_gs = MPI_Wtime() - t_gs;

#if PRINT_LOCALS
        printf("%d: local_x:\n", comm_rank);
        //print_matrix(local_n, 1, local_x);
        print_mat_p(local_n, 1, local_x);
        printf("\n");
#endif

        MPI_Gather(local_x, local_n, MPI_DOUBLE, x, local_n, MPI_DOUBLE, ROOT, comm);

#if PRINT_GLOBALS
        if (comm_rank == ROOT) {
            printf("%d: x:\n", comm_rank);
            //print_matrix(n, 1, x);
            print_mat_p(n, 1, x);
            printf("\n");
        }
#endif

        free(local_a_cols);
#endif

        //if (comm_rank == ROOT)
        //    printf("%5d  %f  %f  %3d  %f  %f  %3d\n", n, t_jac, t_jac / it_jac, it_jac, t_gs, t_gs / it_gs, it_gs);
    }

    if (comm_rank == ROOT) {
        free(a);
        free(b);
        free(x);
    }

    MPI_Finalize();

    return 0;
}
