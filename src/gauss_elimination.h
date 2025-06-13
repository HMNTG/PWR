#define N 1000 // Größe des Gleichungssystems

void fill_matrix(double A[N][N+1]);
int gauss_elimination_par();


int max_col(double** a, int k, int n);
void exchange_row(double** a, double* b, int r, int k, int n);

/**
 * @param a nxn matrix
 * @param b n vector
*/
double* gauss_sequential (double** a, double* b, int n);
