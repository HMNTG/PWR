#define N 1000 // Größe des Gleichungssystems

void fill_matrix(double A[N][N+1]);
int gauss_elimination_par();

/**
 * @brief returns row idx of max col for gauss_sequential
 * @param a matrix (nxn)
 * @param k gauss iteration step
 * @param n matrix dimension
*/
int max_col(double** a, int k, int n);

/**
 * @brief swaps rows r and k in matrix a and vector b
 * @param a matrix (nxn)
 * @param b result vector
 * @param r first idx to swap
 * @param k second idx to swap
 * @param n matrix / vector dimension
*/
void exchange_row(double** a, double* b, int r, int k, int n);

/**
 * @param a nxn matrix
 * @param b n vector
*/
double* gauss_sequential(double** a, double* b, int n);
