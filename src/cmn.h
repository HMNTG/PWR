#include <stdio.h>

/**
 * @brief print command line args
*/
void print_args(int argc, char *argv[]);

/**
 * @brief every mpi thread prints out its rank and the section name.
 * @param name section name
*/
void dbg_sec_mpi(char* name);

/**
 * @brief print c style matrix
 * @param b pointer double matrix of size nxm
 * @param n rows
 * @param m columns
*/
void print_mat(double** a, int n, int m);

/**
 * @brief print c style array of double
 * @param b pointer double array of size n
 * @param n size of the array
*/
void print_vec(double* b, int n);
