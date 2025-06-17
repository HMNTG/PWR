#include <stdio.h>
#include <stdlib.h>

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
 * @param n rows
 * @param m columns
 * @param a pointer double matrix of size nxm
*/
void print_mat_p(int n, int m, double* a);

/**
 * @brief print c style matrix
 * @param n rows
 * @param m columns
 * @param a pointer of pointers double matrix of size nxm
*/
void print_mat_pp(int n, int m,double** a);

/**
 * @brief print c style array of double
 * @param n size of the array
 * @param b pointer double array of size n
*/
void print_vec(int n, double* b);
