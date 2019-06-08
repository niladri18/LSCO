#include <iostream>
#include <stdio.h>
#include <complex>
#include <cmath>
#include <math.h>
using namespace std;


float **mkmatrixf(int row, int col);
void printmatrixf(char const *comment, float **A, int row, int col);
double **mkmatrixd(int row, int col);
void printmatrixd(char const *comment, double **A, int row, int col);
float *mkvectorf(int num);
double *mkvectord(int num);
char *mkvectorc(int NN);
