#include "matrix.h"
#include <cstdlib>

float **mkmatrixf(int row, int col){
	int i;
	float **matrix;
	matrix = (float **) malloc ( row * sizeof(float*) );
	for(i=0; i<row; i++)    matrix[i] = (float *) malloc ( col * sizeof(float) );
	return matrix;

}

void printmatrixf(char const *comment, float **A, int row, int col){
	printf("<:: %s ::> (%d,%d)\n",comment,row,col);
	for (int i = 0; i < row; i++){
		for (int j = 0; j < col; j++){
			printf("%f  ",A[i][j]);
		}
		printf("\n");
	}
}

double **mkmatrixd(int row, int col){
	int i;
	double **matrix;
	matrix = (double **) malloc ( row * sizeof(double*) );
	for(i=0; i<row; i++)    matrix[i] = (double *) malloc ( col * sizeof(double) );
	return matrix;

}

void printmatrixd(char const *comment, double **A, int row, int col){
	printf("<:: %s ::> (%d,%d)\n",comment,row,col);
	for (int i = 0; i < row; i++){
		for (int j = 0; j < col; j++){
			printf("%f  ",A[i][j]);
		}
		printf("\n");
	}
}

float *mkvectorf(int num){
	float *vector;
	vector = (float *) malloc ( num * sizeof(float) );
	return vector;
}

double *mkvectord(int num){
	double *vector;
	vector = (double *) malloc ( num * sizeof(double) );
	return vector;
}

char *mkvectorc(int NN){
        char *vector;
        vector = (char *) malloc ( NN * sizeof(char) );
        //if( vector == NULL )    my_alloc_error("mkvectorc");
        return vector;
}
