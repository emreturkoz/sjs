#include "sjsSolverHelpers.h"

/// Constructor
sjsSolverHelpers::sjsSolverHelpers(){

}

/// Destructor
sjsSolverHelpers::~sjsSolverHelpers(){

}

/// Fill a vector with zeroes
void sjsSolverHelpers::VectorZeroFill(int size, double* vec){
	for(int i = 0; i<size; i++){
		vec[i] = 0;
	}

}

/// Copy the content of vec1 into vec2
void sjsSolverHelpers::TwoVectorCopy(int size, double* vec1, double* vec2){
	for(int i = 0; i<size; i++){
		vec2[i] = vec1[i];
	}

}

/// Subtracts second vec from the first vec: first_vec - second_vec = result
void sjsSolverHelpers::TwoVectorSubtraction(int size, double* first_vec, double* second_vec, double* result){
	for(int i = 0; i<size; i++){
		result[i] = first_vec[i] - second_vec[i];
	}

}


/// Calculate the "two-norm" of a given vector
double sjsSolverHelpers::VectorTwoNorm(int size, double* a){
	double norm = 0;
	for(int i = 0; i<size; i++){
		norm = norm + a[i]*a[i];
	}

	return sqrt((double)norm);
}

// Matrix vector multiplication for banded matrices
void sjsSolverHelpers::BandedMatrixVectorMultiplication(double** A, double* x, double* b, int xdim, int problemDimension, int band){
	int size;
	size = problemDimension;
	VectorZeroFill(size, b);
	for(int i =0; i<size; i++){
		if(band==5){
			//**********First row of the matrix**********//
			if(i==0){
				b[i] = b[i] + x[i]*A[2][i] + x[i+1]*A[3][i] + x[i+xdim]*A[4][i];
			}
			//**********Last row of the matrix**********//
			if(i==size-1){
				b[i] = b[i] + x[i]*A[2][i] + x[i-1]*A[1][i] + x[i-xdim]*A[0][i];	

			}
			//**********Second part of the matrix**********//
			if(i<xdim && i>0){
				b[i] = b[i] + x[i]*A[2][i] + x[i+1]*A[3][i] + x[i+xdim]*A[4][i] + x[i-1]*A[1][i];
			}
			//**********Fourth part of the matrix**********//
			if(i>size-xdim && i<size-1){
				b[i] = b[i] + x[i]*A[2][i] + x[i+1]*A[3][i] + x[i-xdim]*A[0][i] + x[i-1]*A[1][i];	
			}
			//**********Interior part of the matrix**********//
			if(i>=xdim && i<size-xdim){
				b[i] = b[i] + x[i]*A[2][i] + x[i+1]*A[3][i] + x[i-xdim]*A[0][i] + x[i-1]*A[1][i] + x[i+xdim]*A[4][i];	
			}


		}
		else if (band == 3){
			//**********First row of the matrix**********//
			if(i==0){
				b[i] = b[i] + x[i]*A[1][i] + x[i+1]*A[2][i];
			}
			//**********Last row of the matrix**********//
			if(i==size-1){
				b[i] = b[i] + x[i]*A[1][i] + x[i-1]*A[0][i];	

			}
			//**********Interior part of the matrix**********//
			if(i>0 && i<size-1){
				b[i] = b[i] + x[i]*A[1][i] + x[i-1]*A[0][i] + x[i+1]*A[2][i];	
			}						

		}

	}

}

	