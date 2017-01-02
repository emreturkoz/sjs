#ifndef _SJSSOLVERHELPERS_H_
#define _SJSSOLVERHELPERS_H_


#include "sjsCommon.h"

class sjsSolverHelpers{

public:
	sjsSolverHelpers();
	~sjsSolverHelpers();

	// Fill a vector with zeroes
	void VectorZeroFill(int size, double* vec);

	// Copy the content of vec1 into vec2
	void TwoVectorCopy(int size, double* vec1, double* vec2);

	// Subtracts second vec from the first vec: first_vec - second_vec = result
	void TwoVectorSubtraction(int size, double* first_vec, double* second_vec, double* result);

	// Calculate the "two-norm" of a given vector
	double VectorTwoNorm(int size, double* a);

	// Matrix vector multiplication for banded matrices
	void BandedMatrixVectorMultiplication(double** A, double* x, double* b, int xdim, int problemDimension, int band);

	




private:

};





#endif // _SJSSOLVERHELPERS_H_