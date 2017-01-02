#ifndef _SJSSOLVERS_H_
#define _SJSSOLVERS_H_


#include "sjsCommon.h"
#include "sjsSolverHelpers.h"


class sjsSolvers{
public:

	// Constructor
	sjsSolvers();

	// Destructor
	~sjsSolvers();

	// Set the coefficient matrix
	void SetCoefficientMatrix(double **A, bool printToFile);

	// Set the RHS
	void SetRHS(double* b, bool printToFile);

	// Set number of grid points
	void SetProblemDimension(int size);

	// Set the solution vector for the current iteration
	void SetSolutionVector(double* x);

	// Get the solution vector for the current iteration
	double* GetSolutionVector();

	// Set number of points in axial direction
	void SetNumNodesAxial(int nodes);

	// Set coefficient matrix band count
	void SetBand(int band);

	// Set the tolereance for the residual to converge under
	void SetTolerance(double tol);

	// Solve using the Jacobi Method
	void SolveJacobi();

	// Solve using the Gauss-Seidel Method
	void SolveGaussSeidel();

	// Print the result vector
	void WriteResultToFile();

	// Generate the gnuplot file
	void Write2DContourToGnuplot(double daxial, double dother);


private:

	double** m_A; // coefficient matrix
	double* m_b;  // RHS matrix
	double* m_x;  // solution vector

	int m_problemDimension; // dimension of A (number of grid points) 
	int m_band; // number of bands in A

	int m_endIter; // number of maximum iterations
	int m_numNodesAxial; // number of nodes in axial direction 

	double m_tolerance; // tolerance value for the residual

	sjsSolverHelpers *m_helperFunctions;

	double* m_x_after; // solution at next time step
	double* m_dummy;   // dummy vector to be used in varius routines
	double* m_dummy2;  // second dummy vector to be used in various routines 


};




#endif // _SJSSOLVERS_H_
