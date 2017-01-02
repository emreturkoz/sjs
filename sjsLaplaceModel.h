#ifndef _SJSLAPLACEMODEL_H_
#define _SJSLAPLACEMODEL_H_

#include "sjsCommon.h"
#include "sjsSolvers.h"

/*
This is the model for the laplace's equation in 2D cartesian coordinate system.
Implemented for solver check. Easy convergence and easy implementation.
*/

class sjsLaplaceModel {
public:
	sjsLaplaceModel(); // Constructor
	~sjsLaplaceModel(); // Destructor

	// Set the number of nodes in x direction
	void SetNumNodesX(int nodeCount);

	// Set the number of nodes in y direction
	void SetNumNodesY(int nodeCount);

	// Get the number of nodes in x direction
	int GetNumNodesX();

	// Get the number of nodes in y direction
	int GetNumNodesY();

	// Boolean states for Neumann boundary condition
	void SetNeumannBCTop(bool state);
	void SetNeumannBCBottom(bool state);
	void SetNeumannBCLeft(bool state);
	void SetNeumannBCRight(bool state);


	// Drichlet BC values
	void SetDrichletTopValue(double value);
	void SetDrichletBottomValue(double value);
	void SetDrichletLeftValue(double value);
	void SetDrichletRightValue(double value);


	// Set domain length in x direction
	void SetDomainLengthX(double length);


	// Set domain length in y direction
	void SetDomainLengthY(double length);

	// Build the banded sparse matrix through second order finite differencing
	void BuildMatrix();

	// Solve the linear system
	void Solve();




private:

	int m_numNodesX; // Number of nodes in X direction
	int m_numNodesY; // Number of nodes in Y direction

	bool m_neumannTop;    // boolean for Neumann bc: top
	bool m_neumannBottom; // boolean for Neumann bc: bottom
	bool m_neumannLeft;   // boolean for Neumann bc: left
	bool m_neumannRight;  // boolean for Neumann bc: right

	double m_drichletTop;     // drichlet bc for top - valid if m_neumann is false
	double m_drichletBottom;  // drichlet bc for bottom - valid if m_neumann is false
	double m_drichletLeft;    // drichlet bc for left - valid if m_neumann is false
	double m_drichletRight;	  // drichlet bc for right - valid if m_neumann is false 

	double m_xlength; // length of the domain in x direction
	double m_ylength; // length of the domain in y direction

	double** m_A; // coefficient matrix
	double* m_b; // right hand side
	double* m_x; // solution vector 

	sjsSolvers *m_solver;

};


#endif // _SJSLAPLACEMODEL_H_