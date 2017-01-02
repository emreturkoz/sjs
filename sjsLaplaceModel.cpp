/*
This is the model for the laplace's equation in 2D cartesian coordinate system.
Implemented for solver check. Easy convergence and easy implementation.
*/

#include "sjsLaplaceModel.h"

/// Constructor
sjsLaplaceModel::sjsLaplaceModel(){

	// Initially all the Neumann BC is false, Left BC is 1.0, all the other BC's are 0.0

	m_neumannTop = false;
	m_neumannBottom = false;
	m_neumannLeft = false;
	m_neumannRight = false;

	m_drichletTop = 0.0;
	m_drichletBottom = 0.0;
	m_drichletLeft = 1.0;
	m_drichletRight = 0.0;

	// Initiate the solver
	m_solver = new sjsSolvers();



}

/// Destructor
sjsLaplaceModel::~sjsLaplaceModel(){

}

/// Set the number of nodes in x direction
void sjsLaplaceModel::SetNumNodesX(int nodeCount){
	m_numNodesX = nodeCount;
}

/// Set the number of nodes in y direction
void sjsLaplaceModel::SetNumNodesY(int nodeCount){
	m_numNodesY = nodeCount;
}

/// Get the number of nodes in x direction
int sjsLaplaceModel::GetNumNodesX(){
	return m_numNodesX;
}

/// Get the number of nodes in y direction
int sjsLaplaceModel::GetNumNodesY(){
	return m_numNodesY;
}

/// Set the boolean for Neumann boundary condition: Top
void sjsLaplaceModel::SetNeumannBCTop(bool state){
	m_neumannTop = state;
}

/// Set the boolean for Neumann boundary condition: Bottom
void sjsLaplaceModel::SetNeumannBCBottom(bool state){
	m_neumannBottom = state;

}

/// Set the boolean for Neumann boundary condition: Left
void sjsLaplaceModel::SetNeumannBCLeft(bool state){
	m_neumannLeft = state;
}

/// Set the boolean for Neumann boundary condition: Right
void sjsLaplaceModel::SetNeumannBCRight(bool state){
	m_neumannRight = state;
}


/// Set Drichlet BC value: Top
void sjsLaplaceModel::SetDrichletTopValue(double value){
	m_drichletTop = value;
}

/// Set Drichlet BC value: Bottom
void sjsLaplaceModel::SetDrichletBottomValue(double value){
	m_drichletBottom = value;
}

/// Set Drichlet BC value: Left
void sjsLaplaceModel::SetDrichletLeftValue(double value){
	m_drichletLeft = value;

}

/// Set Drichlet BC value: Right
void sjsLaplaceModel::SetDrichletRightValue(double value){
	m_drichletRight = value;
}

/// Set domain length in x direction
void sjsLaplaceModel::SetDomainLengthX(double length){
	m_xlength = length;

}

/// Set domain length in y direction
void sjsLaplaceModel::SetDomainLengthY(double length){
	m_ylength = length;
}


/// Build the banded sparse matrix through second order finite differencing
void sjsLaplaceModel::BuildMatrix(){
	// initiate the matrices
	// 5-point stencil for 2D Cartesian coordinate system
	m_A = AllocateDynamic2DArray<double>(5,m_numNodesX*m_numNodesY);
	m_b = AllocateDynamicVector<double>(m_numNodesX*m_numNodesY);
	m_x = AllocateDynamicVector<double>(m_numNodesX*m_numNodesY);

	// calculate the distance between each point
	double dx, dy;
	dx = m_xlength/((double)m_numNodesX - 1);
	dy = m_ylength/((double)m_numNodesY - 1);


	// for loop integers
	int i,j,k;

	//******boundary conditions******//
	
	// Bottom boundary
	j = 0;
	for(i=0; i<m_numNodesX; i++){
		k = j*m_numNodesX+i;

		if(m_neumannBottom){
			m_A[2][k] = 1.0;
			m_A[4][k] = -1.0;
			m_b[k] = 0.0;
		}
		else{
			m_A[2][k] = 1.0;
			m_b[k] = m_drichletBottom;
		}

	}

	// Top boundary
	j = m_numNodesY-1;
	for(i=0; i<m_numNodesX; i++){
		k = j*m_numNodesX+i;

		if(m_neumannTop){
			m_A[2][k] = 1.0;
			m_A[0][k] = -1.0;
			m_b[k] = 0.0;
		}
		else{
			m_A[2][k] = 1.0;
			m_b[k] = m_drichletTop;
		}

	}	

	// Left Boundary
	i = 0;
	for (j=0; j<m_numNodesY; j++){
		k = j*m_numNodesX+i;

		if(m_neumannLeft){
			m_A[2][k]=1.0;
			m_A[3][k]=-1.0;
			m_b[k] = 0.0;
		}
		else{
			m_A[2][k] = 1.0;
			m_b[k] = m_drichletLeft;
		}

	}


	// Right Boundary
	i = m_numNodesX-1;
	for (j=0; j<m_numNodesY; j++){
		k = j*m_numNodesX + i;

		if(m_neumannRight){
			m_A[2][k]=1.0;
			m_A[1][k]=-1.0;
			m_b[k] = 0.0;	
		}
		else{
			m_A[2][k] = 1.0;
			m_b[k] = m_drichletRight; 
		}

	}


	// Inside the domain
	for (j=1; j<m_numNodesY-1; j++){
		for (i=1; i<m_numNodesX-1; i++){
			k = j*m_numNodesX + i;
			m_A[2][k] = -2.0/pow(dx,2.0) - 2.0/pow(dy,2.0);
			m_A[1][k] = 1.0/pow(dx,2.0);
			m_A[3][k] = 1.0/pow(dx,2.0);
			m_A[0][k] = 1.0/pow(dy,2.0);
			m_A[4][k] = 1.0/pow(dy,2.0);
			m_b[k] = 0.0;
		}

	}

}



/// Solve the linear system
void sjsLaplaceModel::Solve(){

	m_solver->SetProblemDimension(m_numNodesX*m_numNodesY);
	
	m_solver->SetBand(5);
	m_solver->SetCoefficientMatrix(m_A, true);
	m_solver->SetRHS(m_b, true);
	m_solver->SetNumNodesAxial(m_numNodesX);

	m_solver->SetTolerance(0.0001);
	//m_solver->SolveJacobi();
	m_solver->SolveGaussSeidel();

	m_solver->WriteResultToFile();

	m_solver->Write2DContourToGnuplot(m_xlength/((double)m_numNodesX - 1), m_ylength/((double)m_numNodesY - 1));



}






