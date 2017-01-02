
#include "sjsSolvers.h"

/// Constructor
sjsSolvers::sjsSolvers(){

	// Default settings
	m_helperFunctions = new sjsSolverHelpers();
	m_endIter = 50000;
	m_tolerance = 1e-4;

}

/// Destructor
sjsSolvers::~sjsSolvers(){

}

/// Set the coefficient matrix - second entry is true if the matrix is to be printed
void sjsSolvers::SetCoefficientMatrix(double **A, bool printToFile){
	m_A = A;

	if(printToFile){
		std::ofstream debug;
		debug.open("Amatrix.txt");
		for (int i = 0; i<m_problemDimension; i++){
			for (int j = 0; j<m_band; j++){
				debug<<m_A[j][i]<<"	";
			}
			debug<<std::endl;
		}
		debug.close();

	}

}


/// Set the solution vector for the current iteration
void sjsSolvers::SetSolutionVector(double* x){
	for(int i=0; i<m_problemDimension; i++){
		m_x[i] = x[i];
	}
}

/// Get the solution vector for the current iteration
double* sjsSolvers::GetSolutionVector(){
	return m_x;
	//delete [] m_x;
}


/// Set the RHS
void sjsSolvers::SetRHS(double* b, bool printToFile){
	m_b = b;

	if(printToFile){
		std::ofstream debug;
		debug.open("bmatrix.txt");
		for (int i = 0; i<m_problemDimension; i++){
			debug<<m_b[i]<<std::endl;
		}
		debug.close();

	}	

}

/// Set number of grid points
void sjsSolvers::SetProblemDimension(int size){
	m_problemDimension = size;
	// fill the solver array with zeros
	m_x = AllocateDynamicVector<double>(size);

	// fill these dummy vectors that will be used by the solver routines
	m_x_after = AllocateDynamicVector<double>(size);
	m_dummy = AllocateDynamicVector<double>(size);
	m_dummy2 = AllocateDynamicVector<double>(size);
}

/// Set number of points in axial direction
void sjsSolvers::SetNumNodesAxial(int nodes){
	m_numNodesAxial = nodes;

}

/// Set coefficient matrix band count
void sjsSolvers::SetBand(int band){
	m_band = band;
}

/// Set the tolereance for the residual to converge under
void sjsSolvers::SetTolerance(double tol){
	m_tolerance = tol;

}

/// Solve using the Jacobi Method
void sjsSolvers::SolveJacobi(){

	int num_iter = 0;

	double residual = 4.0; // an initial residual value
	//double *x_after = AllocateDynamicVector<double>(m_problemDimension);
	//double *dummy = AllocateDynamicVector<double>(m_problemDimension);
	//double *dummy2 = AllocateDynamicVector<double>(m_problemDimension);

	m_helperFunctions->VectorZeroFill(m_problemDimension, m_x_after);
	m_helperFunctions->VectorZeroFill(m_problemDimension, m_dummy);
	m_helperFunctions->VectorZeroFill(m_problemDimension, m_dummy2);

	//debug<<"Starting the while loop..."<<std::endl;	
	while (residual > m_tolerance){
		if (m_band == 5){
			for(int i =0; i<m_problemDimension; i++){
				//*********First row of the Matrix*********//
				if (i == 0){
					m_x_after[i] = (m_b[i] - m_x[i+1]*m_A[3][i] - m_x[i+m_numNodesAxial]*m_A[4][i])/m_A[2][i];
				}

				//*********Second part the Matrix (m_A[0] is out of bounds)*********//
				if (i>0 && i<m_numNodesAxial){
					m_x_after[i] = (m_b[i] - m_x[i-1]*m_A[1][i] - m_x[i+1]*m_A[3][i] - m_x[i+m_numNodesAxial]*m_A[4][i])/m_A[2][i];
				}
				//*********Interior of the matrix*********//
				if (i>= m_numNodesAxial && i<m_problemDimension-m_numNodesAxial){
					m_x_after[i] = (m_b[i] - m_x[i-1]*m_A[1][i] - m_x[i+1]*m_A[3][i] - m_x[i+m_numNodesAxial]*m_A[4][i] - m_x[i-m_numNodesAxial]*m_A[0][i])/m_A[2][i];
				}
				//*********Fourth part the Matrix (m_A[4] is out of bounds)*********//
				if (i>=m_problemDimension-m_numNodesAxial && i<m_problemDimension-1){
					m_x_after[i] = (m_b[i] - m_x[i-1]*m_A[1][i] - m_x[i+1]*m_A[3][i] - m_x[i-m_numNodesAxial]*m_A[0][i])/m_A[2][i];
				}
				//*********Last row of the Matrix*********//
				if (i == m_problemDimension-1){
					m_x_after[i] = (m_b[i] - m_x[i-1]*m_A[1][i] - m_x[i-m_numNodesAxial]*m_A[0][i])/m_A[2][i];
				}

			}


		}

		if (m_band == 3){
			for (int i = 0; i<m_problemDimension; i++){
				//*********First row of the Matrix*********//
				if (i == 0){
					m_x_after[i] = (m_b[i] - m_x[i+1]*m_A[2][i])/m_A[1][i];
				}
				//*********Interior of the Matrix*********//
				if (i>0 && i<m_problemDimension-1){
					m_x_after[i] = (m_b[i] - m_x[i+1]*m_A[2][i] - m_x[i-1]*m_A[0][i])/m_A[1][i];	
				}
				//*********Last row of the Matrix*********//
				if (i == m_problemDimension-1){
					m_x_after[i] = (m_b[i] - m_x[i-1]*m_A[0][i])/m_A[1][i];	
				}


			}


		}

		m_helperFunctions->VectorZeroFill(m_problemDimension, m_dummy2);
		m_helperFunctions->TwoVectorCopy(m_problemDimension,m_x_after,m_x);
		m_helperFunctions->BandedMatrixVectorMultiplication(m_A, m_x_after, m_dummy, m_numNodesAxial, m_problemDimension, m_band);
		m_helperFunctions->TwoVectorSubtraction(m_problemDimension, m_b, m_dummy, m_dummy2);
		residual = m_helperFunctions->VectorTwoNorm(m_problemDimension, m_dummy2);
		//debug<<residual<<std::endl;
		num_iter++;
		if (num_iter == m_endIter){
			std::cout<<"Jacobi ended due to maxiter. Res: "<<residual<<std::endl;
			break;

		}


	}


}

/// Solve using the Gauss-Seidel Method
void sjsSolvers::SolveGaussSeidel(){

	int num_iter = 0;

	double residual = 4.0; // an initial residual value

	m_helperFunctions->VectorZeroFill(m_problemDimension, m_x_after);
	m_helperFunctions->VectorZeroFill(m_problemDimension, m_dummy);
	m_helperFunctions->VectorZeroFill(m_problemDimension, m_dummy2);

	//debug<<"Starting the while loop..."<<std::endl;	
	while (residual > m_tolerance){
		if (m_band == 5){
			for(int i =0; i<m_problemDimension; i++){
				//*********First row of the Matrix*********//
				if (i == 0){
					m_x_after[i] = (m_b[i] - m_x[i+1]*m_A[3][i] - m_x[i+m_numNodesAxial]*m_A[4][i])/m_A[2][i];
				}

				//*********Second part the Matrix (m_A[0] is out of bounds)*********//
				if (i>0 && i<m_numNodesAxial){
					m_x_after[i] = (m_b[i] - m_x_after[i-1]*m_A[1][i] - m_x[i+1]*m_A[3][i] - m_x[i+m_numNodesAxial]*m_A[4][i])/m_A[2][i];
				}
				//*********Interior of the matrix*********//
				if (i>= m_numNodesAxial && i<m_problemDimension-m_numNodesAxial){
					m_x_after[i] = (m_b[i] - m_x_after[i-1]*m_A[1][i] - m_x[i+1]*m_A[3][i] - m_x[i+m_numNodesAxial]*m_A[4][i] - m_x_after[i-m_numNodesAxial]*m_A[0][i])/m_A[2][i];
				}
				//*********Fourth part the Matrix (m_A[4] is out of bounds)*********//
				if (i>=m_problemDimension-m_numNodesAxial && i<m_problemDimension-1){
					m_x_after[i] = (m_b[i] - m_x_after[i-1]*m_A[1][i] - m_x[i+1]*m_A[3][i] - m_x_after[i-m_numNodesAxial]*m_A[0][i])/m_A[2][i];
				}
				//*********Last row of the Matrix*********//
				if (i == m_problemDimension-1){
					m_x_after[i] = (m_b[i] - m_x_after[i-1]*m_A[1][i] - m_x[i-m_numNodesAxial]*m_A[0][i])/m_A[2][i];
				}

			}


		}

		if (m_band == 3){
			for (int i = 0; i<m_problemDimension; i++){
				//*********First row of the Matrix*********//
				if (i == 0){
					m_x_after[i] = (m_b[i] - m_x[i+1]*m_A[2][i])/m_A[1][i];
				}
				//*********Interior of the Matrix*********//
				if (i>0 && i<m_problemDimension-1){
					m_x_after[i] = (m_b[i] - m_x[i+1]*m_A[2][i] - m_x_after[i-1]*m_A[0][i])/m_A[1][i];	
				}
				//*********Last row of the Matrix*********//
				if (i == m_problemDimension-1){
					m_x_after[i] = (m_b[i] - m_x_after[i-1]*m_A[0][i])/m_A[1][i];	
				}


			}


		}

		m_helperFunctions->VectorZeroFill(m_problemDimension, m_dummy2);
		m_helperFunctions->TwoVectorCopy(m_problemDimension,m_x_after,m_x);
		m_helperFunctions->BandedMatrixVectorMultiplication(m_A, m_x_after, m_dummy, m_numNodesAxial, m_problemDimension, m_band);
		m_helperFunctions->TwoVectorSubtraction(m_problemDimension, m_b, m_dummy, m_dummy2);
		residual = m_helperFunctions->VectorTwoNorm(m_problemDimension, m_dummy2);
		//debug<<residual<<std::endl;
		num_iter++;
		if (num_iter == m_endIter){
			std::cout<<"Gauss-Seidel ended due to maxiter. Res: "<<residual<<std::endl;
			break;

		}


	}	

}

/// Print the result vector
void sjsSolvers::WriteResultToFile(){
	std::ofstream debug;
	debug.open("Solution.txt");
	for (int i = 0; i<m_problemDimension; i++){
		debug<<m_x[i]<<std::endl;
	}
	debug.close();

}

/// Generate the gnuplot file
void sjsSolvers::Write2DContourToGnuplot(double daxial, double dother){
	int ydim = m_problemDimension / m_numNodesAxial ;
	int k = 0;
	std::ofstream debug;
	debug.open("results.dat");
	for (int j = 0; j<ydim; j++){
		for (int i = 0; i<m_numNodesAxial; i++){
			debug<<i*daxial<<" "<<j*dother<<" "<<m_x[k]<<std::endl;
			k++;
		}

	}
	debug.close();

}