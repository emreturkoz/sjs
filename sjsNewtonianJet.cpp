#include "sjsNewtonianJet.h"


/*

	How to use this solver?
	1) Set the number of nodes
	2) Set the axial domain length
	3) Set the time step
	4) Initiate Matrices
	5) Set the Jet Profile
	6) Set top and bottom bc's (whether neumann of not)
	7) Initiate the Jet Profile
	8) Calculate the curvature field
	9) Build H matrix
	
	10) Assign fluid properties (surface tension, kinematic viscosity and density)
	11) Build U matrix

	Then start a loop

	***************Loop***************
	Solve U Matrix
	Solve H Matrix

	Build H Matrix
	Build U Matrix
	***************Loop***************

*/


/// Constructor
sjsNewtonianJet::sjsNewtonianJet(){

	m_Usolver = new sjsSolvers();
	m_Hsolver = new sjsSolvers();
	m_isHtopneumann = false; 	// default: drichlet BC for top H 
	m_isHbottomneumann = false; // default: drichlet BC for top H

	m_uTopDrichlet = 0.0;         // default: not moving boundaries for U
	m_uBottomDrichlet= 0.0;       // default: not moving boundaries for U
	m_timeIterationCount = 0;

	m_r0 = 0.01; // default initial radius (m)
	m_eps = 0.0; // initially the perturbation on the jet profile is zero


	m_isTimeExplicit = false; // default setting for time march
}

/// Destructor
sjsNewtonianJet::~sjsNewtonianJet(){

}


/// Set solver in time
void sjsNewtonianJet::SetIsTimeExplicit(bool state){
	m_isTimeExplicit = state;
}

/// Set the number of nodes in axial direction
void sjsNewtonianJet::SetNumNodesX(int nodeCount){
	m_numNodesX = nodeCount;

}


/// Get the number of nodes in x direction
int sjsNewtonianJet::GetNumNodesX(){
	return m_numNodesX;

}

// Set the number of time iteration count
void sjsNewtonianJet::SetTimeIterationCount(int t){
	m_timeIterationCount = t;
}

/// Set the domain length in axial direction
void sjsNewtonianJet::SetAxialDomainLength(double length){
	m_xlength = length;

	// Initiate the curvature and mesh arrays 
	m_kappa = AllocateDynamicVector<double>(m_numNodesX);
	m_dz = m_xlength/((double)m_numNodesX-1);	

}

/// Set initial radius
void sjsNewtonianJet::SetInitialRadius(double radius){
	m_r0 = radius;
}

/// Set the timestep of the problem
void sjsNewtonianJet::SetTimeStep(double dt){

	double tr = sqrt(m_rho*m_r0*m_r0*m_r0/m_gamma);

	m_dt = dt*tr;


}

/// Get the timestep of the problem
double sjsNewtonianJet::GetTimeStep(){
	return m_dt;
}


/// Set whether top H boundary is neumann
void sjsNewtonianJet::SetHTopNeumann(bool state){
	m_isHtopneumann = state;
}

/// Set whether bottom H boundary is neumann
void sjsNewtonianJet::SetHBottomNeumann(bool state){
	m_isHbottomneumann = state;
}

/// Set whether top U boundary is neumann
void sjsNewtonianJet::SetUTopNeumann(bool state){
	m_isUtopneumann = state;
}

/// Set whether bottom U boundary is neumann
void sjsNewtonianJet::SetUBottomNeumann(bool state){
	m_isUbottomneumann = state;
}


/// Initiate the matrices
void sjsNewtonianJet::InitiateMatrices(){
	m_Ah = AllocateDynamic2DArray<double>(3,m_numNodesX);
	m_bh = AllocateDynamicVector<double>(m_numNodesX);
	m_xh = AllocateDynamicVector<double>(m_numNodesX);
	m_xh_old = AllocateDynamicVector<double>(m_numNodesX);

	m_Au = AllocateDynamic2DArray<double>(3,m_numNodesX-1);
	m_bu = AllocateDynamicVector<double>(m_numNodesX-1);
	m_xu = AllocateDynamicVector<double>(m_numNodesX-1);
	m_xu_old = AllocateDynamicVector<double>(m_numNodesX-1);
}



/// Build the H matrix
void sjsNewtonianJet::BuildHMatrix(){
	
	double v_avg;


	//*******Boundary conditions*******//

	if(m_isHtopneumann){
		m_Ah[1][0] = 1.0;
		m_Ah[2][0] = -1.0;
		m_bh[0] = 0;
	}
	else{
		m_Ah[1][0] = 1.0;
		m_bh[0] = m_xh[0];
	}


	if (m_isHbottomneumann){
		m_Ah[1][m_numNodesX-1] = 1.0;
		m_Ah[0][m_numNodesX-1] = -1.0;
		m_bh[m_numNodesX-1] = 0;
	}
	else{
		m_Ah[1][m_numNodesX-1] = 1.0;
		m_bh[m_numNodesX-1] = m_xh[m_numNodesX-1];		
	}



	// Interior of the matrix
	for (int i=1; i<m_numNodesX-1; i++){
		v_avg = 0.50*(m_xu_old[i-1] + m_xu_old[i]);
		m_Ah[0][i] = -v_avg*0.5/m_dz;
		m_Ah[2][i] = v_avg*0.5/m_dz;
		m_Ah[1][i] = (1.0/m_dt) + 0.5*(m_xu_old[i] - m_xu_old[i-1])/m_dz;
		m_bh[i] = m_xh_old[i]/m_dt;
	}

}



/// Set jet profile perturbation
void sjsNewtonianJet::SetJetProfilePerturbation(double eps){
	m_eps = eps; 
}

/// Calculate curvature field
void sjsNewtonianJet::CalculateCurvatureField(){
	double dh_dzz, dh_dz;


	for(int i = 1; i<m_numNodesX-1; i++){
		if (m_isTimeExplicit){
			dh_dzz = (m_hexp[i+1] - 2*m_hexp[i] + m_hexp[i-1])/(m_dz*m_dz);
			dh_dz =  (m_hexp[i+1] - m_hexp[i-1])/(2*m_dz);						
			m_kappa[i] = 1.0/(m_hexp[i]*sqrt(1+dh_dz*dh_dz)) - dh_dzz/pow(1+dh_dz*dh_dz,1.5);			
		}
		else{
			dh_dzz = (m_xh[i+1] - 2*m_xh[i] + m_xh[i-1])/(m_dz*m_dz);
			dh_dz =  (m_xh[i+1] - m_xh[i-1])/(2*m_dz);
			m_kappa[i] = 1.0/(m_xh[i]*sqrt(1+dh_dz*dh_dz)) - dh_dzz/pow(1+dh_dz*dh_dz,1.5);			
		}
		
	}


	//m_kappa[0] = 1.0/m_xh[0];
	//m_kappa[m_numNodesX-1] = 1.0/m_xh[m_numNodesX-1];

	m_kappa[0] = m_kappa[1];
	m_kappa[m_numNodesX-1] = m_kappa[m_numNodesX-2];


	/*
	std::ofstream debug;
	debug.open("CurvatureProfile.txt");
	for(int i = 0; i<m_numNodesX; i++){
		debug<<i*m_dz<<" "<<m_kappa[i]<<std::endl;
	}

	debug.close();
	*/


}



/// Initiate Jet Profile
void sjsNewtonianJet::InitiateJetProfile(){
	double k = (2*PI/m_xlength);
	//double epsilon = 0.01;


	for(int i = 0; i<m_numNodesX; i++){
		if(m_isTimeExplicit){
			m_hexp[i] = m_r0*(1+m_eps*cos(k*i*m_dz));
			m_hexp_old[i] = m_hexp[i];
		}
		else{
			m_xh[i] = m_r0*(1+m_eps*cos(k*i*m_dz));
			m_xh_old[i] = m_xh[i];	
		}
		
	}


	/*
	std::ofstream debug;
	debug.open("JetProfile.txt");
	for(int i = 0; i<m_numNodesX; i++){
		debug<<i*m_dz<<" "<<m_xh[i]<<std::endl;
	}

	debug.close();
	*/


}


/// Assign the surface tension of the fluid
void sjsNewtonianJet::SetSurfaceTension(double gamma){
	m_gamma = gamma;
}

/// Assign the density of the fluid
void sjsNewtonianJet::SetFluidDensity(double value){
	m_rho = value;
}

/// Assign the kinematic viscosity of the fluid
void sjsNewtonianJet::SetKinematicViscosity(double value){
	m_nu = value;
}

/// Set if the results will be written to a file
void sjsNewtonianJet::SetIsWriteResultsToFile(bool write_results_to_file){
	m_isResultsWriteFile = write_results_to_file;
}

/// Set the time step interval for the results to be written to a file (valid if write_results_to_file = true)
void sjsNewtonianJet::SetResultsIntervalFile(int results_interval_file){
	m_resultsFileTimestepInterval = results_interval_file;
}

/// Set the time step interval for the results to be written on the screen
void sjsNewtonianJet::SetResultsIntervalScreen(int results_interval_screen){
	m_resultsScreenTimestepInterval = results_interval_screen;
}

/// Set solver maximum timestep
void sjsNewtonianJet::SetSolverMaxTimestep(int solver_max_timestep){
	m_maxtimestep = solver_max_timestep;
}

/// Get if the results will be written to a file
bool sjsNewtonianJet::GetIsWriteResultsToFile(){
	return m_isResultsWriteFile;
}

/// Get the time step interval for the results to be written to a file (valid if write_results_to_file = true)
int sjsNewtonianJet::GetResultsIntervalFile(){
	return m_resultsFileTimestepInterval;
}

/// Get the time step interval for the results to be written on the screen
int sjsNewtonianJet::GetResultsIntervalScreen(){
	return m_resultsScreenTimestepInterval;
}

/// Get solver maximum timestep
int sjsNewtonianJet::GetSolverMaxTimestep(){
	return m_maxtimestep;
}


/// Build U Matrix
void sjsNewtonianJet::BuildUMatrix(){
	double dV = m_dz;
	double Fe, Fw; // coefficients at the velocity formula
	double dkappa; // gradient of the curvature
	double h, h2; // effective h^2 value averaged on velocity grid
	double dh2_dz; // the derivative of the h2 
	double Dw, De, Dc; // other coefficients coming from the viscous term


	this->CalculateCurvatureField();


	// Interior of the matrix
	for (int i=1; i<m_numNodesX-2; i++){

		Fe = 0.5*(m_xu_old[i+1] + m_xu_old[i])*dV/m_dz;
		Fw = 0.5*(m_xu_old[i-1] + m_xu_old[i])*dV/m_dz;
		dkappa = (m_kappa[i+1] - m_kappa[i])/m_dz;

		h = 0.5*(m_xh[i+1] + m_xh[i]);
		h2 = h*h;
		dh2_dz = (m_xh[i+1]*m_xh[i+1]  - m_xh[i]*m_xh[i])/m_dz;

		De = -(3.0*m_nu/h2)*(0.5*h2/(m_dz*m_dz)) - ((3.0*m_nu/h2)*(0.5*dh2_dz/m_dz)); 
		Dw = -(3.0*m_nu/h2)*(0.5*h2/(m_dz*m_dz)) + ((3.0*m_nu/h2)*(0.5*dh2_dz/m_dz));
		Dc = (3.0*m_nu/h2)*(h2/(m_dz*m_dz));


		m_Au[0][i]= - 0.5*Fw + Dw;
		m_Au[2][i]= 0.5*Fe + De;
		m_Au[1][i]= (0.5*Fe - 0.5*Fw) + (dV/m_dt) + Dc;
		m_bu[i] = (m_xu_old[i]*dV/m_dt) - (dV*m_gamma*dkappa/m_rho);

	}

	// Boundary conditions
	if (m_isUtopneumann){
		m_Au[1][0] = 1.0;
		m_Au[2][0] = -1.0;
		m_bu[0] = 0;		
	}
	else{
		m_Au[1][0] = 1.0;
		m_bu[0] = m_uTopDrichlet;
	}

	if (m_isUbottomneumann){
		m_Au[1][m_numNodesX-2] = 1.0;
		m_Au[0][m_numNodesX-2] = -1.0;
		m_bu[m_numNodesX-2] = 0;		
	}
	else{
		m_Au[1][m_numNodesX-2] = 1.0;
		m_bu[m_numNodesX-2] = m_uBottomDrichlet;
	}

}


/// Inititate U solver
void sjsNewtonianJet::InitiateUSolver(){
	m_Usolver->SetProblemDimension(m_numNodesX-1);
	m_Usolver->SetBand(3);
	m_Usolver->SetNumNodesAxial(m_numNodesX-1);
	m_Usolver->SetTolerance(0.0001);


}

/// Initiate H solver
void sjsNewtonianJet::InitiateHSolver(){
	m_Hsolver->SetProblemDimension(m_numNodesX);
	m_Hsolver->SetBand(3);
	m_Hsolver->SetNumNodesAxial(m_numNodesX);
	m_Hsolver->SetTolerance(0.0001);


}


/// Solve U Matrix
void sjsNewtonianJet::SolveUMatrix(){
	m_Usolver->SetCoefficientMatrix(m_Au, false);
	m_Usolver->SetRHS(m_bu, false);
	m_Usolver->SetSolutionVector(m_xu);

	m_Usolver->SolveGaussSeidel();


	m_xu = m_Usolver->GetSolutionVector();

}

/// Solve the H matrix
void sjsNewtonianJet::SolveHMatrix(){
	m_Hsolver->SetCoefficientMatrix(m_Ah, false);
	m_Hsolver->SetRHS(m_bh, false);
	m_Hsolver->SetSolutionVector(m_xh);

	m_Hsolver->SolveGaussSeidel();

	m_xh = m_Hsolver->GetSolutionVector();



	//m_solver->WriteResultToFile();

}


/// Solve U explicitly
void sjsNewtonianJet::SolveUExplicit(){

	this->CalculateCurvatureField();

	double f,diff_term, h2, dv2dz2, dv_dz, dh_dz;
	double inertia_term, pressure_term;
	int i;

	//#pragma omp parallel
	//{
		//#pragma omp for
		for(i=1; i<m_numNodesX-1; i++){
			h2 = m_hexp_old[i]*m_hexp_old[i];
			dv2dz2 = (m_uexp_old[i+1] - 2.0*m_uexp_old[i] + m_uexp_old[i-1])/(m_dz*m_dz); 		
			dv_dz = 0.5*(m_uexp_old[i+1] - m_uexp_old[i-1])/m_dz;
			dh_dz = 0.5*(m_hexp_old[i+1] - m_hexp_old[i-1])/m_dz; 
			inertia_term = m_uexp_old[i]*0.5*(m_uexp_old[i+1] - m_uexp_old[i-1])/m_dz;
			pressure_term = -(m_gamma/m_rho)*0.5*(m_kappa[i+1] - m_kappa[i-1])/m_dz;
			diff_term = h2*(dv2dz2) + dv_dz*2.0*m_hexp_old[i]*dh_dz;  
			f = -inertia_term + pressure_term + 3.0*diff_term*m_nu/h2 ;
			m_uexp[i] = m_uexp_old[i] + m_dt*f;
		}
	//}


		// Boundary conditions
	if (m_isUtopneumann){
		m_uexp[0] = m_uexp[1]; 	
	}
	else{
		m_uexp[0] = m_uTopDrichlet;
	}

	if (m_isUbottomneumann){
		m_uexp[m_numNodesX-1] = m_uexp[m_numNodesX-2]; 
	}
	else{
 		m_uexp[m_numNodesX-1] = m_uBottomDrichlet;
	}


}

/// Solve H explicitly
void sjsNewtonianJet::SolveHExplicit(){

	double f;

	for(int i = 1; i<m_numNodesX-1; i++){
		f = (-m_hexp_old[i]/2.0)*0.5*(m_uexp_old[i+1] - m_uexp_old[i-1])/m_dz - (m_uexp_old[i]*0.5*(m_hexp_old[i+1] - m_hexp_old[i-1])/m_dz);
		m_hexp[i] = m_hexp_old[i] + m_dt*f;
	}

	//*******Boundary conditions*******//

	if(m_isHtopneumann){
		m_hexp[0] = m_hexp[1];
	}

	if (m_isHbottomneumann){
		m_hexp[m_numNodesX-1] = m_hexp[m_numNodesX-2];
	}


}

/// Initiate explicit solvers
void sjsNewtonianJet::InitiateExplicitSolvers(){
	// Initiate both of the solution arrays
	m_uexp = AllocateDynamicVector<double>(m_numNodesX);
	m_hexp = AllocateDynamicVector<double>(m_numNodesX);

	// Initiate the arrays for the previous time step
	m_uexp_old = AllocateDynamicVector<double>(m_numNodesX);
	m_hexp_old = AllocateDynamicVector<double>(m_numNodesX);



}

/// Explicit time march
void sjsNewtonianJet::ExplicitTimeMarch(int timestep){

	// copy the previous time step
	for(int i = 0; i<m_numNodesX; i++){
		m_hexp_old[i] = m_hexp[i];
		m_uexp_old[i] = m_uexp[i];
	}

	m_timeIterationCount = timestep;

	if (m_timeIterationCount % m_resultsFileTimestepInterval == 0){
		char fileName[80];
		char fileName2[80];
		sprintf( fileName, "u-timestep%i.txt", m_timeIterationCount );
		sprintf( fileName2, "h-timestep%i.txt", m_timeIterationCount );

		std::ofstream debug;
		std::ofstream debug2;
		debug.open(fileName);
		debug2.open(fileName2);

		for(int i=0; i<m_numNodesX; i++){
			debug<<m_uexp_old[i]<<std::endl;
			debug2<<m_hexp_old[i]<<std::endl;
		}
		debug.close();
		debug2.close();
	}


}



/// Set previous h field
void sjsNewtonianJet::SetPreviousThicknessField(){
	for(int i =0; i<m_numNodesX; i++){
		m_xh_old[i] = m_xh[i];
	}

	if (m_timeIterationCount % m_resultsFileTimestepInterval == 0){
		char fileName[80];
		sprintf( fileName, "h-timestep%i.txt", m_timeIterationCount );

		std::ofstream debug;
		debug.open(fileName);

		for(int i=0; i<m_numNodesX; i++)
			debug<<m_xh_old[i]<<std::endl;

		debug.close();
	}



}

/// Set previous u field
void sjsNewtonianJet::SetPreviousVelocityField(){
	for(int i =0; i<m_numNodesX; i++){
		m_xu_old[i] = m_xu[i];
	}

	if (m_timeIterationCount % m_resultsFileTimestepInterval == 0){
		char fileName[80];
		sprintf( fileName, "u-timestep%i.txt", m_timeIterationCount );

		std::ofstream debug;
		debug.open(fileName);

		for(int i=0; i<m_numNodesX; i++)
			debug<<m_xu_old[i]<<std::endl;

		debug.close();
	}


}










