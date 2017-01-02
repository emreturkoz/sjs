#include "sjsViscoelasticJet.h"


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
sjsViscoelasticJet::sjsViscoelasticJet(){

	m_Usolver = new sjsSolvers();
	m_Hsolver = new sjsSolvers();

	m_srrSolver = new sjsSolvers();
	m_szzSolver = new sjsSolvers();

	m_isHtopneumann = false; 	// default: drichlet BC for top H 
	m_isHbottomneumann = false; // default: drichlet BC for top H

	m_uTopDrichlet = 0.0;         // default: not moving boundaries for U
	m_uBottomDrichlet= 0.0;       // default: not moving boundaries for U
	m_timeIterationCount = 0;

	m_r0 = 0.01; // default initial radius (m)

	m_modelNo = 1; // default model no: 1 - Oldroyd-B


}

/// Destructor
sjsViscoelasticJet::~sjsViscoelasticJet(){

}


/// Set the number of nodes in axial direction
void sjsViscoelasticJet::SetNumNodesX(int nodeCount){
	m_numNodesX = nodeCount;

}


/// Get the number of nodes in x direction
int sjsViscoelasticJet::GetNumNodesX(){
	return m_numNodesX;

}

// Set the number of time iteration count
void sjsViscoelasticJet::SetTimeIterationCount(int t){
	m_timeIterationCount = t;
}

/// Set the domain length in axial direction
void sjsViscoelasticJet::SetAxialDomainLength(double length){
	m_xlength = length;

}

/// Set initial radius
void sjsViscoelasticJet::SetInitialRadius(double radius){
	m_r0 = radius;
}

/// Set the timestep of the problem
void sjsViscoelasticJet::SetTimeStep(double dt){

	double tr = sqrt(m_rho*m_r0*m_r0*m_r0/m_gamma);

	m_dt = dt*tr;


}

/// Get the timestep of the problem
double sjsViscoelasticJet::GetTimeStep(){
	return m_dt;
}


/// Set whether top H boundary is neumann
void sjsViscoelasticJet::SetHTopNeumann(bool state){
	m_isHtopneumann = state;
}

/// Set whether bottom H boundary is neumann
void sjsViscoelasticJet::SetHBottomNeumann(bool state){
	m_isHbottomneumann = state;
}

/// Set whether top U boundary is neumann
void sjsViscoelasticJet::SetUTopNeumann(bool state){
	m_isUtopneumann = state;
}

/// Set whether bottom U boundary is neumann
void sjsViscoelasticJet::SetUBottomNeumann(bool state){
	m_isUbottomneumann = state;
}


/// Initiate the matrices
void sjsViscoelasticJet::InitiateMatrices(){

	// Initiate matrices & vectors for jet thickness, H

	m_Ah = AllocateDynamic2DArray<double>(3,m_numNodesX);
	m_bh = AllocateDynamicVector<double>(m_numNodesX);
	m_xh = AllocateDynamicVector<double>(m_numNodesX);
	m_xh_old = AllocateDynamicVector<double>(m_numNodesX);

	// Initiate matrices & vectors for axial jet velocity, U
	m_Au = AllocateDynamic2DArray<double>(3,m_numNodesX-1);
	m_bu = AllocateDynamicVector<double>(m_numNodesX-1);
	m_xu = AllocateDynamicVector<double>(m_numNodesX-1);
	m_xu_old = AllocateDynamicVector<double>(m_numNodesX-1);

	// Curvature
	m_kappa = AllocateDynamicVector<double>(m_numNodesX);

	// Initiate matrices & vectors for sigma_rr
	m_Arr  = AllocateDynamic2DArray<double>(3,m_numNodesX);
	m_brr = AllocateDynamicVector<double>(m_numNodesX);
	m_xrr = AllocateDynamicVector<double>(m_numNodesX);
	m_xrr_old = AllocateDynamicVector<double>(m_numNodesX);

	// Initiate matrices & vectors for sigma_zz
	m_Azz = AllocateDynamic2DArray<double>(3,m_numNodesX);
	m_bzz = AllocateDynamicVector<double>(m_numNodesX);
	m_xzz = AllocateDynamicVector<double>(m_numNodesX);
	m_xzz_old = AllocateDynamicVector<double>(m_numNodesX);

	m_dz = m_xlength/((double)m_numNodesX-1);

}



/// Build the H matrix
void sjsViscoelasticJet::BuildHMatrix(){
	
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





/// Calculate curvature field
void sjsViscoelasticJet::CalculateCurvatureField(){
	double dh_dzz, dh_dz;


	for(int i = 1; i<m_numNodesX-1; i++){
		dh_dzz = (m_xh[i+1] - 2*m_xh[i] + m_xh[i-1])/(2*m_dz*m_dz);
		dh_dz =  (m_xh[i+1] - m_xh[i-1])/(2*m_dz);
		m_kappa[i] = 1.0/(m_xh[i]*sqrt(1+dh_dzz*dh_dzz)) - dh_dzz/pow(1+dh_dz*dh_dz,1.5);
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
void sjsViscoelasticJet::InitiateJetProfile(){
	double k = (2*PI/m_xlength);
	double epsilon = 0.01;

	for(int i = 0; i<m_numNodesX; i++){
		m_xh[i] = m_r0*(1+epsilon*cos(k*i*m_dz));
		m_xh_old[i] = m_xh[i];
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
void sjsViscoelasticJet::SetSurfaceTension(double gamma){
	m_gamma = gamma;
}

/// Assign the density of the fluid
void sjsViscoelasticJet::SetFluidDensity(double value){
	m_rho = value;
}

/// Assign the kinematic viscosity of the fluid
void sjsViscoelasticJet::SetKinematicViscosity(double value){
	m_nu = value;
}




/// Build U Matrix
void sjsViscoelasticJet::BuildUMatrix(){
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
void sjsViscoelasticJet::InitiateUSolver(){
	m_Usolver->SetProblemDimension(m_numNodesX-1);
	m_Usolver->SetBand(3);
	m_Usolver->SetNumNodesAxial(m_numNodesX-1);
	m_Usolver->SetTolerance(0.0001);


}

/// Initiate H solver
void sjsViscoelasticJet::InitiateHSolver(){
	m_Hsolver->SetProblemDimension(m_numNodesX);
	m_Hsolver->SetBand(3);
	m_Hsolver->SetNumNodesAxial(m_numNodesX);
	m_Hsolver->SetTolerance(0.0001);


}


/// Solve U Matrix
void sjsViscoelasticJet::SolveUMatrix(){
	m_Usolver->SetCoefficientMatrix(m_Au, false);
	m_Usolver->SetRHS(m_bu, false);
	m_Usolver->SetSolutionVector(m_xu);

	m_Usolver->SolveGaussSeidel();

	m_xu = m_Usolver->GetSolutionVector();

}

/// Solve the H matrix
void sjsViscoelasticJet::SolveHMatrix(){
	m_Hsolver->SetCoefficientMatrix(m_Ah, false);
	m_Hsolver->SetRHS(m_bh, false);
	m_Hsolver->SetSolutionVector(m_xh);

	m_Hsolver->SolveGaussSeidel();

	m_xh = m_Hsolver->GetSolutionVector();



	//m_solver->WriteResultToFile();

}


/// Set previous h field
void sjsViscoelasticJet::SetPreviousThicknessField(){
	for(int i =0; i<m_numNodesX; i++){
		m_xh_old[i] = m_xh[i];
	}

	if (m_timeIterationCount % 50000 == 0){
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
void sjsViscoelasticJet::SetPreviousVelocityField(){
	for(int i =0; i<m_numNodesX; i++){
		m_xu_old[i] = m_xu[i];
	}

	if (m_timeIterationCount % 50000 == 0){
		char fileName[80];
		sprintf( fileName, "u-timestep%i.txt", m_timeIterationCount );

		std::ofstream debug;
		debug.open(fileName);

		for(int i=0; i<m_numNodesX; i++)
			debug<<m_xu_old[i]<<std::endl;

		debug.close();
	}


}


/// Set model no
void sjsViscoelasticJet::SetModelNo(int no){
	m_modelNo = no;


}

/// Set polymer viscosity
void sjsViscoelasticJet::SetPolymerViscosity(double visco){
	m_polymerViscosity = visco;
}

/// Set relaxation time
void sjsViscoelasticJet::SetRelaxationTime(double time){
	m_relaxationtime = time;
}


/// Build sigma_rr matrix
void sjsViscoelasticJet::BuildSigmaRRMatrix(){


}

/// Build sigma_zz matrix
void sjsViscoelasticJet::BuildSigmaZZMatrix(){

	double v_avg; // average velocity 

	// Interior of the matrix
	for (int i=1; i<m_numNodesX-1; i++){
		v_avg = 0.50*(m_xu_old[i-1] + m_xu_old[i]);
		m_Azz[0][i] = -v_avg*m_relaxationtime*0.5/m_dz;
		m_Azz[2][i] = v_avg*m_relaxationtime*0.5/m_dz;
		m_Azz[1][i] = 1.0 + (m_relaxationtime/m_dt) + 2.0*(m_xu_old[i] - m_xu_old[i-1])/m_dz;
		m_bzz[i] = (m_relaxationtime*m_xzz_old[i]/m_dt) + 2.0*m_polymerViscosity*(m_xu_old[i] - m_xu_old[i-1])/m_dz;
	}

}







