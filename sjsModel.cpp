#include "sjsModel.h"


/// Constructor
sjsModel::sjsModel(){
	// initiating all the models in the code
	m_laplaceModel = new sjsLaplaceModel();
	m_newtonianJet = new sjsNewtonianJet();
	m_viscoelasticJet = new sjsViscoelasticJet();

}

/// Destructor
sjsModel::~sjsModel(){


}

//**********Functions for the Laplace Model**********//

/// Set the number of nodes in x direction for the Laplace model
void sjsModel::LaplaceSetNumNodesX(int nodes){
	m_laplaceModel->SetNumNodesX(nodes);

}

/// Set the number of nodes in y direction for the Laplace model
void sjsModel::LaplaceSetNumNodesY(int nodes){
	m_laplaceModel->SetNumNodesY(nodes);
}

/// Set domain length in x direction for the Laplace Model
void sjsModel::LaplaceSetDomainLengthX(double length){
	m_laplaceModel->SetDomainLengthX(length);
}

/// Set domain length in y direction for the Laplace Model
void sjsModel::LaplaceSetDomainLengthY(double length){
	m_laplaceModel->SetDomainLengthY(length);
}

/// Build the banded sparse matrix through second order finite differencing for the Laplace model
void sjsModel::LaplaceBuildMatrix(){
	m_laplaceModel->BuildMatrix();
}

///  Solve the linear system for the Laplace Model
void sjsModel::LaplaceSolve(){
	m_laplaceModel->Solve();

}

//**********Functions for the Newtonian Slender Jet Model**********//

/// Set the number of nodes in axial direction for Newtonian Slender Jet model
void sjsModel::NewtonianSetNumNodesX(int nodeCount){
	m_newtonianJet->SetNumNodesX(nodeCount);
}

/// Get the number of nodes in x direction for Newtonian Slender Jet model
int sjsModel::NewtonianGetNumNodesX(){
	return m_newtonianJet->GetNumNodesX();
}

/// Set the number of time iteration count
void sjsModel::NewtonianSetTimeIterationCount(int t){
	m_newtonianJet->SetTimeIterationCount(t);
}

/// Set the domain length in axial direction for Newtonian Slender Jet model
void sjsModel::NewtonianSetAxialDomainLength(double length){
	m_newtonianJet->SetAxialDomainLength(length);
}

/// Set the timestep of the problem for Newtonian Slender Jet model
void sjsModel::NewtonianSetTimeStep(double dt){
	m_newtonianJet->SetTimeStep(dt);
}

/// Get the timestep of the problem for Newtonian Slender Jet model
double sjsModel::NewtonianGetTimeStep(){
	return m_newtonianJet->GetTimeStep();
}

/// Set whether top H boundary is neumann for Newtonian Slender Jet model
void sjsModel::NewtonianSetHTopNeumann(bool state){
	m_newtonianJet->SetHTopNeumann(state);
}

/// Set whether bottom H boundary is neumann for Newtonian Slender Jet model
void sjsModel::NewtonianSetHBottomNeumann(bool state){
	m_newtonianJet->SetHBottomNeumann(state);
}

/// Set whether top U boundary is neumann
void sjsModel::NewtonianSetUTopNeumann(bool state){
	m_newtonianJet->SetUTopNeumann(state);
}

/// Set whether bottom U boundary is neumann
void sjsModel::NewtonianSetUBottomNeumann(bool state){
	m_newtonianJet->SetUBottomNeumann(state);
}

/// Set initial radius
void sjsModel::NewtonianSetInitialRadius(double radius){
	m_newtonianJet->SetInitialRadius(radius);
}


/// Build the H matrix for Newtonian Slender Jet model
void sjsModel::NewtonianBuildHMatrix(){
	m_newtonianJet->BuildHMatrix();
}

/// Solve the H matrix
void sjsModel::NewtonianSolveHMatrix(){
	m_newtonianJet->SolveHMatrix();
}

/// Set jet profile perturbation
void sjsModel::NewtonianSetJetProfilePerturbation(double eps){
	m_newtonianJet->SetJetProfilePerturbation(eps);
}

/// Initiate the matrices for Newtonian Slender Jet model
void sjsModel::NewtonianInitiateMatrices(){
	m_newtonianJet->InitiateMatrices();
}

/// Calculate curvature field for Newtonian Slender Jet model
void sjsModel::NewtonianCalculateCurvatureField(){
	m_newtonianJet->CalculateCurvatureField();
}

/// Initiate Jet Profile for Newtonian Slender Jet model
void sjsModel::NewtonianInitateJetProfile(){
	m_newtonianJet->InitiateJetProfile();
}


/// Assign the surface tension of the fluid
void sjsModel::NewtonianSetSurfaceTension(double gamma){
	m_newtonianJet->SetSurfaceTension(gamma);
}

/// Assign the density of the fluid
void sjsModel::NewtonianSetFluidDensity(double value){
	m_newtonianJet->SetFluidDensity(value);
}

/// Assign the kinematic viscosity of the fluid
void sjsModel::NewtonianSetKinematicViscosity(double value){
	m_newtonianJet->SetKinematicViscosity(value);
}

/// Set previous h field
void sjsModel::NewtonianSetPreviousThicknessField(){
	m_newtonianJet->SetPreviousThicknessField();
}

/// Set previous u field
void sjsModel::NewtonianSetPreviousVelocityField(){
	m_newtonianJet->SetPreviousVelocityField();
}

/// Inititate U solver
void sjsModel::NewtonianInitiateUSolver(){
	m_newtonianJet->InitiateUSolver();
}

/// Initiate H solver
void sjsModel::NewtonianInitiateHSolver(){
	m_newtonianJet->InitiateHSolver();
}

/// Build U Matrix
void sjsModel::NewtonianBuildUMatrix(){
	m_newtonianJet->BuildUMatrix();
}

/// Solve U Matrix
void sjsModel::NewtonianSolveUMatrix(){
	m_newtonianJet->SolveUMatrix();
}


/// Set solver in time
void sjsModel::NewtonianSetIsTimeExplicit(bool state){
	m_newtonianJet->SetIsTimeExplicit(state);

	if (state){
		this->NewtonianInitiateExplicitSolvers();
	}
	else{
		this->NewtonianInitiateMatrices();
	}
}

/// Solve U explicitly
void sjsModel::NewtonianSolveUExplicit(){
	m_newtonianJet->SolveUExplicit();
}

/// Solve H explicitly
void sjsModel::NewtonianSolveHExplicit(){
	m_newtonianJet->SolveHExplicit();
}

/// Initiate explicit solvers
void sjsModel::NewtonianInitiateExplicitSolvers(){
	m_newtonianJet->InitiateExplicitSolvers();
}

/// Explicit time march
void sjsModel::NewtonianExplicitTimeMarch(int timestep){
	m_newtonianJet->ExplicitTimeMarch(timestep);

}

