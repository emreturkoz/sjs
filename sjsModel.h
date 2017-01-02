#ifndef _SJSMODEL_H_
#define _SJSMODEL_H_


#include "sjsLaplaceModel.h"
#include "sjsNewtonianJet.h"
#include "sjsViscoelasticJet.h"

class sjsModel {

public:
	sjsModel();
	~sjsModel();

	//**********Functions for the Laplace Model**********//
	
	// Set the number of nodes in x direction for the Laplace model
	void LaplaceSetNumNodesX(int nodes);

	// Set the number of nodes in y direction for the Laplace model
	void LaplaceSetNumNodesY(int nodes);

	// Set domain length in x direction for the Laplace Model
	void LaplaceSetDomainLengthX(double length);

	// Set domain length in y direction for the Laplace Model
	void LaplaceSetDomainLengthY(double length);

	// Build the banded sparse matrix through second order finite differencing for the Laplace model
	void LaplaceBuildMatrix();

	// Solve the linear system for the Laplace Model
	void LaplaceSolve();


	//**********Functions for the Newtonian Slender Jet Model**********//

	// Set the number of nodes in axial direction
	void NewtonianSetNumNodesX(int nodeCount);

	// Get the number of nodes in x direction
	int NewtonianGetNumNodesX();

	// Set the number of time iteration count
	void NewtonianSetTimeIterationCount(int t);

	// Set the domain length in axial direction
	void NewtonianSetAxialDomainLength(double length);

	// Set the timestep of the problem
	void NewtonianSetTimeStep(double dt);

	// Get the timestep of the problem
	double NewtonianGetTimeStep();

	// Set whether top H boundary is neumann
	void NewtonianSetHTopNeumann(bool state);

	// Set whether bottom H boundary is neumann
	void NewtonianSetHBottomNeumann(bool state);

	// Set whether top U boundary is neumann
	void NewtonianSetUTopNeumann(bool state);

	// Set whether bottom U boundary is neumann
	void NewtonianSetUBottomNeumann(bool state);

	// Set initial radius
	void NewtonianSetInitialRadius(double radius);	

	// Build the H matrix
	void NewtonianBuildHMatrix();

	// Solve the H matrix
	void NewtonianSolveHMatrix();	

	// Set jet profile perturbation
	void NewtonianSetJetProfilePerturbation(double eps);

	// Initiate the matrices
	void NewtonianInitiateMatrices();

	// Calculate curvature field
	void NewtonianCalculateCurvatureField();

	// Initiate Jet Profile
	void NewtonianInitateJetProfile();

	// Assign the surface tension of the fluid
	void NewtonianSetSurfaceTension(double gamma);

	// Assign the density of the fluid
	void NewtonianSetFluidDensity(double value);

	// Assign the kinematic viscosity of the fluid
	void NewtonianSetKinematicViscosity(double value);

	// Set previous h field
	void NewtonianSetPreviousThicknessField();

	// Set previous u field
	void NewtonianSetPreviousVelocityField();

	// Inititate U solver
	void NewtonianInitiateUSolver();

	// Initiate H solver
	void NewtonianInitiateHSolver();

	// Build U Matrix
	void NewtonianBuildUMatrix();

	// Solve U Matrix
	void NewtonianSolveUMatrix();

	// Set solver in time
	void NewtonianSetIsTimeExplicit(bool state);

	// Solve U explicitly
	void NewtonianSolveUExplicit();

	// Solve H explicitly
	void NewtonianSolveHExplicit();

	// Initiate explicit solvers
	void NewtonianInitiateExplicitSolvers();

	// Explicit time march
	void NewtonianExplicitTimeMarch(int timestep); 



private:
	sjsLaplaceModel *m_laplaceModel;
	sjsNewtonianJet *m_newtonianJet;
	sjsViscoelasticJet *m_viscoelasticJet;

};


#endif // _SJSMODEL_H_
