#ifndef _SJSNEWTONIANJET_H_
#define _SJSNEWTONIANJET_H_


/*

	Slender Jet Solver 

	Top boundary is at z=0

	Bottom boundary is initially at z = end of the domain
	But if the profile is adaptive, that can be changed


*/


#include "sjsCommon.h"
#include "sjsSolvers.h"

class sjsNewtonianJet{

public:

	// Constructor
	sjsNewtonianJet();
	// Destructor
	~sjsNewtonianJet();

	// Set the number of nodes in axial direction
	void SetNumNodesX(int nodeCount);

	// Get the number of nodes in x direction
	int GetNumNodesX();

	// Set the number of time iteration count
	void SetTimeIterationCount(int t);

	// Set the domain length in axial direction
	void SetAxialDomainLength(double length);

	// Set the timestep of the problem
	void SetTimeStep(double dt);

	// Get the timestep of the problem
	double GetTimeStep();

	// Set whether top H boundary is neumann
	void SetHTopNeumann(bool state);

	// Set whether bottom H boundary is neumann
	void SetHBottomNeumann(bool state);

	// Set whether top U boundary is neumann
	void SetUTopNeumann(bool state);

	// Set whether bottom U boundary is neumann
	void SetUBottomNeumann(bool state);

	// Build the H matrix
	void BuildHMatrix();

	// Initiate the matrices
	void InitiateMatrices();

	// Calculate curvature field
	void CalculateCurvatureField();
	void CalculateCurvatureFieldParallel();

	// Set jet profile perturbation
	void SetJetProfilePerturbation(double eps);

	// Initiate Jet Profile
	void InitiateJetProfile();

	// Solve the H matrix
	void SolveHMatrix();

	// Set previous h field
	void SetPreviousThicknessField();

	// Set previous u field
	void SetPreviousVelocityField();

	// Set initial radius
	void SetInitialRadius(double radius);

	// Assign the surface tension of the fluid
	void SetSurfaceTension(double gamma);

	// Assign the density of the fluid
	void SetFluidDensity(double value);

	// Assign the kinematic viscosity of the fluid
	void SetKinematicViscosity(double value);

	// Build U Matrix
	void BuildUMatrix();

	// Solve U Matrix
	void SolveUMatrix();

	// Inititate U solver
	void InitiateUSolver();

	// Initiate H solver
	void InitiateHSolver();

	// Set solver in time
	void SetIsTimeExplicit(bool state);

	// Set whether solver works in parallel
	void SetIsSolverParallel(bool state);

	// Get whether solver works in parallel
	bool GetIsSolverParallel();

	// Set if the results will be written to a file
    void SetIsWriteResultsToFile(bool write_results_to_file);
 	bool GetIsWriteResultsToFile();

  	// Set the time step interval for the results to be written to a file (valid if write_results_to_file = true)
    void SetResultsIntervalFile(int results_interval_file);
    int GetResultsIntervalFile();

    // Set the time step interval for the results to be written on the screem 
    void SetResultsIntervalScreen(int results_interval_screen);
  	int GetResultsIntervalScreen();

  	// Set solver maximum timestep
    void SetSolverMaxTimestep(int solver_max_timestep);
    int GetSolverMaxTimestep();

	// Solve U explicitly
	void SolveUExplicit();
	void SolveUExplicitParallel();

	// Solve H explicitly
	void SolveHExplicit();
	void SolveHExplicitParallel();

	// Initiate explicit solvers
	void InitiateExplicitSolvers();

	// Explicit time march
	void ExplicitTimeMarch(int timestep); 


private:

	int m_numNodesX; // Number of nodes in X direction
	int m_timeIterationCount; // iteration count 

	double m_xlength; // length of the domain in x direction
	double m_dt; // time step
	double m_dz; // mesh size - in axial direction


	double m_htopboundary; // h value at the top boundary
	double m_hbottomboundary; // h value at the bottom boundary

	bool m_isHtopneumann;
	bool m_isHbottomneumann;
	bool m_isUtopneumann;
	bool m_isUbottomneumann;

	bool m_isTimeExplicit; // true if the solution is explicit in time

	bool m_isSolverParallel; // true if te solver works in parallel

	bool m_isResultsWriteFile; // true if output files will be generated
	int m_resultsFileTimestepInterval; // timestep interval for results file to be generated
	int m_resultsScreenTimestepInterval; // timestep interval for results to be printed on the screen
	int m_maxtimestep; // maximum timestep for the solver

	double** m_Ah; // coefficient matrix for jet thickness
	double* m_bh; // right hand side for jet thickness
	double* m_xh; // solution vector for jet thickness
	double* m_xh_old; // solution vector for jet thickness at previous time step

	double** m_Au; // coefficient matrix for axial velocity
	double* m_bu; // right hand side for axial velocity
	double* m_xu; // solution vector for axial velocity
	double* m_xu_old; // solution vector for axial velocity at previous time step

	double* m_kappa; // curvature field - calculated using h field

	double* m_uexp; // explicit u solution
	double* m_hexp; // explicit h solution
	double* m_uexp_old;  // explicit u solution at previous time step
	double* m_hexp_old;  // explicit h solution at previous time step



	double m_gamma; // surface tension value
	double m_rho; // density of the fluid
	double m_nu; // kinematic viscosity of the fluid

	double m_uTopDrichlet; // top bc for velocity
	double m_uBottomDrichlet; // bottom bc for velocity

	double m_r0; // initial radius
	double m_eps; // perturbation on the initial jet profile

	sjsSolvers *m_Hsolver;
	sjsSolvers* m_Usolver;


	// arrays to be used in parallel solver
	double* m_h2; 
	double* m_dv2dz2; 
	double* m_dv_dz;
	double* m_dh_dz;
	double* m_inertia_term; 
	double* m_pressure_term; 
	double* m_diff_term; 
	double* m_f;
	
};



#endif // _SJSNEWTONIANJET_H_
