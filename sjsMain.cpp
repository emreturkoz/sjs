
/*
	Started on December 16, 2016
	Developed by Emre Turkoz

	This is initially developed to solve slender jet equations in 1D
	SJS stands for Slender Jet Solver
	
*/


#include "sjsCommon.h"
#include "sjsModel.h"



// Execute this if the routine is implicit in time
int ImplicitRoutine(sjsModel* implicit_model){
  int max_timestep = implicit_model->NewtonianGetSolverMaxTimestep();
  int screen_interval = implicit_model->NewtonianGetResultsIntervalFile();
  double current_timestep;



  implicit_model->NewtonianBuildHMatrix();
  implicit_model->NewtonianBuildUMatrix();

  implicit_model->NewtonianInitiateUSolver();
  implicit_model->NewtonianInitiateHSolver();

  // The loop starts here
  for(int i = 0; i<max_timestep; i++){

    implicit_model->NewtonianSetTimeIterationCount(i);

    implicit_model->NewtonianSolveUMatrix();
    implicit_model->NewtonianSolveHMatrix();

    implicit_model->NewtonianSetPreviousThicknessField();
    implicit_model->NewtonianSetPreviousVelocityField();

    implicit_model->NewtonianBuildHMatrix();
    implicit_model->NewtonianBuildUMatrix();
    
    if(i % screen_interval == 0){
      current_timestep = (double)i*implicit_model->NewtonianGetTimeStep();
      std::cout<<std::setprecision(5)<<current_timestep<<std::endl;
    }

  }



  return 0;
}


// Execute this if the routine is explicit in time
int ExplicitRoutine(sjsModel* explicit_model){

  int max_timestep = explicit_model->NewtonianGetSolverMaxTimestep();
  int screen_interval = explicit_model->NewtonianGetResultsIntervalScreen();
  bool is_parallel = explicit_model->NewtonianGetIsSolverParallel();
  double current_timestep;




  // The loop starts here
  
  if (is_parallel == false){
    for(int i = 0; i<max_timestep; i++){
      explicit_model->NewtonianSolveHExplicit();
      explicit_model->NewtonianSolveUExplicit();
      explicit_model->NewtonianExplicitTimeMarch(i);

      if(i % screen_interval == 0){
        current_timestep = (double)i*explicit_model->NewtonianGetTimeStep();
        std::cout<<i<<" "<<std::setprecision(5)<<current_timestep<<std::endl;
      }

    }
  }
  else{
    for(int i = 0; i<max_timestep; i++){
      explicit_model->NewtonianSolveHExplicitParallel();
      explicit_model->NewtonianSolveUExplicitParallel();
      explicit_model->NewtonianExplicitTimeMarch(i);

      if(i % screen_interval == 0){
        current_timestep = (double)i*explicit_model->NewtonianGetTimeStep();
        std::cout<<i<<" "<<std::setprecision(5)<<current_timestep<<std::endl;
      }

    }

  }



  return 0;
}

int main(){
  
  sjsModel *model = new sjsModel();

  // List of parameters to describe the problem and to be sent to the model
  
  int axial_node_count = 1400;
  double axial_domain_length = 0.04; // [m]
  double initial_jet_radius = 0.005; // [m]
  double jet_profile_perturbation = 0.08; // relative -> perturbation/radius
  double surface_tension = 0.0729; // N/m -> 72.9 mN/m is for water
  double density = 1000; // [kg/m3] density of the fliuid -> 1000 kg/m3 is for water
  double kinematic_viscosity = 1.0e-6;
  double time_step = 3.0e-5; // relative time step. this can be multiplied either with Rayleigh time step, or diffusive time step

  bool is_time_explicit = true;
  bool is_solver_parallel = true;


  // boundary condition
  bool u_top_neumann = true;
  bool u_bottom_neumann = true;
  bool h_top_neumann = true;
  bool h_bottom_neumann = true;


  // post-processing parameters
  bool write_results_to_file = true;
  int results_interval_file = 10000;
  int results_interval_screen = 10000;
  int solver_max_timestep = 200000;

  


  //******Commands for the solution of Laplace's equation******//
  //model->LaplaceSetNumNodesX(25);
  //model->LaplaceSetNumNodesY(25);
  //model->LaplaceSetDomainLengthX(1.0);
  //model->LaplaceSetDomainLengthY(1.0);
  //model->LaplaceBuildMatrix();
  //model->LaplaceSolve();


  // Commands for the solution of Newtonian Slender Jet Equations
  model->NewtonianSetNumNodesX(axial_node_count);
  model->NewtonianSetAxialDomainLength(axial_domain_length);

  model->NewtonianSetInitialRadius(initial_jet_radius);
  model->NewtonianSetJetProfilePerturbation(jet_profile_perturbation);
  
  model->NewtonianSetIsSolverParallel(is_solver_parallel);
  
  // This is also where the arrays are intitated.
  model->NewtonianSetIsTimeExplicit(is_time_explicit); 

  model->NewtonianInitateJetProfile();
  model->NewtonianCalculateCurvatureField();

  model->NewtonianSetSurfaceTension(surface_tension);
  model->NewtonianSetFluidDensity(density);
  model->NewtonianSetKinematicViscosity(kinematic_viscosity);
  model->NewtonianSetTimeStep(time_step);

  model->NewtonianSetUTopNeumann(u_top_neumann);
  model->NewtonianSetUBottomNeumann(u_bottom_neumann);
  model->NewtonianSetHTopNeumann(h_top_neumann);
  model->NewtonianSetHBottomNeumann(h_bottom_neumann);

  model->NewtonianSetIsWriteResultsToFile(write_results_to_file);
  model->NewtonianSetResultsIntervalFile(results_interval_file);
  model->NewtonianSetResultsIntervalScreen(results_interval_screen);
  model->NewtonianSetSolverMaxTimestep(solver_max_timestep);


  // Informative output
  double t_rayleigh = sqrt(density*initial_jet_radius*initial_jet_radius*initial_jet_radius/surface_tension);
  double t_diffusion = initial_jet_radius*initial_jet_radius/kinematic_viscosity;
  std::cout<<"The Rayleigh time scale for this problem: "<<t_rayleigh<<" seconds"<<std::endl;
  std::cout<<"The diffusion time scale is: "<<t_diffusion<<" seconds"<<std::endl; 
  double u_rayleigh = initial_jet_radius/t_rayleigh;
  double t_cfl = (axial_domain_length/(axial_node_count-1))/u_rayleigh;
  std::cout<<"The CFL condition time scale is: "<<t_cfl<<std::endl;


  //ImplicitRoutine(model);
  ExplicitRoutine(model);





  std::cout<<"All is well"<<std::endl;
  return 0;
}
