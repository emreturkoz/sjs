
/*
	Started on December 16, 2016
	Developed by Emre Turkoz

	This is initially developed to solve slender jet equations in 1D
	SJS stands for Slender Jet Solver
	
*/


#include "sjsCommon.h"
#include "sjsModel.h"


int ImplicitRoutine(sjsModel* implicit_model){
  int max_timestep = 4000000;  // 10000000;
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
    
    if(i % 4000 == 0){
      current_timestep = (double)i*implicit_model->NewtonianGetTimeStep();
      std::cout<<std::setprecision(5)<<current_timestep<<std::endl;
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
  double jet_profile_perturbation = 0.05; // relative -> perturbation/radius
  double surface_tension = 0.0729; // N/m -> 72.9 mN/m is for water
  double density = 1000; // [kg/m3] density of the fliuid -> 1000 kg/m3 is for water
  double kinematic_viscosity = 1.0e-6;
  double time_step = 2.0e-5; // relative time step. this can be multiplied either with Rayleigh time step, or diffusive time step

  bool is_time_explicit = true;


  // boundary condition
  bool u_top_neumann = true;
  bool u_bottom_neumann = true;
  bool h_top_neumann = true;
  bool h_bottom_neumann = true;




  


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

  //ImplicitRoutine(model);






  std::cout<<"All is well"<<std::endl;
  return 0;
}
