all:
	g++ sjsMain.cpp sjsNewtonianJet.cpp sjsViscoelasticJet.cpp sjsLaplaceModel.cpp sjsSolvers.cpp sjsModel.cpp sjsSolverHelpers.cpp -fopenmp -lrt -o main
