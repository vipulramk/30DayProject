#include <iostream>
#include "FEM1D.h"
#include <vector>
#include <ctime>

int main() {

	int start_clock = clock();

	std::cout << "Execution started at: " << start_clock << std::endl;

	// Spatial paramater
	FEM1DDD sol; // FEM<Spatial Dimensions><Dir-Dir or Dir-Nue>

	sol.obj_start = 0; // Object start point

	sol.obj_end = 0.1; // Object end point

	sol.area = 1e-4; // Cross-sectional area

	sol.elasticity = 1e11; // Modulus of elasticity

	sol.forfun = 1e11; // body force

	sol.bc1 = 0.; // left end dirichlet

	sol.bc2 = 0.001; // right end dirichlet

	sol.conv_criteria = 0.00001; // convergence criteria

	sol.basisorder = 2; // Order of Basis Function

	char unimesh = 'y'; // uniform mesh = yes
//	char unimesh = 'n'; // uniform mesh = no

	if (unimesh == 'y') {
		sol.numofele = 200; // user paramater for mesh
		sol.Gen_Mesh();
		sol.Dis_Mesh();
		sol.Gen_kvec();
		sol.Dis_kvec();
		sol.Gen_forf();
		sol.Gen_gvec();
//		sol.Sol_GaSi();
//		sol.Sol_ludc();
		sol.Sol_Guel();
//		sol.Dis_MtLU();
		sol.Dis_dvec();
	}

	int stop_clock = clock();

	std::cout << "Execution finished at: " << stop_clock << std::endl;

	double time_clock = (stop_clock - start_clock)/double(CLOCKS_PER_SEC)*1000;

	std::cout << "Execution time: " << time_clock << "ms" << std::endl;

return 0;
}
