#include <iostream>
#include "FEM1D.h"
#include <vector>

int main() {

	// Spatial paramater
	FEM1DDD sol; // 1 for 1D; 2 for 2D; 3 for 3D;

	sol.obj_start = 0; // Object start point

	sol.obj_end = 0.1; // Object end point

	sol.area = 1e-4; // Cross-sectional area

	sol.elasticity = 1e11; // Modulus of elasticity

	sol.forfun = 1e11; // body force

	sol.bc1 = 0.; // left end dirichlet

	sol.bc2 = 0.001; // right end dirichlet

	sol.conv_criteria = 0.00000001; // convergence criteria

	sol.basisorder = 1; // Order of Basis Function

	char unimesh = 'y'; // uniform mesh = yes
//	char unimesh = 'n'; // uniform mesh = no

	if (unimesh == 'y') {
		sol.numofele = 5; // user paramater for mesh
		sol.Gen_Mesh();
		sol.Dis_Mesh();
		sol.Gen_kvec();
		sol.Dis_kvec();
		sol.Mat_ludc();
//		sol.Dis_MtLU();
		sol.Gen_forf();
		sol.Gen_gvec();
		sol.Sol_GaSi();
	}

}