#include <iostream>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <vector>

class FEM1DDD {
public:

//	FEM1D ();
//	~FEM1D ();

	// Member Functions:
	void Gen_Mesh();
	void Dis_Mesh();
	void Gen_kvec();
	void Dis_kvec();
	void Mat_ludc();
	void Dis_MtLU();
	void Gen_forf();
	void Gen_gvec();
	void Sol_GaSi();

	// Member Objects:
	double obj_start, obj_end, lenofele;
	double area, elasticity, forfun, bc1, bc2, conv_criteria;
	int numofele;
	int basisorder;
	int klen,kwid;

	// Member vectors:
	std::vector<double> h;
	std::vector<std::vector<double> > k_vector;
	std::vector<std::vector<double> > l_vector;
	std::vector<std::vector<double> > u_vector;
	std::vector<double> f_vector;
	std::vector<double> g_vector;
	std::vector<double> d_vector;

};

void FEM1DDD::Gen_Mesh() {

	lenofele = obj_end - obj_start;

	for (int i=0; i<numofele; i++) {
		h.push_back(lenofele/numofele);
	}

	return;
}

void FEM1DDD::Dis_Mesh() {

	std::cout << "The mesh is as follows: " << std::endl;

	for (int i=0; i<h.size(); i++) {
		std::cout << h[i] << "\t";
	}
	std::cout << std::endl;

	return;
}

void FEM1DDD::Gen_kvec() {

	if (basisorder == 1) { //--------------------------1st order basis function

		klen = numofele-1;
		kwid = numofele-1;

  		k_vector.resize(klen);
		for (int i = 0; i < klen; ++i){
			k_vector[i].resize(kwid);
		}   

		for (int i=0; i<klen; i++) {

			for (int j=0; j<kwid; j++) {
				
				int check = std::abs(i - j);

				//std::cout<<check<<std::endl;

				switch (check) {
					case 0:
					k_vector[i][j] = elasticity*((1/h[i]) + (1/h[i+1]));
					break;
					case 1:
					k_vector[i][j] = elasticity*(-1/h[(i+j+1)/2]);
//					std::cout<<"element: " << (i+j+1)/2 << "\t" << i  << "\t" << j << std::endl;
					break;
					default:
					k_vector[i][j] = 0;
					break;
				}
			}			
		}
	}

	if (basisorder == 2) { //--------------------------2nd order basis function

		klen = (2*numofele)-1;
		kwid = (2*numofele)-1;

  		k_vector.resize(klen);
		for (int i = 0; i < klen; ++i){
			k_vector[i].resize(kwid);
		}   

		for (int i=0; i<klen; i++) {

			for (int j=0; j<kwid; j++) {
				
				int check = std::abs(i - j);

				//std::cout<<check<<std::endl;
				int k = -1;

				switch (check) {
					case 0:
					if (fmod(i,2) == 0)	k_vector[i][j] = elasticity*8/(3*h[i/2]);
					else if (fmod(i,2) == 1) k_vector[i][j] = elasticity*7*((1/h[(i-1)/2]) + (1/h[(i+1)/2]))/6;
//					if (fmod(i,2) == 0)	k_vector[i][j] = 8/(3*h[i/2]);
//					else if (fmod(i,2) == 1) k_vector[i][j] = 7*((1/h[(i-1)/2]) + (1/h[(i+1)/2]))/6;
					break;
					case 1:
					for (k = -1; (2+4*(k)) < i+j; k++) {

					}
					k_vector[i][j] = -4*elasticity/(3*h[k]);
//					k_vector[i][j] = -4/(3*h[k+1]);
//					std::cout<<"element: " << k << "\t" << i  << "\t" << j << std::endl;
					break;
					case 2:
					if (fmod(j,2) == 1)	{
						k_vector[i][j] = elasticity*1/(6*h[(i+j)/4]);
//						if (fmod(j,2) == 1)	k_vector[i][j] = 1./6;
//						std::cout<<"element 1/6: " << (i+j)/4 << "\t" << i  << "\t" << j << std::endl;
					}
					break;
					default:
					k_vector[i][j] = 0;
					break;
				}
			}		
		}
	}
	return;
}

void FEM1DDD::Dis_kvec() {

	std::cout << "The Stiffness matrix is as follows: " << std::endl;

	std::cout << std::setprecision(2) << std::fixed;

	for (int i=0; i<klen; i++) {

		for (int j=0; j<kwid; j++) {
			
//			k_vector[i][j] = k_vector[i][j]*elasticity/h;
			std::cout << k_vector[i][j] << "\t";
//			k_vector[i][j] = k_vector[i][j]*elasticity/h[0];
		}

		std::cout<< std::endl;
		
	}

	std::cout << std:: endl;
	return;
}

void FEM1DDD::Mat_ludc() {

  	l_vector.resize(klen);
	for (int i = 0; i < klen; ++i){
		l_vector[i].resize(kwid);
	}

  	u_vector.resize(klen);
	for (int i = 0; i < klen; ++i){
		u_vector[i].resize(kwid);
	}   

	for (int i = 0; i<klen; i++) {
		l_vector[i][0] = k_vector[i][0];
//		std::cout<<"L initialised"<<i<<"0"<<std::endl;
	}

	for (int i = 0; i<klen; i++) {
		for (int j = 0; j<kwid; j++) {
			if (i == 0 && j != 0) {
				u_vector[0][j] = k_vector[0][j]/l_vector[0][0];
//				std::cout<<"U initialised"<<"0"<<j<<std::endl;
			}
			if (i == j) {
				u_vector[i][j] = 1;			
//				std::cout<<"U initialised"<<i<<j<<std::endl;
			}
		}
	}

	for (int i = 1; i<klen; i++) {
		for (int j = i; j<klen; j++) {
			double temp = k_vector[j][i];
			for (int k = 0; k<i; k++) {
				temp = temp - (l_vector[j][k]*u_vector[k][i]);
			}
			l_vector[j][i] = temp;
//			std::cout<<"L initialised"<<j<<i<<std::endl;
		}
		for (int j = i+1; j<klen; j++) {
			double temp = k_vector[i][j];
			for (int k = 0; k<i; k++) {
				temp = temp - (l_vector[i][k]*u_vector[k][j]);
			}
			u_vector[i][j] = temp/l_vector[i][i];
//			std::cout<<"U initialised"<<i<<j<<std::endl;
		}
	}
	return;
}

void FEM1DDD::Dis_MtLU() {

	std::cout << "The Lower triangle matrix is as follows: " << std::endl;

	std::cout << std::setprecision(2) << std::fixed;

	for (int i=0; i<klen; i++) {

		for (int j=0; j<kwid; j++) {
			
			std::cout << l_vector[i][j] << "\t";
		}

		std::cout<< std::endl;
		
	}

	std::cout << std:: endl;

	std::cout << "The Upper triangle matrix is as follows: " << std::endl;

	for (int i=0; i<klen; i++) {

		for (int j=0; j<kwid; j++) {
			
			std::cout << u_vector[i][j] << "\t";
		}

		std::cout<< std::endl;
		
	}

	std::cout << std:: endl;
	return;
}

void FEM1DDD::Gen_forf() {

	std::cout << "The force matrix is as follows: " << std::endl;

	for (int i=0; i<numofele-1; i++) {
		f_vector.push_back(forfun*(h[i]+h[i+1])/2);
	}

	for (int i=0; i<numofele-1; i++) {
		std::cout << f_vector[i] << std::endl;;
	}

	std::cout << std::endl;
	return;
}

void FEM1DDD::Gen_gvec() {

	g_vector.resize(numofele-1);
	g_vector[0] = bc1;
	g_vector[numofele-2] = bc2*elasticity/h[numofele-1];
	std::cout << "The G vector is: " << std::endl;
	for (int i=0; i<numofele-1; i++) {
		std::cout << g_vector[i] << std::endl;
	}
	return;
}

void FEM1DDD::Sol_GaSi() {

	int dlen = numofele-1;
	for (int i=0; i<dlen; i++) {
		d_vector.push_back(0);
	}
	double d_conv[dlen];

	int conv_check = 0;

	for (int i=0; conv_check==0; i++) {
		for (int j=0; j<dlen; j++) {

			d_conv[j] = d_vector[j];

			double expansion = (f_vector[j] + g_vector[j])/k_vector[j][j];

			for (int k=0; k<dlen; k++) {
				if (k!=j) {
					expansion = expansion - (k_vector[j][k]*d_vector[k]/k_vector[j][j]);
				}
			}

			d_vector[j] = expansion;
		}

		std::cout << "Iteration: " << i << "\t";

		std::cout << std::setprecision(6) << std::fixed;

		for (int j=0; j<dlen; j++) {
			std::cout << d_vector[j] << "\t";
		}

		std::cout << std::endl;

		double diff = 0;

		for (int j=0; j<dlen; j++) {

			diff = diff + std::abs(std::pow(d_conv[j],2)-std::pow(d_vector[j],2));

		}

		if (std::sqrt(diff)<conv_criteria) conv_check = 1;

	}
	
	void FEM1DDD::back_LU()	{

		std::vector<double> y();
		y[0]=f_vector[0];
		
		//forward 
			for (int j = 0; j <= numofele; j=j+1)	
			{
			t=0;
				for (int num = 0; num < j; num=num+1)
				{
					t = t + ( L[num][j] * y[num] );
					cout << "cal for t" << "\t" << num << endl;
				}
				y[j]= (1/l_vector[j][j])*(f_vector[j]-t);
			}
		
		
		//backward
			for (int j = numofele; j >= 0; j=j-1)
			{
			t=0;
				for (int num = j; num <= numofele; num=num+1)
				{
					t = t + ( U[j][num+1] * d_vector[num+1] );
					cout << "cal for t" << "\t" << num << endl;
				}			
				d_vector[j]=(1/u_vector[j][j]) * (y[j]-t) ;	
			}



	return;
}
