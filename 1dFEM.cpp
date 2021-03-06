

// 1D Finite Element Solver

//Include files

#include <iostream>
#include <vector>

using namespace std;

#idndef 1DFEM_H
#define 1DFEM_H

class 1dFEM
{
public:
	
	int main()
	{
		cout << "program start" << endl;

		// Inputs
		/* meshing inputs */
		double  length= 0.1 /*[m]*/, element_size_of_first_element ;
		int no_of_elements=5000 ;
		float scaling_factor=1 ;
		double Force= 1e+11 /*[N(m^-4)]*/;
		/*propetries of the fluid*/
		double E= 1e+11 /*[Pa]*/, A= 1e-4 /*[m^2]*/ ;
		/*boundary conditions*/
		double field_value_at_x0= 0, field_value_at_xn= 0.001 ;/*[m]*/
		double flux_at_n= ;
		double dx;

		vector<double> element_size(no_of_elements);

		dx=length/no_of_elements;

		 cout << dx << endl;

		if (scaling_factor = 1)
		{
			for (int i = 0; i <= no_of_elements-1; ++i)
			{
				element_size[i] = dx;
            	cout << element_size[i] << endl;

			}
		}


		if (scaling_factor != 1)
		{
			element_size[1]= dx;
			for (int i=2; i <= no_of_elements; i = i+1)	{
				element_size[i] = scaling_factor * (element_size[i-1]);

			}
		}



		// defineing the stifness matrix
		/*defining the element_size of the matrix*/
		vector<double> middle(no_of_elements-1);
		vector<double> up(no_of_elements-2), down(no_of_elements-1,0.);

		for (int i=0; i <= no_of_elements-2 ; i=i+1)
		{
			middle[i] = (E*A) * ( (1/element_size[i]) + (1/element_size[i+1]) );
            cout << "middle matrix  " << i << "\t" << middle [i] << endl;
		}

		for (int i=0; i <= no_of_elements-3 ; i=i+1)
		{
			up[i] = -1* ( (E*A)/element_size[i+1] );
			down[i+1]=up[i];
            cout << "down matrix  " << i << "\t" << down[i] << endl;

		}


		//RHS(right hand side) = Force terms + Boundary Conditions
		/*Force terms*/

		vector<double> f(no_of_elements-1);


		for (int i = 0; i < no_of_elements-1 ; i=i+1)
		{
			f[i]=( (Force*(element_size[i]+element_size[i+1])*A) / 2 );
			cout << "force  " << i << "\t" <<  f[i] << endl;
		}

		/*Boundary Conditions*/
		vector<double> Boundary(no_of_elements-1,0.);

		Boundary[0]=field_value_at_x0;

		Boundary[no_of_elements-2]=field_value_at_xn;

		/*total RHS*/
		vector<double> RHS(no_of_elements-1);
		for (int i = 0; i <= no_of_elements-2; i++)
		{
			RHS[i] =f[i] + (E*A*Boundary[i]/element_size[i] );
			cout << "rhs values  "  << i << "\t" << RHS[i] << "      b  "  << Boundary[i] << "       f "  << f[i] << endl;

		}


		//LU Decomposition-- Thomas Algorithm
		vector<double> up_new(no_of_elements-2);
		vector<double> Force_new(no_of_elements-1);

    	//FORWARD DECOMPOSITION
		up_new[0]= up[0]/middle[0];
		Force_new[0] = RHS[0]/middle[0];

		cout << "changed up  " << up_new[0] << endl;

		for (int i=1; i < no_of_elements-2 ; i=i+1)
		{
			up_new[i]=up[i]/((middle[i])-(down[i]*up_new[i-1]));
			cout << "changed up  " << i << "\t" << up_new[i] << endl;
		}
		cout << "CHANGED RHS  "  << Force_new[0] << endl;

		for (int i=1; i < no_of_elements-1 ; i=i+1)
		{
			Force_new[i]= ((RHS[i])-(Force_new[i-1]*down[i])) / ((middle[i])-(down[i]*up_new[i-1]));
			cout << "CHANGED RHS  "  << i << "\t" << Force_new[i] << endl;
		}

    	//BACKWORD SUBSTITUTION.
		vector<double> field_value(no_of_elements-1);

		field_value[no_of_elements-2] = Force_new[no_of_elements-2];

		for (int i=no_of_elements-3; i>=0; i=i-1){
			field_value[i] = Force_new[i]-( up_new[i]*field_value[i+1] );
			cout << "field  " << i << "\t" << field_value[i] << endl;
		}



		for (int i = 0; i <= no_of_elements; ++i)
		{
			cout <<"field" << field_value[i] << endl;
		}

		return 0;
	}

	if (flux_at_n= )
	{
		

	}
};
#endif
