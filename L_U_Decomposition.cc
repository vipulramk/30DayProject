//L-D Transformation

#include <iostream>

using namespace std;

int main ()

 {


	// an array with 5 rows and 2 columns.
    int 	n=4;
   double a[n][n] = { 	{2,-1,0,0},
   						{-1,2,-1,0},
   						{0,-1,2,-1},
   						{0,0,-1,2}
   					};

	double 	l[n][n]={0};
	double 	u[n][n]={0};
	
	/*for (int i = 0; i < n; i=i+1)
	{
 		u[i][i]=1;
	}*/
	
	
 	// algorithm for L U Decomposition

	/*for (int i = 0; i < n; i=i+1)
		{
			l[i][0]=a[i][0];
			u[0][i]=a[0][i]/l[0][0];
		}*/

 	for (int i = 0; i < n; i=i+1)
 	{
 		
 		for (int j = 0; j < n; j=j+1)
 		{
 //			cout << "i=" << i << "	" << " j=" << j << endl;
 			/*lower triangular matrix*/

 			if (i>=j)
 			{
 				double temp_l=0;

 					for (int ii = 0; ii < j-1; ii=ii+1)
 					{
 						// for (int jj = 0; jj < i-1; jj=jj+1)
 						// {
 						// temp_l= temp_l + ( l[i][jj]*u[ii][j] );

 						// }
 						temp_l= temp_l + ( l[i][ii]*u[ii][j] );
 					}
 					cout << temp_l << "lower" << endl;
 //					cout <<"	"<< endl;

 					cout << "At L: " << i << " " << j << endl << endl;


 				l[i][j]= a[i][j]- temp_l;
 					if (i==0 && j==1)
 					{
 						cout << "error found" << endl;
 						return 0;
 					}

 					l[0][1]=0;

 			}


 			// /*upper triangular matrix*/

 			if (i==j)
 			{
 				u[i][j]=1;
 				cout  << "upper for i=j" << endl;
 //				cout <<"	"<< endl;
 					cout << "At U: " << i << " " << j << endl << endl;

 			}

 			if ( i<j)
 			 {

 				double temp_u=0;

 					for (int iii = 0; iii < i-1; iii=iii+1)
 					{
 						// for (int jjj = 0; jjj < i-1; jjj=jjj+1)
 						// {
 						// temp_u = temp_u+ ( l[i][jjj]*u[iii][j] );
							temp_u = temp_u+ ( l[i][iii]*u[iii][j] );
 						// }
 					}
 					cout << temp_u << "upper for i<j" << endl;
 //					cout <<"	"<< endl;

 				u[i][j]= (1/l[i][i])*(a[i][j] - temp_u );
 					cout << "At U: " << i << " " << j << endl << endl;

 			}


 		}
 	}


 	// displaying the 'l' and 'u' components

 	cout << "        lower       " << endl;
 	// for lower trianfular matrix

	for(int i=0; i<n; i++)    //This loops on the rows.
	{
		for(int j=0; j<n; j++) //This loops on the columns
		{
			cout << l[i][j]  << "  ";
		}
		cout << endl;
	}

	cout << "               " << endl;

	cout << "        upper       " << endl;
	// for upper triangular matrix

	for(int i=0; i<n; i++)    //This loops on the rows.
	{
		for(int j=0; j<n; j++) //This loops on the columns
		{
			cout << u[i][j]  << "  ";
		}
		cout << endl;
	}

	cout << l[0][n-1] << "extra" << endl;
	cout << l[0][n-2] << "extra" << endl;
	cout << l[0][n-3] << "extra" << endl;


  return 0;

 }
