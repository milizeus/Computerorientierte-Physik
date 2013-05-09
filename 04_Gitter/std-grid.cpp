/*
 * Author: Ludwig Treyer
 * Program Name: gitter.cpp
 * Compilier Notes:
 * 		use Makefile to compilie or clean
 * 		$ make 
 * 		$ make clean
 *
 * Program Description:
 *
 *
 *************************/

//# include <blitz/array.h>
//#include <stdio.h>
//#include <cmath>
#include <iostream>
//BZ_USING_NAMESPACE(blitz)
using namespace std;
/*
 * Parameter
 **************/ 
#define dim 2	// dimension of grid
#define N1 5	// valuation of dim, now dim 2
#define N2 5
#define dir 2	// pos., neg. direction => 2 

int main()
{
	int n1=0,n1p=0,n1m=0, n2=0,n2p=0,n2m=0, k=0,mu=0;
	mu=dir*dim;
	int neib[N1*N2][mu];
/*
	//init neib
	for(int i=0;i<(N1*N2);i++)
		for(int j=0;j<mu;j++)
			neib[i][j]=0;
	for(int i=0;i<(N1*N2);i++)
	{	cout << i << "- ";	
		for(int j=0;j<mu;j++)
			cout << neib[i][j] << " | ";
		cout << endl;
	}
*/
	for(n1=0; n1<=N1-1; n1++){
		n1p = n1 + 1;
		n1m = n1 - 1;
		if (n1p==N1)
			n1p = 0;
		if (n1m==-1)
			n1m = N1-1;

		for(n2=0; n2<=N2-1; n2++){
			n2p = n2 + 1;
			n2m = n2 - 1;
			if (n2p==N2)
				n2p = 0;
			if (n2m==-1)
				n2m = N2-1;
		
			k = n1 + n2*N1;
		
			neib[k][0] = n1 + n2p*N1;
			neib[k][1] = n1p + n2*N1;
			neib[k][2] = n1 + n2m*N1;
			neib[k][3] = n1m + n2*N1;
		}	
	}
	
	for(int i=0; i<=k; i++)
	{	
		cout << i <<"- " << neib[i][0] << "|" << neib[i][1] << "|" 
						 << neib[i][2] << "|" << neib[i][3] << endl;
	}	


	return 0;
}

