
/****************************************************************************/
/*	Projekt:	Runge Kutta vierter Ordnung / Doppelpendel					*/
/*	file:		03_RK4_Doppelpendel.cpp										*/
/*	Autor:		Milionis Philipp / 1010925									*/
/*	Datum:		2013-04-03													*/
/****************************************************************************/

#include <stdio.h>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

/****************************************************************************/
/*	Parameter																*/
/****************************************************************************/

#define			N_ODEs		4						// number of ODE
const double	EPS			=0.001;					// epsilon in s
const double	TIME_MAX	=60.0;			// max for t ss*mm*hh*dd

const double	G			=9.81;
const double	M			=1.0;
const double	L			=1.0;

const double	INI_Y0		=  M_PI*2/3 ;				// phi 1
const double	INI_Y1		= -M_PI*2/3 ;				// phi 2
const double	INI_Y2		=  0.0 ;				// p 1
const double	INI_Y3		=  0.0 ;				// p 2

/****************************************************************************/
/*	Prototypen																*/
/****************************************************************************/

double F_ODE	(double TIME,	double Y[],		int I);
double runge4	(double TIME,	double Y[],		double EPS);
double ENERGY	(double Y[]);
double position	(double Y[],	double POS1[],	double POS2[]);
double poincare	(double Y[],	double POINVAL[]);


/****************************************************************************/
/*	main																	*/
/****************************************************************************/

int main(void){
	
	cout << endl << endl <<
	" ************************************************************" << endl <<
	" Project:	Computerorientierte Physik 03_RK4_Doppelpendel.cpp" << endl <<
	" Author:	MILIONIS Philipp / 1010925" << endl <<
	" Date:		2013-04-03" << endl <<
	" ************************************************************" << endl;
	
	// data file
	ofstream DATAFILE("03_RK4_Doppelpendel.dat");
	DATAFILE.setf(ios_base::scientific,ios_base::floatfield);
	
	double TIME = 0.0, Y[N_ODEs], NRG = 0.0, POS1[2], POS2[2], POINVAL[4];

	Y[0]=INI_Y0;
	Y[1]=INI_Y1;
	Y[2]=INI_Y2;
	Y[3]=INI_Y3;
	position(Y, POS1, POS2);
	poincare(Y, POINVAL);

	DATAFILE	<< "# t" << "\t" << "phi1"<< "\t"<< "phi2"<<"\t"<<"p1"<<"\t"<<"p2"
				<< "\t" << "energy"
				<< "\t" << "x0" <<"\t"<<"y0" <<"\t"<<"x1" <<"\t"<<"y1" <<"\t"<<"x2" <<"\t"<<"y2"
				<< "\t" << "poin-phi1"  << "\t" << "poin-phi2" << "\t" << "poin-p1"  << "\t" << "poin-p2"  << endl;

	DATAFILE	<< TIME << "\t" << Y[0] << "\t" << Y[1] << "\t" << Y[2] << "\t" << Y[3]
				<< "\t" << NRG
				<< "\t" << 0 << "\t"<< 0 << "\t"<< POS1[0] << "\t"<< POS1[1] << "\t"<< POS2[0] << "\t"<< POS2[1]
				<< "\t" << POINVAL[0]<< "\t"<< POINVAL[1] << "\t"<< POINVAL[2]<< "\t"<< POINVAL[3] << endl;

	for (TIME=0; TIME <= TIME_MAX ; TIME+=EPS) {
		
		// Ausgabe Fortschritt in Sekunden (funktioniert nicht ganz sauber, ist mir aber egal)
		double tmp= (double)(int)(TIME) ;
		if (     TIME - tmp < 0.001) {
		cout  << TIME << "s von " << TIME_MAX << "s" << endl ;
		}
		
		runge4(TIME, Y, EPS);
		position(Y, POS1, POS2);
		NRG = ENERGY(Y);
		poincare(Y, POINVAL);

		DATAFILE	<< TIME << "\t" << Y[0] << "\t" << Y[1] << "\t" << Y[2] << "\t" << Y[3]
			<< "\t" << NRG
			<< "\t" << 0 << "\t"<< 0 << "\t"<< POS1[0] << "\t"<< POS1[1] << "\t"<< POS2[0] << "\t"<< POS2[1]
			<< "\t" << POINVAL[0]<< "\t"<< POINVAL[1] << "\t"<< POINVAL[2]<< "\t"<< POINVAL[3] << endl;
	}
	
	cout << endl << " Datenfile 03_RK4_Doppelpendel.dat geschrieben";

	cout << endl << endl << " fertig 03_RK4_Doppelpendel.cpp " << endl << endl;

	return 0;
}

/****************************************************************************/
/* ODEs																		*/
/****************************************************************************/

double  F_ODE(double TIME, double Y[], int I) {

	if (I==0)
		return(
			   1/(M*L*L)*(Y[2]-Y[3]*cos(Y[0]-Y[1]))/
			   (1+sin(Y[0]-Y[1])*sin(Y[0]-Y[1]))
			   );

	if (I==1)
		return(
			   1/(M*L*L)*(2*Y[3]-Y[2]*cos(Y[0]-Y[1]))/
			   (1+sin(Y[0]-Y[1])*sin(Y[0]-Y[1]))
			   );

	if (I==2)
		return(
			   -1/(M*L*L)*(Y[2]*Y[3]*sin(Y[0]-Y[1]))/
			   (1+sin(Y[0]-Y[1])*sin(Y[0]-Y[1]))+
			   1/(M*L*L)*(Y[2]*Y[2]+2*Y[3]*Y[3]-2*Y[2]*Y[3]*cos(Y[0]-Y[1]))/
			   ((1+sin(Y[0]-Y[1])*sin(Y[0]-Y[1]))*(1+sin(Y[0]-Y[1])*sin(Y[0]-Y[1])))*
			   sin(Y[0]-Y[1])*cos(Y[0]-Y[1])-2*M*G*L*sin(Y[0])
			   );

	if (I==3)
		return(
			   1/(M*L*L)*(Y[2]*Y[3]*sin(Y[0]-Y[1]))/
			   (1+sin(Y[0]-Y[1])*sin(Y[0]-Y[1]))-
			   1/(M*L*L)*(Y[2]*Y[2]+2*Y[3]*Y[3]-2*Y[2]*Y[3]*cos(Y[0]-Y[1]))/
			   ((1+sin(Y[0]-Y[1])*sin(Y[0]-Y[1]))*(1+sin(Y[0]-Y[1])*sin(Y[0]-Y[1])))*
			   sin(Y[0]-Y[1])*cos(Y[0]-Y[1])-M*G*L*sin(Y[1])
			   );
}

/****************************************************************************/
/* RK4																		*/
/****************************************************************************/

double runge4(double TIME, double Y[], double EPS) {
	double h=EPS/2.0;
	double T1[N_ODEs], T2[N_ODEs], T3[N_ODEs];
	double K1[N_ODEs], K2[N_ODEs], K3[N_ODEs],K4[N_ODEs];
	int I;

	for (I=0; I < N_ODEs ; I++)
		T1[I] = Y[I]+0.5*	(K1[I] = EPS * F_ODE ( TIME			, Y  , I ));
	
	for (I=0; I < N_ODEs ; I++)
		T2[I] = Y[I]+0.5*	(K2[I] = EPS * F_ODE ( TIME + h		, T1 , I ));
	
	for (I=0; I < N_ODEs ; I++)
		T3[I] = Y[I]+		(K3[I] = EPS * F_ODE ( TIME + h		, T2 , I ));
	
	for (I=0; I < N_ODEs ; I++)
							(K4[I] = EPS * F_ODE ( TIME + EPS	, T3 , I ));

	for (I=0; I < N_ODEs ; I++)
		Y[I] += ( K1[I] + 2*K2[I] + 2*K3[I] + K4[I] ) / 6.0;
}

/****************************************************************************/
/* Energie																	*/
/****************************************************************************/

double ENERGY ( double Y[]){
	return (Y[2]*Y[2]+2*Y[3]*Y[3]-2*Y[2]*Y[3]*cos(Y[0]-Y[1]))
	/(2*M*L*L)/
	(1+sin(Y[0]-Y[1])*sin(Y[0]-Y[1]))+M*G*L*(3-2*cos(Y[0])-cos(Y[1]));
}

/****************************************************************************/
/* Position der Massen														*/
/****************************************************************************/

double position(double Y[], double POS1[],double POS2[]) {
	POS1[0]= L * sin(Y[0]);
	POS2[0]= POS1[0] + L * sin(Y[1]);
	
	POS1[1]= -L * cos(Y[0]);
	POS2[1]= POS1[1]-L * cos(Y[1]);
}

/****************************************************************************/
/* Poincare																	*/
/****************************************************************************/

double poincare(double Y[], double POINVAL[]) {
	if(abs(Y[1])<0.005 and Y[3] > 0.0){
		
		if	( Y[0] < -M_PI )	Y[0] = Y[0] + 2.0 * M_PI;
		if	( Y[0] >  M_PI )	Y[0] = Y[0] - 2.0 * M_PI;
		if	( Y[1] < -M_PI )	Y[1] = Y[1] + 2.0 * M_PI;
		if	( Y[1] >  M_PI )	Y[1] = Y[1] - 2.0 * M_PI;
		
		POINVAL[0] = Y[0];
		POINVAL[1] = Y[1];
		POINVAL[2] = Y[2];
		POINVAL[3] = Y[3];
		
		//cout << endl << " zapp";
	}
}




