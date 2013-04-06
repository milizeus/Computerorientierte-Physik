
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

/****************************************************************/
/*	Parameter													*/
/****************************************************************/

#define			N_ODEs		4						// number of ODE
const double	DELTA_T		=0.001;					// epsilonsize in s
const double	T_MAX		=60.0;					// max for t

const double	g			=9.81;
const double	m			=1.0;
const double	l			=1.0;

const double	INITIAL_Y0	=3*M_PI/4 +0.0 ;		// phi 1
const double	INITIAL_Y1	=4*M_PI/4 +0.0;			// phi 2
const double	INITIAL_Y2	=0.0;					// p 1
const double	INITIAL_Y3	=0.0;					// p 2

/****************************************************************/
/*	Prototypen													*/
/****************************************************************/

double F_ODE	(double time,	double y[],		int i);
double runge4	(double x,		double y[],		double epsilon);
double energie	( double y[]);
double position	(double y[],	double pos1[],	double pos2[]);
double poincare	(double y[],	double poinval[]);


/****************************************************************/
/*	main														*/
/****************************************************************/

int main(void){
	
	cout << endl << endl <<
	" *******************************************************************" << endl <<
	" Project: Computerorientierte Physik /" << " 03_RK4_Doppelpendel.cpp" << endl <<
	" Author:  MILIONIS Philipp / 1010925" << endl <<
	" Date:    2013-03-12" << endl <<
	" *******************************************************************" << endl << endl;
	
	ofstream dataoutput("03_RK4_Doppelpendel.dat");
	dataoutput.setf(ios_base::scientific,ios_base::floatfield);
	double time = 0.0, y[N_ODEs], NRG = 0.0, pos1[2], pos2[2], poinval[4];


	y[0]=INITIAL_Y0;			// initial Y[0]
	y[1]=INITIAL_Y1;			// initial Y[1]
	y[2]=INITIAL_Y2;			// initial Y[1]
	y[3]=INITIAL_Y3;			// initial Y[1]
	position(y, pos1, pos2);
	poincare(y, poinval);


	dataoutput << "# t" << "\t" << "phi1"<< "\t"<< "phi2"<<"\t"<<"p1"<<"\t"<<"p2"
	<< "\t" << "energy"
	<< "\t" << "x0" <<"\t"<<"y0" <<"\t"<<"x1" <<"\t"<<"y1" <<"\t"<<"x2" <<"\t"<<"y2"
	<< "\t" << "poin-phi1"  << "\t" << "poin-phi2" << "\t" << "poin-p1"  << "\t" << "poin-p2"  << endl;

	dataoutput << time << "\t"  << y[0] << "\t"  << y[1] << "\t"  << y[2] << "\t"  << y[3];
	dataoutput << "\t" << NRG;
	dataoutput << "\t" << 0 << "\t"<< 0 << "\t"<< pos1[0] << "\t"<< pos1[1] << "\t"<< pos2[0] << "\t"<< pos2[1];
	dataoutput << "\t" << poinval[0]<< "\t"<< poinval[1] << "\t"<< poinval[2]<< "\t"<< poinval[3] << endl;


	for (time=0; time <= T_MAX ; time+=DELTA_T) {
		runge4(time, y, DELTA_T);
		position(y, pos1, pos2);
		NRG = energie(y);
		poincare(y, poinval);

		dataoutput << time << "\t"  << y[0] << "\t"  << y[1] << "\t"  << y[2] << "\t"  << y[3];
		dataoutput << "\t" << NRG;
		dataoutput << "\t" << 0 << "\t"<< 0 << "\t"<< pos1[0] << "\t"<< pos1[1] << "\t"<< pos2[0] << "\t"<< pos2[1];
		dataoutput << "\t" << poinval[0]<< "\t"<< poinval[1] << "\t"<< poinval[2]<< "\t"<< poinval[3] << endl;
	}

	printf(" Datenfile 03_RK4_Doppelpendel.dat geschrieben");

	cout << endl << endl << " fertig 03_RK4_Doppelpendel.cpp " << endl << endl;

	return 0;
}



/**************************************************************/
/* ODEs                                                       */
/**************************************************************/
double  F_ODE(double time, double y[], int i) {

	if (i==0) return(
					 1/(m*l*l)*(y[2]-y[3]*cos(y[0]-y[1]))/(1+sin(y[0]-y[1])*sin(y[0]-y[1]))
					 );

	if (i==1) return(
					 1/(m*l*l)*(2*y[3]-y[2]*cos(y[0]-y[1]))/(1+sin(y[0]-y[1])*sin(y[0]-y[1]))
					 );

	if (i==2) return(
					 -1/(m*l*l)*(y[2]*y[3]*sin(y[0]-y[1]))/(1+sin(y[0]-y[1])*sin(y[0]-y[1]))+1/(m*l*l)*(y[2]*y[2]+2*y[3]*y[3]-2*y[2]*y[3]*cos(y[0]-y[1]))/((1+sin(y[0]-y[1])*sin(y[0]-y[1]))*(1+sin(y[0]-y[1])*sin(y[0]-y[1])))*sin(y[0]-y[1])*cos(y[0]-y[1])-2*m*g*l*sin(y[0])
					 );

	if (i==3) return(
					 1/(m*l*l)*(y[2]*y[3]*sin(y[0]-y[1]))/(1+sin(y[0]-y[1])*sin(y[0]-y[1]))-1/(m*l*l)*(y[2]*y[2]+2*y[3]*y[3]-2*y[2]*y[3]*cos(y[0]-y[1]))/((1+sin(y[0]-y[1])*sin(y[0]-y[1]))*(1+sin(y[0]-y[1])*sin(y[0]-y[1])))*sin(y[0]-y[1])*cos(y[0]-y[1])-m*g*l*sin(y[1])
					 );
}


/**************************************************************/
/* RK4															*/
/**************************************************************/

double runge4(double time, double y[], double epsilon) {
	double h=epsilon/2.0;
	double t1[N_ODEs], t2[N_ODEs], t3[N_ODEs];                         // tmp
	double k1[N_ODEs], k2[N_ODEs], k3[N_ODEs],k4[N_ODEs];
	int i;

	for (i=0; i < N_ODEs ; i++) t1[i] = y[i]+   0.5*(k1[i] = epsilon * F_ODE ( time            , y  , i ));
	for (i=0; i < N_ODEs ; i++) t2[i] = y[i]+   0.5*(k2[i] = epsilon * F_ODE ( time + h        , t1 , i ));
	for (i=0; i < N_ODEs ; i++) t3[i] = y[i]+       (k3[i] = epsilon * F_ODE ( time + h        , t2 , i ));
	for (i=0; i < N_ODEs ; i++)                     (k4[i] = epsilon * F_ODE ( time + epsilon  , t3 , i ));

	for (i=0; i < N_ODEs ; i++) y[i] += ( k1[i] + 2*k2[i] + 2*k3[i] + k4[i] ) / 6.0;
	
}


/**************************************************************/
/* Energie														*/
/**************************************************************/

double energie ( double y[])
{
    return (y[2]*y[2]+2*y[3]*y[3]-2*y[2]*y[3]*cos(y[0]-y[1]))/(2*m*l*l)/(1+sin(y[0]-y[1])*sin(y[0]-y[1]))+m*g*l*(3-2*cos(y[0])-cos(y[1]));
}


/**************************************************************/
/* Position der Massen											*/
/**************************************************************/

double position(double y[], double pos1[],double pos2[]) {

	pos1[0]= l * sin(y[0]);
	pos2[0]= pos1[0] + l * sin(y[1]);
	
	pos1[1]= -l * cos(y[0]);
	pos2[1]= pos1[1]-l * cos(y[1]);
}

/**************************************************************/
/* Poincare														*/
/**************************************************************/

double poincare(double y[], double poinval[]) {
	
	if(abs(y[1])<0.01 && y[3] > 0.0){
		
		poinval[0]=y[0];
		poinval[1]=y[1];
		poinval[2]=y[2];
		poinval[3]=y[3];
	}
}

