
/********************************************************************************************/
/* Aufgabe: Runge Kutta vierter Ordnung / Doppelpendel                                      */
/* file:    03_RK4_Doppelpendel.cpp                                                             */
/* Autor:   Milionis Philipp / 1010925                                                      */
/* Datum:   2013-04-03                                                                      */
/********************************************************************************************/

#include <stdio.h>
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;

/********************************************************************************************/
/* Parameter                                                                               */
/********************************************************************************************/
#define N_ODEs      4                                       // number of first order equations
#define DELTA_T     0.001                                     // epsilonsize in s
#define T_MAX       2                               // max for t

const double    g=9.81;
double          m=1;
double          l=1;

#define INITIAL_Y0  1.0                                     // phi 1
#define INITIAL_Y1  1.0                                     // phi 2
#define INITIAL_Y2  0.0                                     // p 1
#define INITIAL_Y3  0.0                                     // p 2

/********************************************************************************************/
/* ODEs                                                                                     */
/********************************************************************************************/
double  F_ODE(double time, double y[], int i) {
    
    if (i==0) return(
                     1/(m*l*l)*(y[2]-y[3]*cos(y[0]-y[1]))/(1+sin(y[0]-y[1])*sin(y[0]-y[1]))
                    );
    
    if (i==1) return(
                     1/(m*l*l)*(2*y[3]-y[2]*cos(y[0]-y[1]))/(1+sin(y[0]-y[1])*sin(y[0]-y[1]))
                     );
    
    if (i==2) return(
                     1/(m*l*l)*(y[2]*y[3]*sin(y[0]-y[1]))/(1+sin(y[0]-y[1])*sin(y[0]-y[1]))+1/(m*l*l)*(y[2]*y[2]+2*y[3]*y[3]-2*y[2]*y[3]*cos(y[0]-y[1]))/((1+sin(y[0]-y[1])*sin(y[0]-y[1]))*(1+sin(y[0]-y[1])*sin(y[0]-y[1])))*sin(y[0]-y[1])*cos(y[0]-y[1])-2*m*g*l*sin(y[0])
                     );
    
    if (i==3) return(
                     1/(m*l*l)*(y[2]*y[3]*sin(y[0]-y[1]))/(1+sin(y[0]-y[1])*sin(y[0]-y[1]))-1/(m*l*l)*(y[2]*y[2]+2*y[3]*y[3]-2*y[2]*y[3]*cos(y[0]-y[1]))/((1+sin(y[0]-y[1])*sin(y[0]-y[1]))*(1+sin(y[0]-y[1])*sin(y[0]-y[1])))*sin(y[0]-y[1])*cos(y[0]-y[1])-m*g*l*sin(y[1])
                     );
}

/********************************************************************************************/
/* main start                                                                               */
/********************************************************************************************/

int main(void){
    cout << endl << endl <<
    " ******************************************************" << endl <<
    " Project: Computerorientierte Physik /" << " 03_RK4_Doppelpendel.cpp" << endl <<
    " Author:  MILIONIS Philipp / 1010925" << endl <<
    " Date:    2013-03-12" << endl <<
    " ******************************************************" << endl << endl;
    
    ofstream dataoutput("03_RK4_Doppelpendel.dat");
   // dataoutput.setf(ios_base::scientific,ios_base::floatfield);
    
    double time = 0.0, y[N_ODEs];
    //int j,k;
    void runge4(double x, double y[], double epsilon);              // Runge-Kutta function
    
    y[0]=INITIAL_Y0;                                                // initial Y[0]
    y[1]=INITIAL_Y1;                                                // initial Y[1]
    y[2]=INITIAL_Y2;                                                // initial Y[1]
    y[3]=INITIAL_Y3;                                                // initial Y[1]
    dataoutput  << "# t" << " "  << "x" << " "  << "v" << endl;
    dataoutput  << time << " "  << y[0] << " "  << y[1] << " "  << y[2] << " "  << y[3] << endl;
    
    
    for (time=0; time <= T_MAX+DELTA_T/5.0 ; time+=DELTA_T) {
        runge4(time, y, DELTA_T);
        dataoutput  << time << " "  << y[0] << " "  << y[1] << " "  << y[2] << " "  << y[3] << endl;
    }
    
    printf(" Datenfile 03_RK4_Doppelpendel.dat geschrieben");
    
    cout << endl << endl << " fertig 03_RK4_Doppelpendel.cpp " << endl << endl;
    
    return 0;
}

/********************************************************************************************/
/* RK4                                                                                      */
/********************************************************************************************/
void runge4(double time, double y[], double epsilon) {
    double h=epsilon/2.0;                                              // the midpoint
    double t1[N_ODEs], t2[N_ODEs], t3[N_ODEs];                         // temporary storage arrays
    double k1[N_ODEs], k2[N_ODEs], k3[N_ODEs],k4[N_ODEs];              // for Runge-Kutta
    int i;
    
    for (i=0; i < N_ODEs ; i++) t1[i] = y[i]+   0.5*(k1[i] = epsilon * F_ODE ( time            , y  , i ));
    for (i=0; i < N_ODEs ; i++) t2[i] = y[i]+   0.5*(k2[i] = epsilon * F_ODE ( time + h        , t1 , i ));
    for (i=0; i < N_ODEs ; i++) t3[i] = y[i]+       (k3[i] = epsilon * F_ODE ( time + h        , t2 , i ));
    for (i=0; i < N_ODEs ; i++)                     (k4[i] = epsilon * F_ODE ( time + epsilon  , t3 , i ));
    
    for (i=0; i < N_ODEs ; i++) y[i] += ( k1[i] + 2*k2[i] + 2*k3[i] + k4[i] ) / 6.0;
	return;
}






