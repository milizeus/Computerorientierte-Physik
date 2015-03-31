
/****************************************************************************************/
/* Aufgabe: Runge Kutta vierter Ordnung                                                 */
/* file:    02_RK4_Oszillator.cpp                                                       */
/* Autor:   Milionis Philipp / 1010925                                                  */
/* Datum:   2015-03-12                                                                  */
/****************************************************************************************/

#include <stdio.h>
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;

/****************************************************************************************/
/* Parameter                                                                            */
/****************************************************************************************/
#define N_ODEs      2																		// number of first order equations
#define INITIAL_Y0  0.0																		// Initial Y[0]
#define INITIAL_Y1  1.0																		// Initial Y[1]
#define omega       0.0																		// resonant frequency
#define DELTA_T     0.1																		// epsilonsize in s
//#define T_MAX       60 //*M_PI																// max for t
int T_MAX;
#define damping     0.2																		// Damping coefficient
#define Omega       1.0																		// driving frequency
#define amp         0.0																		// driving amplitude



/****************************************************************************************/
/* main start                                                                           */
/****************************************************************************************/

int main(void){
    cout << endl << endl <<
    " ******************************************************" << endl <<
    " Project: Computerorientierte Physik /" << " 02_RK4_Oszillator.cpp" << endl <<
    " Author:  MILIONIS Philipp / 1010925" << endl <<
    " Date:    2013-03-12" << endl <<
    " ******************************************************" << endl << endl;
    
	cout << "Please enter an integer value for time in [s]: ";
	cin >> T_MAX;
	
	
    ofstream dataoutput("dat//02_RK4_Oszillator.dat");
    dataoutput.setf(ios_base::scientific,ios_base::floatfield);

    double time = 0.0, y[N_ODEs];
    //int j,k;
    void runge4(double x, double y[], double epsilon);										// Runge-Kutta function
    
    y[0]=INITIAL_Y0;																		// initial Y[0]
    y[1]=INITIAL_Y1;																		// initial Y[1]
    dataoutput  << "# t" << " "  << "x" << " "  << "v" << endl;
    dataoutput  << time << " "  << y[0] << " "  << y[1] << endl;
    
    for (time=0; time <= T_MAX+DELTA_T/5.0 ; time+=DELTA_T) {
        runge4(time, y, DELTA_T);
        dataoutput << time << " " << y[0] << " " << y[1] << endl;
    } // end for
    
    printf(" Datenfile 02_RK4_Oszillator.dat geschrieben");
	
    cout << endl << endl << " fertig 02_RK4_Oszillator.cpp " << endl << endl;
    
	dataoutput.close();
	
	return 0;
}

/****************************************************************************************/
/* ODEs                                                                                 */
/****************************************************************************************/
double  F_ODE(double time, double y[], int i) {
    switch (i) {
        case 0:																				// dY[0]/dt
            return( y[1]);
            break;
            
        case 1:																				// dY[1]/dt
            return( - damping*y[1] - omega*omega*y[1] + amp*sin(Omega * time)-y[0] );
            break;
        default:
            break;
    } // end switch
    return 0;
} // end double



/****************************************************************************************/
/* RK4                                                                                  */
/****************************************************************************************/
void runge4(double time, double y[], double epsilon) {
    double h=epsilon/2.0;																	// the midpoint
    double t1[N_ODEs], t2[N_ODEs], t3[N_ODEs];												// temporary storage arrays
    double k1[N_ODEs], k2[N_ODEs], k3[N_ODEs],k4[N_ODEs];									// for Runge-Kutta
    int i;
    
    for (i=0; i < N_ODEs ; i++) t1[i] = y[i]+   0.5*(k1[i] = epsilon * F_ODE ( time            , y  , i ));
    for (i=0; i < N_ODEs ; i++) t2[i] = y[i]+   0.5*(k2[i] = epsilon * F_ODE ( time + h        , t1 , i ));
    for (i=0; i < N_ODEs ; i++) t3[i] = y[i]+       (k3[i] = epsilon * F_ODE ( time + h        , t2 , i ));
    for (i=0; i < N_ODEs ; i++)                     (k4[i] = epsilon * F_ODE ( time + epsilon  , t3 , i ));
    
    for (i=0; i < N_ODEs ; i++) y[i] += ( k1[i] + 2*k2[i] + 2*k3[i] + k4[i] ) / 6.0;
	return;
}





