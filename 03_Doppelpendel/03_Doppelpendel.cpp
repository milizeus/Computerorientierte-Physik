
/********************************************************************************************/
/* Aufgabe: Runge Kutta vierter Ordnung / Doppelpendel                                      */
/* file:    03_Doppelpendel.cpp                                                             */
/* Autor:   Milionis Philipp / 1010925                                                      */
/* Datum:   2013-03-18                                                                      */
/********************************************************************************************/

/*Bauen Sie Ihr Runge Kutta Programm so um, dass damit die Hamiltongleichungen für das Doppelpendel gelöst werden können.
 Berechnen Sie nach jedem Schritt die Gesammtenergie, so dass Sie die Richtigkeit und Rechengenauigkeit überprüfen können.
 Machen Sie einige Plots von θ2 als Funktion von θ1 und von p2 versus p1 für verschiedene Startbedingungen.
 In diesem Projekt soll insbesondere auch das chaotische Verhalten des Doppelpendels analysiert werden. 
 Verwenden Sie dazu die beiden in der Vorlesung besprochenen Methoden, Poincare Schnitte und die Stabilitätsanalyse von Bahnen. */

/********************************************************************************************/
// Globale Anweisungen
/********************************************************************************************/

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

/********************************************************************************************/
/* Parameter                                                                               */
/********************************************************************************************/

const double ll = 1 ;
const double mm = 1 ;
const double gg = 9.81 ;
const double eps = 0.001 ;
const int nsteps = 60000;


// −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
// Prototypen
// −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
double energie( double yy [4]);
void fcalc( double yy[4], double ff[4] ) ;




// −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
// Hauptprogramm
// −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

int main ( )
{
    double tt , ee ;
    double yy [4] , ff[4] , yaux [4] , k1 [4] , k2 [4] , k3 [4] , k4 [4] ;
    cout << endl ;
    cout << "03_Doppelpendel.cpp" << endl ;
    cout << endl ;
    cout << endl ;
    ifstream infile ("pendel2.start") ;
    infile >> yy [0] ;
    infile >> yy [1] ;
    infile >> yy [2] ;
    infile >> yy [3] ;
    infile.close ( ) ;

    ofstream outfile ("pendel2.out") ;
    tt = 0.0 ;
    ee = energie(yy) ;
    outfile.width (16) ;
    outfile.precision (7) ;
    outfile.setf(ios_base::scientific,ios_base::floatfield);
    outfile << tt << " " << ee << " " << yy[0] << " " << yy[1] << " " << yy[2] << " " << yy[3] << endl ;
    
    
    for ( int n = 0 ; n < nsteps ; n++)
    {
        fcalc ( yy , k1 ) ;
        yaux [0] = yy [0] + k1 [0] * eps * 0.5 ;
        yaux [1] = yy [1] + k1 [1] * eps * 0.5 ;
        yaux [2] = yy [2] + k1 [2] * eps * 0.5 ;
        yaux [3] = yy [3] + k1 [3] * eps * 0.5 ;
        fcalc ( yaux , k2 ) ;
        yaux [0] = yy [0] + k2 [0] * eps * 0.5 ;
        yaux [1] = yy [1] + k2 [1] * eps * 0.5 ;
        yaux [2] = yy [2] + k2 [2] * eps * 0.5 ;
        yaux [3] = yy [3] + k2 [3] * eps * 0.5 ;
        fcalc ( yaux , k3 ) ;
        yaux [0] = yy [0] + k3 [0] * eps ;
        yaux [1] = yy [1] + k3 [1] * eps ;
        yaux [2] = yy [2] + k3 [2] * eps ;
        yaux [3] = yy [3] + k3 [3] * eps ;
        fcalc ( yaux , k4 ) ;
        yy [0] = yy [0] + ( k1 [0] + 2.0 * k2 [0] + 2.0 * k3 [0] + k4 [0] )* eps / 6.0 ;
        yy [1] = yy [1] + ( k1 [1] + 2.0 * k2 [1] + 2.0 * k3 [1] + k4 [1] )* eps / 6.0 ;
        yy [2] = yy [2] + ( k1 [2] + 2.0 * k2 [2] + 2.0 * k3 [2] + k4 [2] )* eps / 6.0 ;
        yy [3] = yy [3] + ( k1 [3] + 2.0 * k2 [3] + 2.0 * k3 [3] + k4 [3] )* eps / 6.0 ;
        ee = energie ( yy ) ;
        tt = tt + eps ;
        outfile.width ( 16 ) ;
        outfile.precision ( 7 ) ;
        outfile.setf(ios_base::scientific,ios_base::floatfield);
        outfile << tt << " " << ee << " " << yy [0]
        << " " << yy [1] << " " << yy [2] << " " << yy [3] << endl ;
    }
    outfile.close ( ) ;
    cout << endl ;
    cout << " Fertig! " << endl ;
    cout << endl ;
    return 0 ;

}
    


// −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
// Unterprogramme
// −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−


double energie ( double yy [ 4 ] ) {
    double ee ;
    //  ee = ( yy [ 2 ]  yy [ 2 ] + 2yy [ 3 ]  yy [ 3 ] − 2yy [ 2 ]  yy [ 3 ]  cos ( yy [0]−yy [ 1 ] ) ) / ( 2mm l l  l l ( 1 + s i n ( yy [0]−yy [ 1 ] )  s i n ( yy [0]−yy [ 1 ] ) ) )+ mmgg l l ( 3 − 2 cos ( yy [ 0 ] ) − cos ( yy [ 1 ] ) ) ;
    return ee ;
}


// −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
void fcalc ( double yy [4] , double ff [4] ) {
    double fakt , den;
    fakt = 1 / ( mm * ll * ll );
    den = 1 /( 1 + sin( yy[0]−yy[1] ) * sin( yy[0]−yy[1] ) ) ;
    //    ff [0] = fakt *den *( yy [2]−yy [3] * cos ( yy [0]−yy [1] ) ) ;
    //    ff [1] = fakt *den *( 2*yy [3]−yy [2] * cos ( yy [0]−yy [1] ) ) ;
    //    ff [2] = − fakt *den*yy [2] * yy [3] * sin ( yy [0]−yy [1] )+ fakt *den*den *( yy [2] * yy [2] + 2*yy [3] * yy [3] − 2*yy [2] * yy [3] * cos ( yy [0]−yy [1] ) )* sin ( yy [0]−yy [1] ) *cos ( yy [0]−yy [1] ) − 2*mm* ll *gg* sin ( yy [0] ) ;
    //    ff [3] = fakt *den*yy [2] * yy [3] * sin ( yy [0]−yy [1] ) − fakt *den*den *( yy [2] * yy [2]+ 2*yy [3] * yy [3]− 2*yy [2] * yy [3] * cos ( yy [0]−yy [1] ) )* sin ( yy [0]−yy [1] ) *cos ( yy [0]−yy [1] ) − mm* ll *gg* sin ( yy [1] ) ;
}
// −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−


    
    
    