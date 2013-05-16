#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

ifstream input ("MonteCarlo.txt");
ofstream output ("Analyse2.txt");

//Hauptprogramm---------------------------------------------------------------

int main()
{
    double k0, dk, maxk;
    int steps, equisteps, jumps, g;
    
    string label;
    double mav1, hav1, mav2, hav2, c, x;
    
    input >> label >> g;
    input >> label >> k0;
    input >> label >> dk;
    input >> label >> maxk;
    input >> label >> jumps;
    input >> label >> equisteps;
    input >> label >> steps;
    
    cout << " Parameter:" << endl << endl;
    cout << " Seitenlaenge Gitter..............: " << g << endl;
    cout << " Anfangswert der Kopplung.........: " << k0 << endl;
    cout << " Schrittweite der Kopplung........: " << dk << endl;
    cout << " Maximalwert der Kopplung.........: " << maxk << endl;
    cout << " Anzahl der Spruenge druch den" << endl;
    cout << " Raum der Konfigurationen.........: " << jumps << endl;
    cout << " Anzahl der Schritte zum" << endl;
    cout << " Equilibrieren....................: " << equisteps << endl;
    cout << " Anzahl der Schritte pro Sprung...: " << steps << endl << endl;
    
    int a=g*g;
    int mag[jumps], h[jumps];
    
    cout << " Einlesen der Daten..." << endl;
    
    for (double k=k0; k<maxk; k=k+dk)
    {
        for (int i=0; i<jumps; i++)
        {
            input >> mag[i];
            mag[i]=abs(mag[i]);
            input >> h[i];
        }
        
        mav1=0.0;
        hav1=0.0;
        mav2=0.0;
        hav2=0.0;
        c=0.0;
        x=0.0;
        
        for (int i=0; i<jumps; i++)
        {
            mav1=mav1+mag[i];
            hav1=hav1+h[i];
        }
        
        hav1=hav1/jumps;
        mav1=mav1/jumps;
        
        for (int i=0; i<jumps; i++)
        {
            x=x+(mag[i]-mav1)*(mag[i]-mav1);
            c=c+(h[i]-hav1)*(h[i]-hav1);
        }
        
        c=c/jumps;
        x=x/jumps;
        
        
        
        output << k << " " << c/a << " " << x/a << endl;
    }
        
    cout << " Pruefwerte: " << mag[jumps-1] << " " << h[jumps-1] << endl;
    cin >> k0;
        
    return 0;
}
