#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

ifstream input ("MonteCarlo.txt");
ofstream output ("Analyse.txt");

//Hauptprogramm---------------------------------------------------------------

int main()
{
    double k0, dk, maxk;
    int steps, equisteps, jumps, g;
    
    string label;
    double mav, mf, hav, hf;
    
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
        
        mav=0.0;
        hav=0.0;
        
        for (int i=0; i<jumps; i++)
        {
            mav=mav+mag[i];
            hav=hav+h[i];
        }
        
        mav=mav/jumps;
        hav=hav/jumps;
        
        mf=0.0;
        hf=0.0;
        
        for (int i=0; i<jumps; i++)
        {
            mf=mf+pow(mav-mag[i], 2);
            hf=hf+pow(hav-h[i], 2);
        }
        
        mf=sqrt(mf)/jumps;
        hf=sqrt(hf)/jumps;
        
        output << k << " " << mav/a << " " << mf/a << " " << hav/a << " " << hf/a << endl;
    }
        
    cout << " Pruefwerte: " << mag[jumps-1] << " " << h[jumps-1] << endl;
    cin >> k0;
        
    return 0;
}
        
        
  
