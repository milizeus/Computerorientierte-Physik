#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int g=35;         //Seitenlänge Gitter
const int a=g*g;         //Anzahl Gitterpunkte

ifstream input ("Anfangswerte.txt");
ofstream output ("MonteCarlo.txt");

//Prototypen-----------------------------------------------------------------

void in(double& k0, double& dk, double& maxk, int& jumps, int& euqisteps, int& steps);
void feld(int neben[a][4]);
void hotstart(int spin[a]);
void coldstart(int spin[a]);
void update(int N, int spin[a], int neben[a][4], double k);
void out(int spin[a], int neben[a][4]);


//Hauptprogramm---------------------------------------------------------------

int main()
{
    int neben[a][4];
    int spin[a];
    double k0, dk, maxk;
    int steps, equisteps, jumps, h;
    
    cout << endl << " Program Monte Carlo" << endl << endl;
    cout << " Art der Startkonfiguration (1=coldstart, 2=hotstart)  ";
    cin >> h;
    cout << endl;
    
    in ( k0, dk, maxk, jumps, equisteps, steps );  //Einlesen der Parameter
    feld ( neben );                                //Initialisieren Nachbarnfeld
       
    if (h==1) coldstart(spin);                     //Startkonfiguration
    else hotstart(spin);                                                     
   
    cout << " Verschiedene Kopplungswerte werden untersucht:" << endl;
     
    for ( double k=k0; k<maxk; k=k+dk )          //Läuft über verschiedene Kopplungswerte von k0 bis maxk in dk Schritten                                       
    {
        update ( equisteps, spin, neben, k );     //Equilibriert das Programm mit passend vielen Schritten quisteps bis es im Gleichgewicht ist
        
        cout << " " << k << endl; 
        
        for ( int j=0; j<jumps; j++ )             //Spring durch den Raum aller Konfigurationen                                              
        {
            update ( steps, spin, neben, k );     //Jeder Sprung hat #steps
            out ( spin, neben );                  //Ausgabe 
        }
    }
    
    cout << endl << "Fertig!" << endl;
    
    return 0;      
}

//Unterprograme------------------------------------------------------------------

void in(double& k0, double& dk, double& maxk, int& jumps, int& equisteps, int& steps)
{
     string label;
     int G;
     
     input >> G >> label;
     input >> k0 >> label;
     input >> dk >> label;
     input >> maxk >> label;
     input >> jumps >> label;
     input >> equisteps >> label;
     input >> steps >> label;
     
     cout << " Parameter:" << endl << endl;
     cout << " Seitenlaenge Gitter...............: " << g << endl;
     cout << " Anfangswert der Kopplung.........: " << k0 << endl;
     cout << " Schrittweite der Kopplung........: " << dk << endl;
     cout << " Maximalwert der Kopplung.........: " << maxk << endl;
     cout << " Anzahl der Sprünge druch den" << endl;
     cout << " Raum der Konfigurationen.........: " << jumps << endl;
     cout << " Anzahl der Schritte zum" << endl;
     cout << " Equilibrieren....................: " << equisteps << endl;
     cout << " Anzahl der Schritte pro Sprung...: " << steps << endl << endl;
     
     output << "g.................: " << g << endl;
     output << "k0................: " << k0 << endl;
     output << "dk................: " << dk << endl;
     output << "maxk..............: " << maxk << endl;
     output << "jumps.............: " << jumps << endl;
     output << "equisteps.........: " << equisteps << endl;
     output << "steps.............: " << steps << endl;
}

void feld (int neben[a][4])
{
     int xp, xm, yp, ym, h;
     
     for( int x=0; x<g; x++)
     {
          xp=x+1;
          xm=x-1;
          if(xp==g) xp=0;
          if(xm==-1) xm=g-1;
          
          for ( int y=0; y<g; y++)
          {
              yp=y+1;
              ym=y-1;
              if(yp==g) yp=0;
              if(ym==-1) ym=g-1;
              
              h=x+y*g;
              
              neben[h][0]=x+yp*g;
              neben[h][1]=xp+y*g;
              neben[h][2]=x+ym*g;
              neben[h][3]=xm+y*g;
          }
     }
}

void hotstart (int spin[a])
{
     double h;
     
     for ( int i=0; i<a; i++)
     {
         h=rand()/(RAND_MAX+1.0);
         if( h < 0.5 ) spin[i]=1;
         else spin[i]=-1; 
     }     
}

void coldstart (int spin[a])
{
     for (int i=0; i<a; i++) spin[i]=1;
}

void update ( int N, int spin[a], int neben[a][4], double k)
{
     double rho, z; 
     int h; 
     
     for (int i=0; i<=N; i++)
     {
         for (int j=0; j<a; j++)
         {
             h=0;
             for (int l=0; l<4; l++) h=h+spin[neben[j][l]];
             rho=exp(-2*k*h*spin[j]);
             z=rand()/(RAND_MAX+1.0);
             if (z<rho) spin[j]=-spin[j];
         }
     }
}

void out ( int spin[a], int neben[a][4] )
{
     int mag=0, h=2*a;
     
     for (int i=0; i<a; i++)
     {
         mag=mag+spin[i];
         h=h-spin[i]*(spin[neben[i][0]]+spin[neben[i][1]]);
     }
     
     output << mag << " " << h << endl;
}

     
