#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int g=50;        //Seitenlänge Gitter
const int a=g*g;        //Anzahl Gitterpunkte
const double epsilon=0.001;

//Prototypen-------------------------------------------------

void feld( int neben[a][4] );
void ab ( double phi[a], bool change[a], double rho[a]);

//Hauptprogramm--------------------------------------------------

int main()
{
    int neben[a][4];
    double phi[a], rho[a], fehler=1.0, delta, h;
    bool change[a];
    
    feld(neben);
    ab (phi, change,rho);
    
    while(fehler > epsilon)
    {
                 fehler=0.0;
                 
                 for( int i=0; i<a; i++)
                 {
                      if(change[i])
                      {
                                   h=0.0;
                                   for( int j=0; j<4; j++) h=h+phi[neben[i][j]];
                                   
                                   h=h/4.0-rho[i];
                                   delta=abs(phi[i]-h);
                                   phi[i]=h;
                                   if( delta>fehler) fehler=delta;
                      }
                 }
    }
    
    ofstream output("Poisson.txt");
    
    for( int x=0; x<g; x++)
    {
         for( int y=0; y<g; y++)
         {
                            
              //output.width(10);
              //output.precision(3);
              //output.setf(ios_base::scientific,ios_base::floatfield);
              output << x << " " << y << " " << phi[x+y*g] << endl;
         }
    }
    
     output.close();
}

//Unterprogramme------------------------------------------------

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

void ab (double phi[a], bool change[a], double rho[a])
{
     int h;
     
     for ( int i=0; i<a; i++)
     {
         phi[i]=0.0;
         rho[i]=0.0;
         change[i]=true;
     }
     
     h=g/2+(g/2+2)*g;
     rho[h]=10.0;
     
     h=g/2+(g/2-2)*g;
     rho[h]=-10.0;
}
