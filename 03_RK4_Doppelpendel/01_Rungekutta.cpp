#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

//Prototypen--------------------------------------------------------------------

double energie (vector<double> y);
int rungekutta (vector<double> a, vector<double> &k1, vector<double> &k2, vector<double> &k3, vector<double> &k4);

//Globale Variablen-------------------------------------------------------------

int d=4;                  //Phasenraumdimension
double epsilon=0.001;      //Schrittweite

const double g=9.81;             
double m=4;                   
double l=2;     
            
//Hauptprogramm-----------------------------------------------------------------

int main()
{
	int x=0, steps=60000;
	double h;
	vector<double> k1(d), k2(d), k3(d), k4(d), y(d);	
		
	cout << endl << " Programm Runge-Kutta Doppelpendel." << endl << endl;
	
	while(x==0)
	{
    	cout << " Bitte geben Sie die Startwerte ein. " << endl << endl;
    	
    	for(int i=0; i<d/2; i++)
    	{
                cout << " phi" << i+1 << "= ";
                cin >> y[i];
        }
        
    	for(int i=d/2; i<d; i++)
    	{
                cout << " p" << i-d/2+1 << "= ";
                cin >> y[i];
        }
                      
        ofstream output("03_Ergebnisse.txt");
            
        for(int i=0; i<steps; i++)
        {                             
                output << i;
                
                h=energie (y);
                output << " " << h;
                
                for(int j=0; j<d; j++) output << " " << y[j];
                output << endl;
                 
                rungekutta (y, k1, k2, k3, k4);
                 
                for(int j=0; j<d; j++) y[j] = y[j] + (k1[j] + 2*k2[j] + 2*k3[j] + k4[j])*epsilon/6;
        }
        
        cout << endl << " Fertig! " << endl;
        cout << " 0...Nochmal!" << endl;
        cout << " 1...Schliessen. ";
        cin >> x;
        
    }
}

//Unterprogramme------------------------------------------------------------------

double energie (vector<double> y)
{
    return (y[2]*y[2]+2*y[3]*y[3]-2*y[2]*y[3]*cos(y[0]-y[1]))/(2*m*l*l)/(1+sin(y[0]-y[1])*sin(y[0]-y[1]))+m*g*l*(4-2*cos(y[0])-cos(y[1]));
}

int rungekutta (vector<double> a, vector<double> &k1, vector<double> &k2, vector<double> &k3, vector<double> &k4)
{
    vector<double> y(d);
    
    for(int j=0; j<d; j++)  y[j]=a[j];   
    k1[0]=1/(m*l*l)*(y[2]-y[3]*cos(y[0]-y[1]))/(1+sin(y[0]-y[1])*sin(y[0]-y[1]));
    k1[1]=1/(m*l*l)*(2*y[3]-y[2]*cos(y[0]-y[1]))/(1+sin(y[0]-y[1])*sin(y[0]-y[1]));
    k1[2]=-1/(m*l*l)*(y[2]*y[3]*sin(y[0]-y[1]))/(1+sin(y[0]-y[1])*sin(y[0]-y[1]))+1/(m*l*l)*(y[2]*y[2]+2*y[3]*y[3]-2*y[2]*y[3]*cos(y[0]-y[1]))/((1+sin(y[0]-y[1])*sin(y[0]-y[1]))*(1+sin(y[0]-y[1])*sin(y[0]-y[1])))*sin(y[0]-y[1])*cos(y[0]-y[1])-2*m*g*l*sin(y[0]);
    k1[3]=1/(m*l*l)*(y[2]*y[3]*sin(y[0]-y[1]))/(1+sin(y[0]-y[1])*sin(y[0]-y[1]))-1/(m*l*l)*(y[2]*y[2]+2*y[3]*y[3]-2*y[2]*y[3]*cos(y[0]-y[1]))/((1+sin(y[0]-y[1])*sin(y[0]-y[1]))*(1+sin(y[0]-y[1])*sin(y[0]-y[1])))*sin(y[0]-y[1])*cos(y[0]-y[1])-m*g*l*sin(y[1]);
       
    for(int j=0; j<d; j++)  y[j]=a[j] + epsilon*k1[j]/2;    
    k2[0]=1/(m*l*l)*(y[2]-y[3]*cos(y[0]-y[1]))/(1+sin(y[0]-y[1])*sin(y[0]-y[1]));
    k2[1]=1/(m*l*l)*(2*y[3]-y[2]*cos(y[0]-y[1]))/(1+sin(y[0]-y[1])*sin(y[0]-y[1]));
    k2[2]=-1/(m*l*l)*(y[2]*y[3]*sin(y[0]-y[1]))/(1+sin(y[0]-y[1])*sin(y[0]-y[1]))+1/(m*l*l)*(y[2]*y[2]+2*y[3]*y[3]-2*y[2]*y[3]*cos(y[0]-y[1]))/((1+sin(y[0]-y[1])*sin(y[0]-y[1]))*(1+sin(y[0]-y[1])*sin(y[0]-y[1])))*sin(y[0]-y[1])*cos(y[0]-y[1])-2*m*g*l*sin(y[0]);
    k2[3]=1/(m*l*l)*(y[2]*y[3]*sin(y[0]-y[1]))/(1+sin(y[0]-y[1])*sin(y[0]-y[1]))-1/(m*l*l)*(y[2]*y[2]+2*y[3]*y[3]-2*y[2]*y[3]*cos(y[0]-y[1]))/((1+sin(y[0]-y[1])*sin(y[0]-y[1]))*(1+sin(y[0]-y[1])*sin(y[0]-y[1])))*sin(y[0]-y[1])*cos(y[0]-y[1])-m*g*l*sin(y[1]);   
       
    for(int j=0; j<d; j++)  y[j]=a[j] + epsilon*k2[j]/2;    
    k3[0]=1/(m*l*l)*(y[2]-y[3]*cos(y[0]-y[1]))/(1+sin(y[0]-y[1])*sin(y[0]-y[1]));
    k3[1]=1/(m*l*l)*(2*y[3]-y[2]*cos(y[0]-y[1]))/(1+sin(y[0]-y[1])*sin(y[0]-y[1]));
    k3[2]=-1/(m*l*l)*(y[2]*y[3]*sin(y[0]-y[1]))/(1+sin(y[0]-y[1])*sin(y[0]-y[1]))+1/(m*l*l)*(y[2]*y[2]+2*y[3]*y[3]-2*y[2]*y[3]*cos(y[0]-y[1]))/((1+sin(y[0]-y[1])*sin(y[0]-y[1]))*(1+sin(y[0]-y[1])*sin(y[0]-y[1])))*sin(y[0]-y[1])*cos(y[0]-y[1])-2*m*g*l*sin(y[0]);
    k3[3]=1/(m*l*l)*(y[2]*y[3]*sin(y[0]-y[1]))/(1+sin(y[0]-y[1])*sin(y[0]-y[1]))-1/(m*l*l)*(y[2]*y[2]+2*y[3]*y[3]-2*y[2]*y[3]*cos(y[0]-y[1]))/((1+sin(y[0]-y[1])*sin(y[0]-y[1]))*(1+sin(y[0]-y[1])*sin(y[0]-y[1])))*sin(y[0]-y[1])*cos(y[0]-y[1])-m*g*l*sin(y[1]);    
       
    for(int j=0; j<d; j++)  y[j]=a[j] + epsilon*k3[j];    
    k4[0]=1/(m*l*l)*(y[2]-y[3]*cos(y[0]-y[1]))/(1+sin(y[0]-y[1])*sin(y[0]-y[1]));
    k4[1]=1/(m*l*l)*(2*y[3]-y[2]*cos(y[0]-y[1]))/(1+sin(y[0]-y[1])*sin(y[0]-y[1]));
    k4[2]=-1/(m*l*l)*(y[2]*y[3]*sin(y[0]-y[1]))/(1+sin(y[0]-y[1])*sin(y[0]-y[1]))+1/(m*l*l)*(y[2]*y[2]+2*y[3]*y[3]-2*y[2]*y[3]*cos(y[0]-y[1]))/((1+sin(y[0]-y[1])*sin(y[0]-y[1]))*(1+sin(y[0]-y[1])*sin(y[0]-y[1])))*sin(y[0]-y[1])*cos(y[0]-y[1])-2*m*g*l*sin(y[0]);
    k4[3]=1/(m*l*l)*(y[2]*y[3]*sin(y[0]-y[1]))/(1+sin(y[0]-y[1])*sin(y[0]-y[1]))-1/(m*l*l)*(y[2]*y[2]+2*y[3]*y[3]-2*y[2]*y[3]*cos(y[0]-y[1]))/((1+sin(y[0]-y[1])*sin(y[0]-y[1]))*(1+sin(y[0]-y[1])*sin(y[0]-y[1])))*sin(y[0]-y[1])*cos(y[0]-y[1])-m*g*l*sin(y[1]);   

}

