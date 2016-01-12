#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace std;

void Euler(double* const q, double* const p, const double& dt, double& H);

int main(){

const double dt = 0.05; 		//step size
double t = 0;
const double tEnd = 20*M_PI;
const double steps = tEnd/dt; 	// number of steps
const double e=0.6; 			//eccentricity of the ellipse
double q[2]; 			// position of the second body
double p[2];
double H;

// initial conditions
q[0]=1-e;
q[1]=0;
p[0]=0;
p[1]= sqrt((1+e)/(1-e));

//Symplectic Euler method

H = 0.5*(p[0]*p[0]+p[1]*p[1])-(1/sqrt(q[0]*q[0]+q[1]*q[1]));

ofstream out("data.txt");   
out << t << "\t" << q[0] << "\t" << q[1] << "\t" << p[0] << "\t" << p[1] << "\t" << H << endl;

for (int i=0; i<steps; i++){
Euler(q, p , dt, H);
t += dt;
out << t << "\t" << q[0] << "\t" << q[1] << "\t" << p[0] << "\t" << p[1] << "\t" << H << endl;
}


out.close();
return(0);
}


void Euler(double* const q, double* const p, const double& dt, double& H){
double x=q[0]*q[0]+q[1]*q[1];

 for(int i=0; i<2; i++){
  p[i] -= dt*q[i]/(pow(x,3.0/2.0));
  q[i] += dt*p[i];
 }

H = 0.5*(pow(p[0],2)+pow(p[1],2))-(1/sqrt(x));
}
