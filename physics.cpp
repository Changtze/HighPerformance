#include "physics.h"
#include <array>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <ostream>
#include <cmath>

using namespace std;

class Particle;  // forward declaration of particle class

Physics::Physics(){
	initialise_LJ();
}

void Physics::initialise_LJ() {
	// Epsilon array
	epsilon[0][0] = 3.0;
	epsilon[0][1] = 15.0;
	epsilon[1][0] = 15.0;
	epsilon[1][1] = 60.0;
	
	// Sigma array
	sigma[0][0] = 1.0;
	sigma[0][1] = 2.0;
	sigma[1][0] = 2.0;
	sigma[1][1] = 3.0;
}

double Physics::getEpsilon(int type1, int type2) const {
	return epsilon[type1][type2];
}	

double Physics::getSigma(int type1, int type2) const{
	return sigma[type1][type2];
}

double Physics::computeLJDeriv(int type1, int type2,
							   double distance, array<double, 3>& r, 
							   int dimension) const {
	
	double sig = sigma[type1][type2];	
	double eps = epsilon[type1][type2];

//	double sig12 = pow(sig, 12);
//	double sig6 = pow(sig, 6);
//	double r14 = pow(sqrt(distance), 14);
//	double r8 = pow(sqrt(distance), 8);
	
	// WRITE ANOTHER TABLE TO PRE-COMPUTE SIGMA
	double sig2 = sig*sig;
	double sig4 = sig2*sig2;
	double sig6 = sig2 * sig4;
	double sig8 = sig4*sig4;
	double sig12 = sig8 * sig4;
	
	
	double r4 = distance * distance;
	double r8 = r4 * r4;
	double r14 = r8 * r4 * distance;
	
	double x_ij = r[0];
	double y_ij = r[1];
	double z_ij = r[2];
	
	// pre-calculating sigma and distance ratios (part of optimisation)
	double sig12r14 = sig12 / r14;	
	double sig6r8 = sig6 / r8;
	
	double phi = 0.0;  // initialise LJ derivative potential
	switch (dimension){
		case 0:
			phi = -24.0 * eps * x_ij * ( 2.0 * sig12r14 - sig6r8 );
			break;
		case 1:
			phi = -24.0 * eps * y_ij * ( 2.0 * sig12r14 - sig6r8 );
			break;
		case 2:
			phi = -24.0 * eps * z_ij * ( 2.0 * sig12r14 - sig6r8 );
			break;
	}   
	return phi;	
}


array<double, 3> Physics::calculateForceVector(int type1, int type2, 
											   array<double, 3> dist_vector, 
											   double r_sq) const {
	double Fx = computeLJDeriv(type1, type2, r_sq, dist_vector, 0);  
	double Fy = computeLJDeriv(type1, type2, r_sq, dist_vector, 1);
	double Fz = computeLJDeriv(type1, type2, r_sq, dist_vector, 2);
	
	array<double, 3> F{Fx, Fy, Fz};
	
	return F;
}
													
double Physics::calculateTemp(double KE, int N) {
	return (2.0/(3.0*kB)) * (KE/N);
}
	
double Physics::calculateVelScale(double currentTemp, double userTemp) {
	return sqrt(userTemp/currentTemp);
}
	

					
																					