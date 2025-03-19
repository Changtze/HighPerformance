#include "particle.h"
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;

Particle::Particle():
	posVector({0.0, 0.0, 0.0}),
	velVector({0.0, 0.0, 0.0}),
	forceVector({0.0, 0.0, 0.0}),
	type(0),
	mass(1) {}
	
	
Particle::Particle(const std::array<double, 3>& position,
				   const std::array<double, 3>& velocity,
				   int particleType):
	posVector(position),
	velVector(velocity),
	forceVector({0.0, 0.0, 0.0}),
	type(particleType) {
		switch (particleType) {
			case 1:
				mass = 10;
				break;
			default:
				mass = 1;
				break;
		}
	}
	
void Particle::displayVel() const {
	cout << fixed << setprecision(8);
	cout << "Velocity: ";
	for (int dim = 0; dim < 3; ++dim){
		cout << this->velVector[dim] << " ";
	}
	cout << endl;
}

void Particle::displayForce() const{
	cout << fixed << setprecision(8);
	cout << "Force: ";
	for (int dim = 0; dim < 3; ++dim){
		cout << this->forceVector[dim] << " ";
	}
	cout << endl;
}
		
void Particle::displayPos() const{
	cout << fixed << setprecision(8);
	cout << "Position: ";
	for (int dim = 0; dim < 3; ++dim){
		cout << this->posVector[dim] << " ";
	}
	cout << endl;
}
double Particle::getMass() const{
	return this->mass;
}

int Particle::getType() const {
	return this->type;
}

const std::array<double, 3>& Particle::getPosVector() const{
	return this->posVector;
}

const std::array<double, 3>& Particle::getVelVector() const {
	return this->velVector;
}

const std::array<double, 3>& Particle::getForceVector() const{
	return this->forceVector;
}

void Particle::setPosition(const std::array<double, 3>& newPos) {
	this->posVector = newPos;
}

void Particle::setVelocity(const std::array<double, 3>& newVel) {
	this->velVector = newVel;
}

void Particle::setForce(const std::array<double, 3>& newForce) {
	this->forceVector = newForce;
}

void Particle::zeroForce() {
	this->forceVector = {0.0, 0.0, 0.0};
}

void Particle::applyBC(const std::array<double, 3>& domain) {
	// for each dimension, check BCs
	for (int i = 0; i < 3; ++i){
		if (this->posVector[i] < 0.0) {
			
			this->posVector[i] = -this->posVector[i];
			this->velVector[i] = std::abs(velVector[i]);
		}
		else if (posVector[i] > domain[i]){
			
			this->posVector[i] = 2.0 * domain[i] - this->posVector[i];
			this->velVector[i] = -std::abs(this->velVector[i]);
		}
	}
}

