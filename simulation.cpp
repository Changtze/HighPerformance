#include "simulation.h"
#include <random>
#include <iomanip>
#include <iostream>
#include <ios>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace std;

Simulation::Simulation(double Lx, double Ly, double Lz, 
					   double T_final, double dt) :
	Lx(Lx),
	Ly(Ly),
	Lz(Lz),
	domain({Lx, Ly, Lz}),
	dt(dt),
	currentTemp(currentTemp),
	T_final(T_final),
	nextOutputTime(0.0),
	outputInterval(0.1),
	fixTemp(false),
	T_current(0.0),
	randomSim(false),
	targetTemp(0.0),
	percentType1(0.0),
	N(0),
	R_min(0.5)
	{
	
	vector<Particle> particles;  // vector of particles 
	
	// needed for temperature and kinetic energy
	Physics physics;  // initialise the physics of the system
}	

	void Simulation::setN(int numParticles) {
		this->N = numParticles;
	}

	void Simulation::setFixTemp() {
		this->fixTemp = true;
	}
	
	void Simulation::setTargetTemp(double T0){
		this->targetTemp = T0;
	}
	
	
	void Simulation::TestCase1() {
		particles.clear();  // ensure particles are empty 
		array<double, 3> posVector {10.0, 10.0, 10.0};
		array<double, 3> velVector {0.0, 0.0, 0.0};
		particles.emplace_back(posVector, velVector, 0);  // construct particle in place to avoid overhead
	}

	
	/**
	* @brief initialise test case 2
	*/
	void Simulation::TestCase2(){
		particles.clear();  // ensure particles are empty 
		array<double, 3> posVector {10.0, 10.0, 10.0};
		array<double, 3> velVector {5.0, 2.0, 1.0};
		particles.emplace_back(posVector, velVector, 0);  // construct particle in place to avoid overhead
	}
	
	/**
	* @brief initialise test case 3
	*/
	void Simulation::TestCase3(){
		particles.clear();  // ensure particles are empty 
		
		array<double, 3> posVector1 {8.5, 10.0, 10.0};
		array<double, 3> velVector1 {0.0, 0.0, 0.0};
		
		array<double, 3> posVector2 {11.5, 10, 10};
		array<double, 3> velVector2 {0.0, 0.0, 0.0};
		
		particles.emplace_back(posVector1, velVector1, 0);  // construct particle in place to avoid overhead
		particles.emplace_back(posVector2, velVector2, 0);
	}
	
	/**
	* @brief initialise test case 4		
	*/
	void Simulation::TestCase4(){
		particles.clear();  // ensure particles are empty 
		
		array<double, 3> posVector1 {8.5, 11.5, 10.0};
		array<double, 3> velVector1 {0.5, 0.0, 0.0};
		
		array<double, 3> posVector2 {11.5, 8.5, 10.0};
		array<double, 3> velVector2 {-0.5, 0.0, 0.0};
		
		particles.emplace_back(posVector1, velVector1, 0);  // construct particle in place to avoid overhead
		particles.emplace_back(posVector2, velVector2, 0);
	}
	
	
	/**
	* @brief initialise test case 5
	*/
	void Simulation::TestCase5(){
		particles.clear();  // ensure particles are empty 
		
		array<double, 3> posVector1 {8.5, 11.3, 10};
		array<double, 3> velVector1 {0.5, 0.0, 0.0};
		
		array<double, 3> posVector2 {11.5, 8.7, 10.0};
		array<double, 3> velVector2 {-0.5, 0.0, 0.0};
		
		particles.emplace_back(posVector1, velVector1, 0);  // construct particle in place to avoid overhead
		particles.emplace_back(posVector2, velVector2, 0);		

	}
	
	/**
	* @brief initialise test case 6
	*/
	void Simulation::TestCase6(){
		particles.clear();  // ensure particles are empty 
		
		array<double, 3> posVector1 {8.5, 11.3, 10};
		array<double, 3> velVector1 {0.5, 0.0, 0.0};
		
		array<double, 3> posVector2 {11.5, 8.7, 10.0};
		array<double, 3> velVector2 {-0.5, 0.0, 0.0};
		
		particles.emplace_back(posVector1, velVector1, 1);  // construct particle in place to avoid overhead
		particles.emplace_back(posVector2, velVector2, 1);
	}
	
	
	/**
	* @brief initialise general dynamics
	* 
	* @param N 						number of particles to randomly spawn
	* @param percentage_type_1		percentage of particles with mass 10
	*/
	void Simulation::Random(int numParticles, double percentage_type_1){
		// set --ic-random condition
		this->randomSim = true;
		this->N = numParticles;
		// clear particles first
		particles.clear();
		
		// compute number of type 1 particles
		int N_1 =  (int) (N*percentType1/100.0);
		int N_0 = (int) N - N_1;
		
		random_device getSeed;		// Get seed for RNG
		mt19937 gen(getSeed());
		
		// Get range of values for x, y and z positions
		uniform_real_distribution<double> distX(0.0, Lx);		// X
		uniform_real_distribution<double> distY(0.0, Ly);		// Y
		uniform_real_distribution<double> distZ(0.0, Lz);		// Z
		
		uniform_real_distribution<double> vel(-0.5, 0.05);     // Velocity generator
		
		
		
		int particlesSpawned = 0;
		
		while (particlesSpawned < N) {
			array<double, 3> randomPosition {distX(gen), distY(gen), distZ(gen)};  // random position vector
			array<double, 3> randomVelocity {vel(gen), vel(gen), vel(gen)};  // random velocity vector
			int type = 0;
			if (particlesSpawned < N_1) {
				type = 1;
			}
			
			// check for distance violation
			bool Instability = false;
			
			for (const auto& p : particles) {
				double R_sq = 0.0;
				for (int dimension = 0; dimension < 3; ++dimension){
					double R_dim = randomPosition[dimension] - p.getPosVector()[dimension];
					R_sq += R_dim * R_dim;
				}
				if (R_sq < R_min) {
					Instability = true;
					break;
				}
			}
			
			// Spawn particle after checks are done 
			if (!Instability) {
				particles.emplace_back(randomPosition, randomVelocity, type);
				particlesSpawned++;
			}
			

		}
		
		
	}
	
	
	size_t Simulation::getN() const {
		return N;
	}
	
	
	const vector<Particle>& Simulation::getParticles() const{
		return particles;
	}
	
	
	double Simulation::getCurrentTime(){
		return T_current;
	}
	
	
	double Simulation::getFinalTime(){
		return T_final;
	}
	
	
	double Simulation::getKE(){
		double totalSystemEnergy = 0.0;
		for (const auto& p : particles){
			for (int i = 0; i < 3; ++i){
				totalSystemEnergy += p.getVelVector()[i] * p.getVelVector()[i] * p.getMass();
			}
		}
		return 0.5 * totalSystemEnergy;
	}
	
	
	double Simulation::getTemperature(){
		double KE = Simulation::getKE();
		return physics.calculateTemp(KE, N);
	}
	
	
	void Simulation::setTemperature(double targetTemp){
		this->currentTemp = getTemperature();
		double lambda = physics.calculateVelScale(currentTemp, targetTemp);
	
		for (auto& p : particles){
			array<double, 3> velocity = p.getVelVector();
			for (int i = 0; i < 3; ++i) {
				velocity[i] *= lambda;
			}
			p.setVelocity(velocity);
		}
	}
	
	
	void Simulation::advance(){
		
		// if user passed fixed temperature, re-scale velocities
		if (fixTemp){  
			this->setTemperature(this->targetTemp);
		}	
			
		// STEP 1: Propagate positions
		for (auto& p : this->particles) {
			
			const array<double, 3> currentPos = p.getPosVector();
			const array<double, 3> currentVel = p.getVelVector();
			const array<double, 3> currentForce = p.getForceVector();
			
			double mass = p.getMass();
			
			array<double, 3> newPosVector;			
			array<double, 3> newVelVector;
			
			// Calculate position with forward Euler
			for (int i = 0; i < 3; ++i) {
				newPosVector[i] = currentPos[i] + this->dt * currentVel[i];
				//cout << "acceleration: " << currentForce[i]/mass << endl;		// debugging purposes		
			}
			
			// Update position
			p.setPosition(newPosVector);			
			
			// check for boundary collision
			p.applyBC(this->domain); 
			
			//p.displayPos();
			
		}
		

		// STEP 2: Calculate force vector
		if (particles.size() > 1) {
			for (int i = 0; i < particles.size(); ++i){
				for (int j = i + 1; j < particles.size(); ++j){
					if (i != j) {
						int type1 = particles[i].getType();
						int type2 = particles[j].getType();

						const array<double, 3> posVec1 = particles[i].getPosVector();
						const array<double, 3> posVec2 = particles[j].getPosVector();
						
						
						// calculate distance between two particles
						double r_sq = 0.0;  // initialise squared distance
						array<double, 3> dist_vector;  // distance in a given dimension
						
						for (int dim = 0; dim < 3; ++dim){
							dist_vector[dim] = posVec1[dim] - posVec2[dim];
							r_sq += dist_vector[dim] * dist_vector[dim];
						}
						
						// force exerted on particle i by particle j
						array<double, 3> force_ij = physics.calculateForceVector(type1, type2, dist_vector, r_sq);
						
					
						
						// current force on particle i
						particles[i].zeroForce();
						array<double, 3> force_i = particles[i].getForceVector();
					
						// current force on particle j
						particles[j].zeroForce();
						array<double, 3> force_j = particles[j].getForceVector();
						
						for (int f = 0; f < 3; ++f) {
							force_i[f] += -force_ij[f];
							force_j[f] += force_ij[f];  // Newton's 3rd Law
						}	
				
						particles[i].setForce(force_i);
						particles[j].setForce(force_j);
					}
				}
			}			
		}
	
		
		// STEP 3: Update velocities
		for (auto& p : this->particles) {
			const array<double, 3> currentVel = p.getVelVector();
			const array<double, 3> currentForce = p.getForceVector();
			double mass = p.getMass();
			
			array<double, 3> newVelVector;

			for (int i = 0; i < 3; ++i) {
				newVelVector[i] = currentVel[i] + this->dt * currentForce[i]/mass;
				//cout << "acceleration: " << currentForce[i]/mass << endl;		// debugging purposes		
			}
			
			p.setVelocity(newVelVector);
			
		}

		this->T_current += dt;
		
	}
	
	void Simulation::writeOutput(const array<double, 3>& position, 
								 const array<double, 3>& velocity,
								 double kineticEnergy,
								 double currentTime, 
								 bool createNewFile,
								 bool writeEnergy,
								 int particleNumber){
		
		
		if (createNewFile){  // always output energy
			if (writeEnergy == true){
				ofstream energyData;
				energyData.open("energy.txt", ios_base::trunc | ios_base::out);
				energyData << "t E" << endl;
				energyData.close();
			}
			if (!randomSim){
				ofstream particleData; particleData.open("output.txt", ios_base::trunc | ios_base::out);
				particleData << "t pNum x y z u v w " << endl;
				particleData.close();
				writeOutput(position, velocity, kineticEnergy, currentTime, false, false, particleNumber);
			}
		} else {
			if (writeEnergy == true) {
				ofstream energyData; 
				energyData.open("energy.txt", ios_base::app | ios_base::out);  // output file for kinetic energy
				if (energyData.is_open()){
					energyData << setprecision(8) << fixed;
					energyData << currentTime << " " << this->getKE() << "\n";
				}
			}
			if (!randomSim){
				ofstream particleData; particleData.open("output.txt", ios_base::app | ios_base::out);
				if (particleData.is_open()){
					particleData << setprecision(8) << fixed;
					particleData << currentTime << " " << particleNumber << " "
								 << position[0] << " " << position[1] << " " << position[2] << " "
								 << velocity[0] << " " << velocity[1] << " " << velocity[2] << "\n";  // BACKSLASH (my keyboard doesn't do backslashes)
				}
			} 
		}
	}
	
	void Simulation::simulate(){
		
		double nextOutputTime = 0.0; // tracking output time
		bool newFile = true;  // assuem no output created
		bool writeEnergy = true;  // assume energy has not been written yet
		int numSteps = ceil(this->T_final / this->dt);  // how many steps the simulation should run for
		int writeStep = ceil(this->outputInterval / this->dt);  // how often data should be written
		
		
		// writing initial conditions
		for (int i = 0; i < particles.size(); ++i){
			
			const auto& p = particles[i];  // iterator 
			const array<double, 3>& position = p.getPosVector();
			const array<double, 3>& velocity = p.getVelVector();
			
			writeOutput(position, velocity, 
						this->getKE(), 
						this->T_current, 
						newFile,
						writeEnergy,
						i + 1);  // i + 1 to identify particles with counting numbers
			
			if (i == 0) {
				newFile = false;  // set newFile to false after first write
				writeEnergy = false;  // so energy does not get written more than once

			}
			
		}
		
		nextOutputTime = this->outputInterval;  // increment after initial conditions
		
		
		// main loop
		for (int step = 1; step <= numSteps; ++step) {
			
			this->advance();  // Forward Euler to update position and velocity
			
			// check output condition
			if (step % writeStep == 0) {
				
				for (size_t i = 0; i < particles.size(); ++i){
					const auto& p = particles[i];  // iterator 
					const array<double, 3>& position = p.getPosVector();
					const array<double, 3>& velocity = p.getVelVector();
					if (i == 0) writeEnergy = true;
					
					writeOutput(position, velocity, 
								this->getKE(), 
								this->T_current, 
								newFile,
								writeEnergy,
								i + 1);
					
					if (i == 0) writeEnergy = false;  // prevent energy from being written more than once
					
				}
				nextOutputTime += this->outputInterval;
			}
			
			for (const auto& p: particles){
			//p.displayPos(); // debugging purposes
			}
			
		}
	
	}
	
	


