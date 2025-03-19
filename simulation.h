#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include <string>
#include <array>
#include "particle.h"
#include "physics.h"

class Particle;  // forward declaration of particle class


/**
* @brief This class manages the molecular dynamics simulation
*
* Time stepping, initialisation, force calculations and data output
*/

class Simulation{	
private:
	const std::array<double, 3> domain;  // vector containing domain boundaries
	double Lx;  // x-domain length
	double Ly;  // y-domain length
	double Lz;  // z-domain length
	double dt;  // simulation time step
	double T_current;  // current simulation time
	double T_final; // final time
	size_t N; // number of particles
	std::vector<Particle> particles;  // vector of particles 
	double currentTemp;  
	double targetTemp;  // only if temperature control is enabled
	bool fixTemp;  // indicate if temperature is fixed by user
	bool randomSim;  // check for --ic-random
	double percentType1;  // percentage of type 1 particles
	
	const double R_min;  // minimum distance between two particles
	const double outputInterval;  // watch out for this
	
	double nextOutputTime;  // next time to write to text file
	
	
	
public:
	Physics physics;  // initialise the physics of the system
	/**
	* @brief Parameterised constructor
	* 
	* @param Lx				domain x-length (Angstroms)
	* @param Ly 			domain y-length (Angstroms)
	* @param Lz				domain z-length (Angstroms)
	* @param T_final		final time (arbitrary units)
	* @param Temperature	system temperature (left as -1 if not specified by user)
	*/
	Simulation(double Lx, double Ly, 
			   double Lz, double T_final, 
			   double dt);
			   
	void setN(int numParticles);
	
	/**
	* @brief sets fixTemp to true or false based on the presence of --temp argument
	*/
	void setFixTemp();
	
	/**
	* @brief sets the target temperature provided by --temp
	*/
	void setTargetTemp(double T0);
	
	
	/**
	* @brief initialise test case 1
	*/
	void TestCase1();
	
	/**
	* @brief initialise test case 2
	*/
	void TestCase2();
	
	/**
	* @brief initialise test case 3
	*/
	void TestCase3();
	
	/**
	* @brief initialise test case 4
	*/
	void TestCase4();
	
	/**
	* @brief initialise test case 5
	*/
	void TestCase5();
	
	/**
	* @brief initialise test case 6
	*/
	void TestCase6();
	
	/**
	* @brief initialise general dynamics
	* 
	* @param N 						number of particles to randomly spawn
	* @param percentage_type_1		percentage of particles with mass 10
	*/
	void Random(int N, double percentage_type_1);
	
	/**
	* @brief updates the system by a time step using forward Euler
	*/
	void advance(); 
	
	
	/**
	* @brief runs a simulation for the full specified time
	* 
	* @param T_final				end time of the simulation
	* @param outputInterval			how often data should be outputted
	*/
	void simulate();
	
	/**
	* @brief writes energy and particle data to two files
	* 
	* @param positon				position vector 
	* @param velocity				 	velocity vector
	* @param kineticEnergy			K.E. of system
	* @param currentTime			current time in simulation
	* @param createNewFile			whether to create a new file (i.e. overwrite)
	* @param particleNumber			particle identifier
	*/
	void writeOutput(const std::array<double, 3>& position, 
								 const std::array<double, 3>& velocity,
								 double kineticEnergy,
								 double currentTime, 
								 bool createNewFile,
								 bool writeOutput,
								 int particleNumber);
	
	
	/**
	* @brief get number of particles in system
	*/
	size_t getN() const;
	
	/**
	* @brief get current particles in the system
	* 
	* @return Reference to vector of particles
	*/
	const std::vector<Particle>& getParticles() const;
	
	
	/**
	* @brief gets current time of system
	* 
	* @return current time value
	*/
	double getCurrentTime();
	
	/**
	* @brief gets final time of system
	* 
	* @return final time value
	*/
	double getFinalTime();
	
	/**
	* @brief calculates average kinetic energy of system
	* 
	* @return kinetic energy as given in the handout
	*/
	double getKE();
	
	/**
	* @brief gets temperature of system
	* 
	* @return temperature value of system in Kelvin
	*/
	double getTemperature();
	
	/**
	* @brief sets temperature of system
	*/
	void setTemperature(double targetTemp);
	
	
};


#endif