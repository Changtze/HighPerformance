#ifndef PHYSICS_H
#define PHYSICS_H

#include <array>
class Particle;  // forward delcaration


class Physics {
	
public:
	static constexpr long double kB = 0.8314459920816467;
	double epsilon[2][2];  
	double sigma[2][2];

	/**
	* @brief Class constructor to initialise system physics
	*/
	Physics();  													
	
	
	/**
	* @brief this class function handles the instantiation of Lennard-Jones parameters
	*/
	void initialise_LJ();
	

	/**
	* @brief Retrieves relevant epsilon parameter
	* 
	* @return epsilon_ij
	*/
	double getEpsilon(int type1, int type2) const;
	
	
	/**
	* @brief Retrieves relevant sigma parameter
	* 
	* @return sigma_ij
	*/
	double getSigma(int type1, int type2) const;
	
	
	/**
	* @brief this function handles calculating the Lennard-Jones derivative in a given dimension
	* @param type1		 particle 1 type
	* @param type2		 particle 2 type
	* @param r			 distance vector between particles
	* @param distance	 L2-norm of r (i.e. Euclidean distance between particles)
	* @param dimension   0: x, 1: y, 2: z
	*/
	double computeLJDeriv(int type1, int type2,
						  double distance, std::array<double, 3>& r, 
						  int dimension) const;
	
	
	/**
	* @brief this function calculates the force vector experienced by particle 1 due to particle 2
	* 
	* @param type1		 particle 1 type
	* @param type2		 particle 2 type
	* @param r			 distance vector between particles
	* @param distance	 L2-norm of r (i.e. Euclidean distance between particles)
	* 
	* @return force vector applied to particle 
	*/
	std::array<double, 3> calculateForceVector(int type1, int type2,
											   std::array<double, 3> r, 
											   double distance) const ;
	
	/**
	* @brief calculates temperature of the system
	* 
	* @param avgKE		average kinetic energy of particles
	* @param N			number of particles in system
	* @param physics	physics object
	*/
	double calculateTemp(double avgKE, int N);
	
	
	/**
	* @brief calculates velocity scale lambda
	* 
	* @param currentTemp		current temperature of system
	* @param userTemp			target temperature to fix system at
	*/
	double calculateVelScale(double currentTemp, double userTemp);
	
};

#endif  // PHYSICS_H