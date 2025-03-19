#ifndef PARTICLE_H
#define PARTICLE_H

#include <array>

/**
 * @brief This class represents a single particle in the system
 * 
 * It stores the position and velocity of a particle
 * It also stores the force vector experienced by the particle
 * Attributes include mass and type
*/

class Particle {  
private:
	std::array<double, 3> posVector;
	std::array<double, 3> velVector;
	std::array<double, 3> forceVector;
	int type;
	double mass; 
public:
	/**
	* @brief Default constructor
	*/
	Particle(); 
	
	
	/**
	* @brief Parameterised constructor
	* 
	* @param initPosVector		position vector in (x, y, z) form
	* @param initVelVector		velocity vecotr in (u, v, w) form
	* @param type			Particle type (0 or 1)
	*/
	
	Particle(const std::array<double, 3>& position,
			 const std::array<double, 3>& velocity,
			 int type);
			 
		/**
		 * @brief Displays force vector in (Fx, Fy, Fz) form
		 */
		void displayForce() const;	
			
		/**
		 * @brief Displays velocity vector in (u, v, w) form
		 */
		void displayVel() const;	
		
		/**
		 * @brief Displays position vector in (x, y, z) form
		 */
		void displayPos() const;
		
		/**
		 * @brief Getter function: particle mass
		 * @return particle mass
		 */
		 double getMass() const;
		 
		 
		 /**
		 * @brief Getter function: particle type
		 * @return particle type
		 */
		 int getType() const;
		 
		 
		 /**
		 * @brief Getter function: position vector
		 * @return position vector reference
		 */
		 const std::array<double, 3>& getPosVector() const;
		
		 
		 /**
		 * @brief Getter function: velocity vector
		 * @return velocity vector reference
		 */
		 const std::array<double, 3>& getVelVector() const;
		 
		 /**
		 * @brief Getter function: force vector experienced by particle
		 * @return force vector reference
		 */
		 const std::array<double, 3>& getForceVector() const;
		 
		 
		 /**
		 * @brief Setter function: set new position vector
		 * @param desired position
		 */
		 void setPosition(const std::array<double, 3>& newPos);
		 
		 
		 /**
		 * @brief Setter function: set new velocity vector
		 * @param desired velocity
		 */
		 void setVelocity(const std::array<double, 3>& newVel);
		 
		 
		 /**
		 * @brief Setter function: set new force vector
		 * @param desired force
		 */
		 void setForce(const std::array<double, 3>& newForce);
		 
		 /**
		 * @brief Setter function: set force to 0
		 */
		 void zeroForce();
		 
		 /**
		 * @brief Setter function: apply boundary conditions of reflection
		 * @param domain size
		 */
		 void applyBC(const std::array<double, 3>& domain);
		
};
#endif