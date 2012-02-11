/*
 * Orbit.h
 *
 *  Created on: Feb 10, 2012
 *      Author: chaako
 */

#ifndef ORBIT_H_
#define ORBIT_H_

class Orbit {
public:
	Orbit(Eigen::Vector3d initialPosition, Eigen::Vector3d initialVelocity);
	virtual ~Orbit();
	void integrate(ElectricField eField);
	Eigen::Vector3d initialPosition;
	Eigen::Vector3d initialVelocity;
	Eigen::Vector3d finalPosition;
	Eigen::Vector3d finalVelocity;
};

#endif /* ORBIT_H_ */
