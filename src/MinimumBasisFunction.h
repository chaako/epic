/*
 * MinimumBasisFunction.h
 *
 *  Created on: Apr 28, 2013
 *      Author: chaako
 */

#ifndef MINIMUMBASISFUNCTION_H_
#define MINIMUMBASISFUNCTION_H_

#include "typesAndDefinitions.h"

//class Mesh;
//template<int n> class VelocityAndAcceleration;
//class DriftStepper;
#include "Mesh.h"
//#include "Field.h"
#include "VelocityAndAcceleration.h"
#include "Stepper.h"

class MinimumBasisFunction {
public:
	MinimumBasisFunction(Mesh* inputMesh_ptr,
			boost::array<Eigen::Matrix<double,NDIM,1>, 2> *positionAndVelocity_ptr,
			VelocityAndAcceleration<NDIM> *inputVelocityAndAcceleration_ptr,
			DriftStepper *inputTimestepper_ptr);
	virtual ~MinimumBasisFunction();

	double operator()(double dt);

	Mesh *mesh_ptr;
	// TODO: pointer may not be fastest way to do this
	boost::array<Eigen::Matrix<double,NDIM,1>, 2> *positionAndVelocity_ptr;
	VelocityAndAcceleration<NDIM> *velocityAndAcceleration_ptr;
	// TODO: should have parent Stepper class here
	DriftStepper *timestepper_ptr;

};

#endif /* MINIMUMBASISFUNCTION_H_ */
