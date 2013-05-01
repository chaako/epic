/*
 * MinimumBasisFunction.cpp
 *
 *  Created on: Apr 28, 2013
 *      Author: chaako
 */

#include "epic.h"
#include "MinimumBasisFunction.h"

MinimumBasisFunction::MinimumBasisFunction(Mesh* inputMesh_ptr,
		boost::array<Eigen::Matrix<double,NDIM,1>, 2> *inputPositionAndVelocity_ptr,
		VelocityAndAcceleration<NDIM> *inputVelocityAndAcceleration_ptr,
		DriftStepper *inputTimestepper_ptr) {
	mesh_ptr = inputMesh_ptr;
	positionAndVelocity_ptr = inputPositionAndVelocity_ptr;
	velocityAndAcceleration_ptr = inputVelocityAndAcceleration_ptr;
	timestepper_ptr = inputTimestepper_ptr;
}

MinimumBasisFunction::~MinimumBasisFunction() {
	// TODO Auto-generated destructor stub
}

double MinimumBasisFunction::operator()(double dt) {
	boost::array<Eigen::Matrix<double,NDIM,1>, 2> positionAndVelocity;
	positionAndVelocity[0] = (*positionAndVelocity_ptr)[0];
	positionAndVelocity[1] = (*positionAndVelocity_ptr)[1];
	timestepper_ptr->do_step(boost::ref(*velocityAndAcceleration_ptr),
			positionAndVelocity, 0., dt);
	return mesh_ptr->minimumBasisFunction(positionAndVelocity[0],
			velocityAndAcceleration_ptr->currentRegionIndex);
}
