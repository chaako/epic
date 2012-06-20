/*
 * VelocityAndAcceleration.h
 *
 *  Created on: Jun 4, 2012
 *      Author: chaako
 */

#ifndef VELOCITYANDACCELERATION_H_
#define VELOCITYANDACCELERATION_H_

#include "typesAndDefinitions.h"
//#include <boost/array.hpp>
//#include "iMesh.h"
//#include "Eigen/Dense"

#include "Mesh.h"
#include "Field.h"

template<int N> class VelocityAndAcceleration {
public:
	typedef boost::array<Eigen::Matrix<double,N,1>, 2> state_type;

	VelocityAndAcceleration(PotentialField &inputPotentialField,
			double inputCharge, entHandle initialNode) :
		potentialField(inputPotentialField),
		charge(inputCharge)
	{
		currentPosition=potentialField.mesh_ptr->getCoordinates(initialNode);
		currentVelocity*=0.;
		currentAcceleration*=0.;
		foundTet = false;
		currentElement = potentialField.mesh_ptr->findTet(currentPosition,
				currentPosition, initialNode, &foundTet, false);
//		int dimension=potentialField.mesh_ptr->getEntityDimension(currentElement);
//		assert(dimension==3);
//		if (!foundTet)
//			throw int(OUTSIDE_DOMAIN);
		interpolationOrder = INTERPOLATIONORDER;
	}

	// TODO: should perhaps be consistent about ptrs vs refs
	PotentialField& potentialField;
	double charge;
	entHandle currentElement;
	bool foundTet;
	int interpolationOrder;

	void operator()(const state_type& x, state_type& dxdt) {
		// TODO: need DELTA_LENGTH...make header hierarchy?
		double potential;
//		std::cout << currentPosition.transpose() << std::endl;
		currentVelocity = x[1];
		// TODO: for efficiency should find way to only calc accel when needed
//		if (x[0]!=currentPosition) {
			currentPosition = x[0];
			potential = potentialField.getField(x[0], &currentElement,
					interpolationOrder);
			for (int i=0; i<N; i++) {
				Eigen::Matrix<double,N,1> perturbedPosition = x[0] +
						Eigen::Matrix<double,N,1>::Unit(i)*DELTA_LENGTH;
				double perturbedPotential = potentialField.getField(
						perturbedPosition, &currentElement, interpolationOrder);
				currentAcceleration[i] =
						-charge*(perturbedPotential-potential)/DELTA_LENGTH;
			}
//		}
		dxdt[0] = currentVelocity;
		dxdt[1] = currentAcceleration;
//		std::cout << currentPosition.transpose() << std::endl;
//		std::cout << currentVelocity.transpose() << std::endl;
//		std::cout << currentAcceleration.transpose() << std::endl << std::endl;
	}

	Eigen::Matrix<double,N,1> currentPosition;
	Eigen::Matrix<double,N,1> currentVelocity;
	Eigen::Matrix<double,N,1> currentAcceleration;
};

#endif /* VELOCITYANDACCELERATION_H_ */
