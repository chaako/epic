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
			ElectricField &inputElectricField,
			double inputCharge, entHandle initialNode, vect3d initialVelocity,
			bool inputOnlyUsePotential=true) :
		potentialField(inputPotentialField),
		electricField(inputElectricField),
		charge(inputCharge),
		onlyUsePotential(inputOnlyUsePotential)
	{
		currentPosition=potentialField.mesh_ptr->getCoordinates(initialNode);
		for (int i=0; i<NDIM; i++) {
			currentVelocity[i] = 0.;
			currentAcceleration[i] = 0.;
		}
		// TODO: determine if NaN*0=NaN is a problem with uninitalized memory (it appears so...)
//		currentVelocity*=0.;
//		currentAcceleration*=0.;
		foundTet = false;
		currentRegionIndex=-1;
		try {
			currentElement = potentialField.mesh_ptr->findStartingTet(currentPosition,
					initialVelocity, initialNode);
			currentRegionIndex =
					potentialField.mesh_ptr->indicesOfEntities[currentElement];
			foundTet = true;
		} catch (int numberOfRegionsWithinTolerance) {
			if (numberOfRegionsWithinTolerance==0) {
				// TODO: come up with better handling of this
				currentElement = initialNode;
			} else {
				throw;
			}
		}
//		int dimension=potentialField.mesh_ptr->getEntityDimension(currentElement);
//		assert(dimension==3);
//		if (!foundTet)
//			throw int(OUTSIDE_DOMAIN);
		// TODO: Handle this more transparently
		if (onlyUsePotential) {
			interpolationOrder = INTERPOLATIONORDER;
		} else {
			interpolationOrder = 1;
		}
		currentPotential = potentialField.getField(initialNode);
	}

	// TODO: should perhaps be consistent about ptrs vs refs
	PotentialField& potentialField;
	ElectricField& electricField;
	double charge;
	bool onlyUsePotential;
	entHandle currentElement;
	int currentRegionIndex;
	bool foundTet;
	int interpolationOrder;

	void operator()(const state_type& x, state_type& dxdt, double t=-1.) {
		// TODO: need DELTA_LENGTH...make header hierarchy?
		double potential;
//		std::cout << currentPosition.transpose() << std::endl;
		currentVelocity = x[1];
		// TODO: for efficiency should find way to only calc accel when needed
//		if (x[0]!=currentPosition) {
			currentPosition = x[0];
			if (onlyUsePotential) {
				potentialField.evalFieldAndDeriv(&potential, &currentAcceleration,
						x[0], &currentRegionIndex,
						interpolationOrder);
				currentAcceleration *= -1.;
				currentPotential = potential;
			} else {
				electricField.evalField(&currentAcceleration,
						x[0], &currentRegionIndex,
						interpolationOrder);
				potentialField.evalField(&currentPotential,
						x[0], &currentRegionIndex,
						interpolationOrder);
			}
			currentElement = potentialField.mesh_ptr->
					entitiesVectors[iBase_REGION][currentRegionIndex];
//			// TODO: Determine if better to include drift in potential
//			//       (conductor shields it, thermal velocities could be different, etc.)
//			currentAcceleration -= E;
			currentAcceleration *= charge;
//			potential = potentialField.getField(x[0], &currentElement,
//					interpolationOrder);
//			for (int i=0; i<N; i++) {
//				Eigen::Matrix<double,N,1> perturbedPosition = x[0] +
//						Eigen::Matrix<double,N,1>::Unit(i)*DELTA_LENGTH;
//				double perturbedPotential = potentialField.getField(
//						perturbedPosition, &currentElement, interpolationOrder);
//				currentAcceleration[i] =
//						-charge*(perturbedPotential-potential)/DELTA_LENGTH;
//			}
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
	double currentPotential;
};

#endif /* VELOCITYANDACCELERATION_H_ */
