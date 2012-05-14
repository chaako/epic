/*
 * Orbit.cpp
 *
 *  Created on: Feb 10, 2012
 *      Author: chaako
 */

#include "epic.h"
#include "Orbit.h"

//Orbit::Orbit(Eigen::Vector3d inputPosition, Eigen::Vector3d inputVelocity,
//		iBase_EntityHandle inputElement) {
//}

Orbit::Orbit(Mesh *inputMesh_ptr, iBase_EntityHandle inputNode,
		Eigen::Vector3d inputVelocity, double inputCharge) {
	mesh_ptr = inputMesh_ptr;
	initialNode = inputNode;
	// TODO: perhaps find a random adjacent tet rather than give node
	currentElement = inputNode;
	initialPosition = inputMesh_ptr->getCoordinates(inputNode);
	initialVelocity = inputVelocity;
	charge = inputCharge;
}

Orbit::~Orbit() {
	// TODO Auto-generated destructor stub
}

void Orbit::integrate(ElectricField& electricField,
		PotentialField& potentialField,
		Field<int>& faceTypeField, CodeField& vertexTypeField, FILE *outFile) {
	Eigen::Vector3d eField = electricField.getField(initialNode);
	int nVertices=4;
	std::vector<Eigen::Vector3d> eFields(nVertices), vertexVectors(nVertices);
	std::vector<iBase_EntityHandle> vertices(nVertices);
	// TODO: need some clever way to set tMax and/or detect trapped orbits
	double dt=std::min(0.01,0.01/initialVelocity.norm()), tMax=100;
	Eigen::Vector3d currentPosition = initialPosition;
	Eigen::Vector3d currentVelocity = initialVelocity;
	// TODO: shouldn't hard-code quasi-neutral operation
	double phiSurface = -4;
	vertexType = vertexTypeField.getField(initialNode);
	Eigen::Vector3d vertexNormalVector;
	Eigen::Vector3d initialNormalVelocity;
	if (vertexType==4 && charge<0.) {
		vertexNormalVector =
				mesh_ptr->getVertexNormalVector(initialNode, faceTypeField);
		Eigen::Vector3d coords = mesh_ptr->getCoordinates(initialNode);
		initialNormalVelocity =
				currentVelocity.dot(vertexNormalVector)*vertexNormalVector;
		// TODO: could do something like below to get distribution at surface
//		currentVelocity += sqrt(2.*charge*phiSurface)*vertexNormalVector;
	}
	// Don't integrate orbit if doesn't have enough energy to escape potential
	// TODO: this should be refined
//	if (0.5*pow(initialVelocity.norm(),2.)+0.22 <
	if ((0.5*pow(initialVelocity.norm(),2.) +
			charge*potentialField.getField(initialNode)) < 0.) {
		negativeEnergy = true;
		tMax = 0.;
	} else {
		negativeEnergy = false;
	}
	// For second order leap-frog, offset position from velocity in time
	currentPosition -= currentVelocity*dt/2.;
	// TODO: treat inwards electrons in a better way?
	if (vertexType==4 &&
			0.5*pow(initialNormalVelocity.norm(),2.)<charge*phiSurface &&
			currentVelocity.dot(vertexNormalVector)<0.) {
		currentVelocity -= 2.*initialNormalVelocity;
	}
	bool inNewTet = true;
	int nSteps=0, nNewTet=0;
//	std::cout << "Initial radius=" << currentPosition.norm() << std::endl;
	bool isTet=false;
	bool firstStep=true;
	for (double t=0; t<tMax; t+=dt) {
		nSteps++;
		Eigen::Vector3d previousPosition = currentPosition;
		currentPosition += dt*currentVelocity;
//		clock_t startClock = clock(); // timing
		if (!firstStep)
			inNewTet = !mesh_ptr->checkIfInTet(currentPosition, vertexVectors);
		firstStep=false;
//		clock_t endClock = clock(); // timing
//		extern_checkIfInNewTet += endClock-startClock; // timing
		if (inNewTet) {
			nNewTet++;
			bool foundTet=false;
//			if  (vertexVectors.size()==nVertices)
//				isTet = true;
			clock_t startClock = clock(); // timing
			iBase_EntityHandle previousElement = currentElement;
			currentElement = mesh_ptr->findTet(previousPosition,
					currentPosition, previousElement, &foundTet, isTet);
			clock_t endClock = clock(); // timing
			extern_findTet += endClock-startClock; // timing
			// TODO: should handle failure to find tet in some way
			if (foundTet==false) {
				iBase_EntityHandle faceCrossed = mesh_ptr->findFaceCrossed(
						previousElement, previousPosition, currentPosition);
				// TODO: grazing orbits don't enter domain in first time-step
				int faceType;
				Eigen::Vector3d normalVector, normalVelocity;
				if (faceCrossed!=NULL) {
					faceType = faceTypeField.getField(faceCrossed);
					normalVector = mesh_ptr->getNormalVector(faceCrossed,
									previousPosition);
					normalVelocity =
							currentVelocity.dot(normalVector)*normalVector;
					finalPotential = 0;
					vertices = mesh_ptr->getVertices(faceCrossed);
					for (int i=0; i<vertices.size(); i++) {
						// TODO: should use point where left domain here
						finalPotential += 1./3.*potentialField.getField(vertices[i]);
					}
					// TODO: shouldn't hard-code boundary code
				} else {
					faceType = 0;
				}
				finalFaceType = faceType;
				if (faceType==4 && 0.5*pow(normalVelocity.norm(),2.)<charge*phiSurface) {
					foundTet = true;
					currentElement = previousElement;
					// TODO: resetting position isn't quite right
					currentPosition = previousPosition;
					// TODO: revisit this in magnetized case
					currentVelocity -= 2.*normalVelocity;
				} else {
					break;
				}
			}
			isTet=true;

			vertices = mesh_ptr->getVertices(currentElement);
			assert( vertexVectors.size()==vertices.size() &&
					eFields.size()==vertices.size() );
			for (int i=0; i<vertices.size(); i++) {
				vertexVectors[i] = mesh_ptr->getCoordinates(vertices[i]);
				eFields[i] = electricField.getField(vertices[i]);
			}
		}
//		if (t>=tMax-dt)
//			std::cout << "orbit reached tMax" << std::endl;

		assert(vertexVectors.size()==nVertices);
		std::vector<double> vertexWeights = mesh_ptr->getVertexWeights(currentPosition,
				vertexVectors);
		Eigen::Vector3d currentAcceleration(0.,0.,0.);
		assert(eFields.size()==vertexWeights.size());
		for (int i=0; i<vertexWeights.size(); i++) {
			currentAcceleration += charge*eFields[i]*vertexWeights[i];
		}
//		currentAcceleration = -currentPosition/pow(currentPosition.norm(),3.);
		double eFieldR = currentAcceleration.dot(currentPosition)/
				currentPosition.norm();
		currentVelocity += dt*currentAcceleration;
//		double potential = potentialField.getField(currentPosition, currentElement);
//		Eigen::Vector3d velocityAtPosition = currentVelocity - 1./2.*dt*currentAcceleration;
//		double energy = 1./2.*pow(velocityAtPosition.norm(),2.) + charge*potential;
//		if (outFile) {
//			fprintf(outFile, "%f %f %f %f\n", currentPosition[0], currentPosition[1],
//					currentPosition[2], energy);
//		}
	}
//	std::cout << "Final radius=" << currentPosition.norm() << std::endl;
	finalPosition = currentPosition;
	// TODO: correct for time-step offset?
	finalVelocity = currentVelocity;
}

void Orbit::integrate(PotentialField& potentialField, ElectricField& electricField,
		Field<int>& faceTypeField, CodeField& vertexTypeField, FILE *outFile) {
	// TODO: need some clever way to set tMax and/or detect trapped orbits
	double dt=std::min(0.01,0.01/initialVelocity.norm()), tMax=100;
	Eigen::Vector3d currentPosition = initialPosition;
	Eigen::Vector3d currentVelocity = initialVelocity;
	// TODO: shouldn't hard-code quasi-neutral operation
	double phiSurface = -4;
	vertexType = vertexTypeField.getField(initialNode);
	Eigen::Vector3d vertexNormalVector;
	Eigen::Vector3d initialNormalVelocity;
	if (vertexType==4 && charge<0.) {
		vertexNormalVector =
				mesh_ptr->getVertexNormalVector(initialNode, faceTypeField);
		Eigen::Vector3d coords = mesh_ptr->getCoordinates(initialNode);
		initialNormalVelocity =
				currentVelocity.dot(vertexNormalVector)*vertexNormalVector;
		// TODO: could do something like below to get distribution at surface
//		currentVelocity += sqrt(2.*charge*phiSurface)*vertexNormalVector;
	}
	// Don't integrate orbit if doesn't have enough energy to escape potential
	// TODO: this should be refined
//	if (0.5*pow(initialVelocity.norm(),2.)+0.22 <
	if ((0.5*pow(initialVelocity.norm(),2.) +
			charge*potentialField.getField(initialNode)) < 0.) {
		negativeEnergy = true;
		tMax = 0.;
	} else {
		negativeEnergy = false;
	}
	// For second order leap-frog, offset position from velocity in time
	currentPosition -= currentVelocity*dt/2.;
	// TODO: treat inwards electrons in a better way?
	if (vertexType==4 &&
			0.5*pow(initialNormalVelocity.norm(),2.)<charge*phiSurface &&
			currentVelocity.dot(vertexNormalVector)<0.) {
		currentVelocity -= 2.*initialNormalVelocity;
	}
	int nSteps=0;
	double potential=0.;
	bool endLoop=false;
	bool foundTet=false;
	for (double t=0.; t<tMax; t+=dt) {
		nSteps++;
		Eigen::Vector3d previousPosition = currentPosition;
		currentPosition += dt*currentVelocity;
		iBase_EntityHandle previousElement = currentElement;
		Eigen::Vector3d currentAcceleration(0.,0.,0.);
		// TODO: replace this hack
		if (t==0.) {
			currentElement = mesh_ptr->findTet(previousPosition,
					currentPosition, initialNode, &foundTet, false);
			// TODO: figure out why sometimes throw here (grazing orbits?)
//			if (!foundTet)
//				throw;
			int dimension=mesh_ptr->getEntityDimension(currentElement);
			if (dimension!=iBase_REGION)
				foundTet = false;
		}
//		if (!foundTet)
//			break;
		if (foundTet) {
			try {
				// TODO: set order through input parameter?
				int interpolationOrder = 1;
//				currentAcceleration = charge*
//						electricField.getField(currentPosition, &currentElement);
				potential = potentialField.getField(currentPosition, &currentElement,
						interpolationOrder);
				if (isnan(potential))
					potential = potentialField.getField(currentPosition,
							&currentElement, 1);
				// TODO: hard-coding dimension here...
				for (int i=0; i<3; i++) {
					Eigen::Vector3d perturbedPosition = currentPosition +
							Eigen::Vector3d::Unit(i)*DELTA_LENGTH;
					double perturbedPotential = potentialField.getField(
							perturbedPosition, &currentElement, interpolationOrder);
					// TODO: track down why perturbedPotential sometimes is NaN
					if (isnan(perturbedPotential)) {
						perturbedPotential = potentialField.getField(
								perturbedPosition, &currentElement, 1);
					}
//					if (isnan(potential) || isnan(perturbedPotential)) {
//						std::cout << potential << " " << perturbedPotential <<
//						" " << currentPosition.transpose() << std::endl;
//						throw;
//					}
					currentAcceleration[i] =
							-charge*(perturbedPotential-potential)/DELTA_LENGTH;
				}
			} catch (int signal) {
				switch (signal) {
				case OUTSIDE_DOMAIN:
					foundTet = false;
					break;
				default:
					// TODO: handle other exceptions?
					throw;
					break;
				}
			} catch (...) {
				// TODO: handle other exceptions?
				throw;
			}
		}

		if (foundTet==false) {
			iBase_EntityHandle faceCrossed = mesh_ptr->findFaceCrossed(
					previousElement, previousPosition, currentPosition);
			// TODO: grazing orbits don't enter domain in first time-step
			int faceType;
			Eigen::Vector3d normalVector(0.,0.,0.), normalVelocity(0.,0.,0.);
			if (faceCrossed!=NULL) {
				faceType = faceTypeField.getField(faceCrossed);
				normalVector = mesh_ptr->getNormalVector(faceCrossed,
								previousPosition);
				normalVelocity =
						currentVelocity.dot(normalVector)*normalVector;
				finalPotential = 0.;
				std::vector<iBase_EntityHandle> vertices =
						mesh_ptr->getVertices(faceCrossed);
				for (int i=0; i<vertices.size(); i++) {
					// TODO: should use point where left domain here
					finalPotential += 1./3.*potentialField.getField(vertices[i]);
				}
				// TODO: shouldn't hard-code boundary code
			} else {
				faceType = 0;
			}
			finalFaceType = faceType;
			if (faceType==4 && 0.5*pow(normalVelocity.norm(),2.)<charge*phiSurface) {
				currentElement = previousElement;
				foundTet = true;
				// TODO: resetting position isn't quite right
				currentPosition = previousPosition;
				// TODO: revisit this in magnetized case
				currentVelocity -= 2.*normalVelocity;
			} else {
				endLoop=true;
			}
		}

		double eFieldR = currentAcceleration.dot(currentPosition)/
				currentPosition.norm();
		currentVelocity += dt*currentAcceleration;
		Eigen::Vector3d velocityAtPosition = currentVelocity - 1./2.*dt*currentAcceleration;
		double energy = 1./2.*pow(velocityAtPosition.norm(),2.) + charge*potential;
		if (outFile) {
			fprintf(outFile, "%f %f %f %f\n", currentPosition[0], currentPosition[1],
					currentPosition[2], energy);
		}
		if (endLoop)
			break;
	}
//	std::cout << "Final radius=" << currentPosition.norm() << std::endl;
	finalPosition = currentPosition;
	// TODO: correct for time-step offset?
	finalVelocity = currentVelocity;
}
