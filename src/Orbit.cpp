/*
 * Orbit.cpp
 *
 *  Created on: Feb 10, 2012
 *      Author: chaako
 */

#include "epic.h"
#include "Orbit.h"

//Orbit::Orbit(vect3d inputPosition, vect3d inputVelocity,
//		entHandle inputElement) {
//}

Orbit::Orbit(Mesh *inputMesh_ptr, entHandle inputNode,
		vect3d inputVelocity, double inputCharge) {
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
	vect3d eField = electricField.getField(initialNode);
	int nVertices=4;
	vector<vect3d> eFields(nVertices), vertexVectors(nVertices);
	vector<entHandle> vertices(nVertices);
	// TODO: need some clever way to set tMax and/or detect trapped orbits
	double dt=min(0.01,0.01/initialVelocity.norm()), tMax=100;
	vect3d currentPosition = initialPosition;
	vect3d currentVelocity = initialVelocity;
	// TODO: shouldn't hard-code quasi-neutral operation
	double phiSurface = -4;
	vertexType = vertexTypeField.getField(initialNode);
	vect3d vertexNormalVector;
	vect3d initialNormalVelocity;
	if (vertexType==4 && charge<0.) {
		vertexNormalVector =
				mesh_ptr->getVertexNormalVector(initialNode, faceTypeField);
		vect3d coords = mesh_ptr->getCoordinates(initialNode);
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
//	cout << "Initial radius=" << currentPosition.norm() << endl;
	bool isTet=false;
	bool firstStep=true;
	for (double t=0; t<tMax; t+=dt) {
		nSteps++;
		vect3d previousPosition = currentPosition;
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
			entHandle previousElement = currentElement;
			currentElement = mesh_ptr->findTet(previousPosition,
					currentPosition, previousElement, &foundTet, isTet);
			clock_t endClock = clock(); // timing
			extern_findTet += endClock-startClock; // timing
			// TODO: should handle failure to find tet in some way
			if (foundTet==false) {
				entHandle faceCrossed = mesh_ptr->findFaceCrossed(
						previousElement, previousPosition, currentPosition);
				// TODO: grazing orbits don't enter domain in first time-step
				int faceType;
				vect3d normalVector, normalVelocity;
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
//			cout << "orbit reached tMax" << endl;

		assert(vertexVectors.size()==nVertices);
		vector<double> vertexWeights = mesh_ptr->getVertexWeights(currentPosition,
				vertexVectors);
		vect3d currentAcceleration(0.,0.,0.);
		assert(eFields.size()==vertexWeights.size());
		for (int i=0; i<vertexWeights.size(); i++) {
			currentAcceleration += charge*eFields[i]*vertexWeights[i];
		}
//		currentAcceleration = -currentPosition/pow(currentPosition.norm(),3.);
		double eFieldR = currentAcceleration.dot(currentPosition)/
				currentPosition.norm();
		currentVelocity += dt*currentAcceleration;
//		double potential = potentialField.getField(currentPosition, currentElement);
//		vect3d velocityAtPosition = currentVelocity - 1./2.*dt*currentAcceleration;
//		double energy = 1./2.*pow(velocityAtPosition.norm(),2.) + charge*potential;
//		if (outFile) {
//			fprintf(outFile, "%f %f %f %f\n", currentPosition[0], currentPosition[1],
//					currentPosition[2], energy);
//		}
	}
//	cout << "Final radius=" << currentPosition.norm() << endl;
	finalPosition = currentPosition;
	// TODO: correct for time-step offset?
	finalVelocity = currentVelocity;
}

void Orbit::integrate(PotentialField& potentialField, ElectricField& electricField,
		Field<int>& faceTypeField, CodeField& vertexTypeField, FILE *outFile) {
	// TODO: need some clever way to set tMax and/or detect trapped orbits
//	double dt=min(0.005,0.005/initialVelocity.norm()), tMax=100;
	double dt=min(0.01,0.01/initialVelocity.norm()), tMax=100;
//	double dt=min(0.02,0.02/initialVelocity.norm()), tMax=100;
	vect3d currentPosition = initialPosition;
	vect3d currentVelocity = initialVelocity;
	// TODO: shouldn't hard-code quasi-neutral operation
	double phiSurface = -4;
	vertexType = vertexTypeField.getField(initialNode);
	vect3d vertexNormalVector;
	vect3d initialNormalVelocity;
	if (vertexType==4 && charge<0.) {
		vertexNormalVector =
				mesh_ptr->getVertexNormalVector(initialNode, faceTypeField);
		vect3d coords = mesh_ptr->getCoordinates(initialNode);
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

	boost::array<Eigen::Matrix<double,NDIM,1>, 2> positionAndVelocity;
	boost::array<Eigen::Matrix<double,NDIM,1>, 2> positionAndVelocityOut;
//	positionAndVelocity[0] = currentPosition;
//	positionAndVelocity[1] = currentVelocity;
	VelocityVerletStepper<NDIM> timestepper;
//	boost::numeric::odeint::runge_kutta4<boost::array<vect3d,2> >
//		timestepper;
	VelocityAndAcceleration<NDIM> velocityAndAcceleration(potentialField,
			charge, initialNode);
	foundTet = velocityAndAcceleration.foundTet;

//	// For second order leap-frog, offset position from velocity in time
//	currentPosition -= currentVelocity*dt/2.;
	for (double t=0.; t<tMax; t+=dt) {
		nSteps++;
		assert(!isnan(currentPosition.norm()));
		vect3d previousPosition = currentPosition;
		vect3d previousVelocity = currentVelocity;
		entHandle previousElement = currentElement;
//		currentPosition += dt*currentVelocity;
//		vect3d currentAcceleration(0.,0.,0.);
//		// TODO: replace this hack
//		if (t==0.) {
//			currentElement = mesh_ptr->findTet(previousPosition,
//					currentPosition, initialNode, &foundTet, false);
//			// TODO: figure out why sometimes throw here (grazing orbits?)
////			if (!foundTet)
////				throw;
//			int dimension=mesh_ptr->getEntityDimension(currentElement);
//			if (dimension!=iBase_REGION)
//				foundTet = false;
//		}
////		if (!foundTet)
////			break;
//		if (foundTet) {
			try {
				positionAndVelocity[0] = currentPosition;
				positionAndVelocity[1] = currentVelocity;
				timestepper.do_step(boost::ref(velocityAndAcceleration), positionAndVelocity, t, dt);
//				// TODO: figure out why inout arg leads to problems
//				//		 (Eigen not liking certain optimization tricks in odeint?.
				//		 Actually, just missing initialization)
//				timestepper.do_step(boost::ref(velocityAndAcceleration),
//						positionAndVelocity, t, positionAndVelocityOut, dt);
//				positionAndVelocity[0] = positionAndVelocityOut[0];
//				positionAndVelocity[1] = positionAndVelocityOut[1];
//				timestepper.do_step(velocityAndAcceleration, positionAndVelocity, t, dt);
				currentPosition = positionAndVelocity[0];
				currentVelocity = positionAndVelocity[1];
				currentElement = velocityAndAcceleration.currentElement;
				foundTet = mesh_ptr->checkIfInTet(currentPosition, currentElement);
				if (!foundTet) {
					currentElement = mesh_ptr->findTet(previousPosition, currentPosition,
							currentElement, &foundTet);
				}
//				// TODO: set order through input parameter?
//				int interpolationOrder = INTERPOLATIONORDER;
////				currentAcceleration = charge*
////						electricField.getField(currentPosition, &currentElement);
//				potential = potentialField.getField(currentPosition, &currentElement,
//						interpolationOrder);
//				if (isnan(potential))
//					potential = potentialField.getField(currentPosition,
//							&currentElement, 1);
//				// TODO: hard-coding dimension here...
//				for (int i=0; i<3; i++) {
//					vect3d perturbedPosition = currentPosition +
//							vect3d::Unit(i)*DELTA_LENGTH;
//					double perturbedPotential = potentialField.getField(
//							perturbedPosition, &currentElement, interpolationOrder);
//					// TODO: track down why perturbedPotential sometimes is NaN
////					cout << "calculation of error term succeeded." <<
////							i << " " << currentElement << endl;
//					if (isnan(perturbedPotential)) {
//						perturbedPotential = potentialField.getField(
//								perturbedPosition, &currentElement, 1);
//						cout << "calculation of error term failed." <<
//								i << " " << currentElement << endl;
//					}
////					if (isnan(potential) || isnan(perturbedPotential)) {
////						cout << potential << " " << perturbedPotential <<
////						" " << currentPosition.transpose() << endl;
////						throw;
////					}
//					currentAcceleration[i] =
//							-charge*(perturbedPotential-potential)/DELTA_LENGTH;
//				}
			} catch (int signal) {
				switch (signal) {
				case OUTSIDE_DOMAIN:
					foundTet = false;
					// TODO: not guaranteed that outside domain (since different than stepper)
					// TODO: Revisit this in magnetized case
					currentPosition += dt*currentVelocity;
//					currentPosition = positionAndVelocity[0];
//					currentVelocity = positionAndVelocity[1];
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
//		}

		if (foundTet==false) {
			// TODO: if near vertex could leave domain through non-boundary face
			//       by crossing sliver of other tet in one time-step
			entHandle faceCrossed = mesh_ptr->findFaceCrossed(
					previousElement, previousPosition, currentPosition);
			// TODO: grazing orbits don't enter domain in first time-step
			int faceType;
			vect3d normalVector(0.,0.,0.), normalVelocity(0.,0.,0.);
			if (faceCrossed!=NULL) {
				faceType = faceTypeField.getField(faceCrossed);
				normalVector = mesh_ptr->getNormalVector(faceCrossed,
								previousPosition);
				normalVelocity =
						currentVelocity.dot(normalVector)*normalVector;
				finalPotential = 0.;
				vector<entHandle> vertices =
						mesh_ptr->getVertices(faceCrossed);
				for (int i=0; i<vertices.size(); i++) {
					// TODO: should use point where left domain here
					finalPotential += 1./3.*potentialField.getField(vertices[i]);
				}
				// TODO: shouldn't hard-code boundary code
			} else {
				faceType = 0;
			}
//			if (faceType==0) {
//				cout << previousElement << " " <<
//						faceCrossed << " " << currentPosition.norm() <<
//						" " << previousPosition.norm() << endl;
//			}
			finalFaceType = faceType;
			if (faceType==4 && 0.5*pow(normalVelocity.norm(),2.)<charge*phiSurface) {
				currentElement = previousElement;
				foundTet = true;
				// TODO: resetting position isn't quite right
				currentPosition = previousPosition;
				// TODO: revisit this in magnetized case
				// TODO: calc this from previousVelocity?
				currentVelocity -= 2.*normalVelocity;
			} else {
//				if (faceType==0)
//					cout << "last face crossed was an interior one" << endl;
				endLoop=true;
			}
		}

//		if (nSteps==1 && !foundTet) {
//			cout << currentPosition.transpose() << " " <<
//					currentVelocity.transpose() << " " <<
//					currentPosition.dot(currentVelocity) << endl;
//		}

//		double eFieldR = currentAcceleration.dot(currentPosition)/
//				currentPosition.norm();
////		 // TODO: comment out this
////		currentAcceleration = -charge*currentPosition/
////				pow(currentPosition.norm(),3.);
//		currentVelocity += dt*currentAcceleration;
//		vect3d velocityAtPosition = currentVelocity - 1./2.*dt*currentAcceleration;
//		double energy = 1./2.*pow(velocityAtPosition.norm(),2.) + charge*potential;
		if (foundTet) {
			potential = potentialField.getField(currentPosition,
					&velocityAndAcceleration.currentElement,
					velocityAndAcceleration.interpolationOrder);
		} else {
			potential = 0.;
		}
		double energy = 1./2.*pow(currentVelocity.norm(),2.) + charge*potential;
		if (outFile) {
//			fprintf(outFile, "%f %f %f %p\n", currentPosition[0], currentPosition[1],
//					currentPosition[2], (void*)currentElement);
			fprintf(outFile, "%f %f %f %f\n", currentPosition[0], currentPosition[1],
					currentPosition[2], energy);
//					currentPosition[2], eFieldR);
		}
		if (endLoop)
			break;
	}
//	cout << "Final radius=" << currentPosition.norm() << " nSteps=" <<
//			nSteps << "faceType =" << finalFaceType << endl;
	finalPosition = currentPosition;
	// TODO: correct for time-step offset?
	finalVelocity = currentVelocity;
}
