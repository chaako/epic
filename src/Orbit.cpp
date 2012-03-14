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
		Eigen::Vector3d inputVelocity) {
	mesh_ptr = inputMesh_ptr;
	initialNode = inputNode;
	// TODO: perhaps find a random adjacent tet rather than give node
	currentElement = inputNode;
	initialPosition = inputMesh_ptr->getCoordinates(inputNode);
	initialVelocity = inputVelocity;
}

Orbit::~Orbit() {
	// TODO Auto-generated destructor stub
}

void Orbit::integrate(ElectricField& electricField,
		PotentialField& potentialField, FILE *outFile) {
	Eigen::Vector3d eField = electricField.getField(initialNode);
	int nVertices=4;
	std::vector<Eigen::Vector3d> eFields(nVertices), vertexVectors(nVertices);
	std::vector<iBase_EntityHandle> vertices(nVertices);
	// TODO: need some clever way to set tMax and/or detect trapped orbits
	double dt=std::min(0.01,0.01/initialVelocity.norm()), tMax=100;
	Eigen::Vector3d currentPosition = initialPosition;
	Eigen::Vector3d currentVelocity = initialVelocity;
	// Don't integrate orbit if doesn't have enough energy to escape potential
	// TODO: this should be refined
//	if (0.5*pow(initialVelocity.norm(),2.)+0.22 <
	if ((0.5*pow(initialVelocity.norm(),2.) +
			potentialField.getField(initialNode)) < 0.) {
		negativeEnergy = true;
		tMax = 0.;
	} else {
		negativeEnergy = false;
	}
	// For second order leap-frog, offset position from velocity in time
	currentPosition -= currentVelocity*dt/2.;
	bool inNewTet = true;
	int nSteps=0, nNewTet=0;
//	std::cout << "Initial radius=" << currentPosition.norm() << std::endl;
	bool isTet=false;
	bool firstStep=true;
	for (double t=0; t<tMax; t+=dt) {
		nSteps++;
		currentPosition += dt*currentVelocity;
//		clock_t startClock = clock(); // timing
		if (!firstStep)
			inNewTet = !checkIfInTet(currentPosition, vertexVectors);
		firstStep=false;
//		clock_t endClock = clock(); // timing
//		extern_checkIfInNewTet += endClock-startClock; // timing
		if (inNewTet) {
			nNewTet++;
			bool foundTet=false;
//			if  (vertexVectors.size()==nVertices)
//				isTet = true;
			clock_t startClock = clock(); // timing
			currentElement = mesh_ptr->findTet(currentPosition,
					currentElement, &foundTet, isTet);
			clock_t endClock = clock(); // timing
			extern_findTet += endClock-startClock; // timing
			// TODO: should handle failure to find tet in some way
			if (foundTet==false)
				break;
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
		std::vector<double> vertexWeights = getVertexWeights(currentPosition,
				vertexVectors);
		Eigen::Vector3d currentAcceleration(0.,0.,0.);
		assert(eFields.size()==vertexWeights.size());
		for (int i=0; i<vertexWeights.size(); i++) {
			currentAcceleration += eFields[i]*vertexWeights[i];
		}
//		currentAcceleration = -currentPosition/pow(currentPosition.norm(),3.);
		double eFieldR = currentAcceleration.dot(currentPosition)/
				currentPosition.norm();
		currentVelocity += dt*currentAcceleration;
		if (outFile) {
			fprintf(outFile, "%f %f %f %f\n", currentPosition[0], currentPosition[1],
					currentPosition[2], eFieldR);
		}
	}
//	std::cout << "Final radius=" << currentPosition.norm() << std::endl;
	finalPosition = currentPosition;
	// TODO: correct for time-step offset?
	finalVelocity = currentVelocity;
}
